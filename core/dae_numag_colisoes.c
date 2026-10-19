/* dae_numag_colisoes.c — colisões como defasagem, por trajetórias. O modelo,
 * os fatores e o desdobramento estão em dae_numag.h, seção "colisões e
 * coerência". */
#include "dae_numag.h"

#include "dae_rng.h"
#include "dae_traj.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

double dae_numag_salto_por_km(double gamma_colisao_eV)
{
  return 0.5 * gamma_colisao_eV * DAE_NU_T_POR_KM;
}

/* gamma (1/km) em s: linear entre os nós da trajetória, segura nas bordas. */
static double col_gamma(const dae_numag_trajetoria *T, const dae_numag_colisao *c,
                        int32_t i, double s)
{
  double t;
  if (!c->gamma_km) return c->gamma0_km;
  t = (s - T->s_km[i]) / (T->s_km[i + 1] - T->s_km[i]);
  if (t < 0.0) t = 0.0;
  if (t > 1.0) t = 1.0;
  return c->gamma_km[i] + t * (c->gamma_km[i + 1] - c->gamma_km[i]);
}

static int32_t col_intervalo(const dae_numag_trajetoria *T, double s)
{
  int32_t lo = 0, hi = T->n - 1;
  if (s <= T->s_km[0]) return 0;
  if (s >= T->s_km[T->n - 1]) return T->n - 2;
  while (hi - lo > 1) {
    const int32_t m = lo + (hi - lo) / 2;
    if (T->s_km[m] <= s) lo = m; else hi = m;
  }
  return lo;
}

/* Anda de s até no máximo `lim` consumindo a taxa acumulada `*resto`. Devolve
   o ponto onde ela zera (salto) ou `lim` (sem salto, `*resto` diminuído).
   Exato para gamma linear em cada intervalo da tabela. */
static double col_avanca(const dae_numag_trajetoria *T, const dae_numag_colisao *c,
                         double s, double lim, double *resto, int *saltou)
{
  *saltou = 0;
  if (!c->gamma_km) {
    const double g = c->gamma0_km;
    if (g > 0.0 && s + *resto / g <= lim) {
      s += *resto / g; *resto = 0.0; *saltou = 1;
      return s;
    }
    *resto -= g * (lim - s);
    return lim;
  }
  while (s < lim) {
    const int32_t i = col_intervalo(T, s);
    const double fim = (T->s_km[i + 1] < lim && T->s_km[i + 1] > s) ? T->s_km[i + 1] : lim;
    const double ga = col_gamma(T, c, i, s), gb = col_gamma(T, c, i, fim);
    const double area = 0.5 * (ga + gb) * (fim - s);
    if (area >= *resto && *resto >= 0.0) {
      const double k = (gb - ga) / (fim - s);
      double x;
      if (fabs(k) * (fim - s) < 1e-12 * (fabs(ga) + fabs(gb)) || k == 0.0)
        x = (ga > 0.0) ? *resto / ga : (fim - s);
      else
        x = (-ga + sqrt(ga * ga + 2.0 * k * *resto)) / k;   /* k/2 x^2 + ga x = resto */
      if (x < 0.0) x = 0.0;
      if (x > fim - s) x = fim - s;
      *resto = 0.0; *saltou = 1;
      return s + x;
    }
    *resto -= area;
    s = fim;
  }
  return lim;
}

static void col_fisico(const double *re, const double *im, double *pr, double *pi)
{
  dae_nu_projeta(re, im, 2, 1, pr, pi);
}

static void col_mergulho(const double *pr, const double *pi, double *re, double *im)
{
  const double c = 1.0 / sqrt(2.0);
  int a;
  for (a = 0; a < 2; ++a) {
    re[a] = pr[a] * c; im[a] = pi[a] * c;
    re[2 + a] = pi[a] * c; im[2 + a] = -pr[a] * c;
  }
}

/* Uma trajetória de s0 a s1, registrando nas `nt` amostras. */
static dae_status col_trajetoria(const dae_numag_params *P, const dae_numag_trajetoria *T,
                                 double s0, double s1, const dae_numag_colisao *c,
                                 double tol, dae_rng *g, const double *pr0, const double *pi0,
                                 int32_t nt, double *amostra_pr, double *amostra_cre,
                                 double *amostra_cim, int32_t *n_saltos)
{
  double re[4], im[4], pr[2], pi[2], s = s0, resto;
  int32_t k;
  dae_status st;

  col_mergulho(pr0, pi0, re, im);
  *n_saltos = 0;
  resto = -log(1.0 - dae_rng_uniform(g));          /* u em (0, 1] */
  for (k = 0; k < nt; ++k) {
    const double alvo = s0 + (s1 - s0) * (double)(k + 1) / (double)nt;
    while (s < alvo) {
      int saltou;
      const double prox = col_avanca(T, c, s, alvo, &resto, &saltou);
      if (prox > s) {
        st = dae_numag_propaga_tabela(P, T, s, prox, tol, re, im, NULL);
        if (st != DAE_OK) return st;
      }
      s = prox;
      if (saltou) {
        /* Salto: P_j com probabilidade |psi_j|^2, fase de psi_j mantida. */
        double nl, nr;
        col_fisico(re, im, pr, pi);
        nl = pr[0] * pr[0] + pi[0] * pi[0];
        nr = pr[1] * pr[1] + pi[1] * pi[1];
        if (dae_rng_uniform(g) * (nl + nr) < nl) {
          const double m = sqrt(nl);
          pr[0] /= m; pi[0] /= m; pr[1] = 0.0; pi[1] = 0.0;
        } else {
          const double m = sqrt(nr);
          pr[1] /= m; pi[1] /= m; pr[0] = 0.0; pi[0] = 0.0;
        }
        col_mergulho(pr, pi, re, im);
        ++*n_saltos;
        resto = -log(1.0 - dae_rng_uniform(g));
      }
    }
    col_fisico(re, im, pr, pi);
    amostra_pr[k] = pr[1] * pr[1] + pi[1] * pi[1];
    /* psi_L psi_R^* */
    amostra_cre[k] = pr[0] * pr[1] + pi[0] * pi[1];
    amostra_cim[k] = pi[0] * pr[1] - pr[0] * pi[1];
  }
  return DAE_OK;
}

dae_status dae_numag_ensemble(const dae_numag_params *P,
                              const dae_numag_trajetoria *T,
                              double s0, double s1,
                              const dae_numag_colisao *col, double tol,
                              int32_t n_traj, uint64_t semente_base,
                              const double *pr0, const double *pi0,
                              int32_t nt, double *p_r, double *p_r_var,
                              double *coh_re, double *coh_im,
                              double *saltos_medio)
{
  double *a_pr = NULL, *a_cre = NULL, *a_cim = NULL;
  double *soma = NULL, *soma2 = NULL, *scre = NULL, *scim = NULL;
  double saltos = 0.0;
  int32_t idx, k;
  dae_status st = DAE_OK;

  if (!P || !T || !col || !pr0 || !pi0 || n_traj < 1 || nt < 1 || T->n < 2 ||
      !(s1 > s0) || !(tol > 0.0) || !(col->gamma0_km >= 0.0))
    return DAE_ERR_PARAM;
  a_pr = (double *)malloc((size_t)nt * sizeof(double));
  a_cre = (double *)malloc((size_t)nt * sizeof(double));
  a_cim = (double *)malloc((size_t)nt * sizeof(double));
  soma = (double *)calloc((size_t)nt, sizeof(double));
  soma2 = (double *)calloc((size_t)nt, sizeof(double));
  scre = (double *)calloc((size_t)nt, sizeof(double));
  scim = (double *)calloc((size_t)nt, sizeof(double));
  if (!a_pr || !a_cre || !a_cim || !soma || !soma2 || !scre || !scim) {
    st = DAE_ERR_ALLOC;
    goto fim;
  }

  /* Uma trajetória por vez, na ordem do índice: as somas não dependem de
     quem terminou primeiro (CONVENTIONS.md 11.6). */
  for (idx = 0; idx < n_traj && st == DAE_OK; ++idx) {
    dae_rng g;
    int32_t ns = 0;
    dae_rng_seed(&g, dae_traj_semente(semente_base, idx));
    st = col_trajetoria(P, T, s0, s1, col, tol, &g, pr0, pi0, nt, a_pr, a_cre, a_cim, &ns);
    if (st != DAE_OK) break;
    /* Welford, na ordem do índice: média e soma dos quadrados dos DESVIOS.
       A forma soma2 - n m^2 perde dígitos por cancelamento (medido: variância
       3,6e-17 para 50 valores iguais); Welford dá zero exato nesse caso. */
    for (k = 0; k < nt; ++k) {
      const double d = a_pr[k] - soma[k];
      soma[k] += d / (double)(idx + 1);
      soma2[k] += d * (a_pr[k] - soma[k]);
      scre[k] += a_cre[k];
      scim[k] += a_cim[k];
    }
    saltos += (double)ns;
  }
  if (st == DAE_OK) {
    const double n = (double)n_traj;
    for (k = 0; k < nt; ++k) {
      if (p_r) p_r[k] = soma[k];                  /* soma[] já é a média */
      if (p_r_var) p_r_var[k] = (n > 1.0) ? soma2[k] / (n - 1.0) : 0.0;
      if (coh_re) coh_re[k] = scre[k] / n;
      if (coh_im) coh_im[k] = scim[k] / n;
    }
    if (saltos_medio) *saltos_medio = saltos / n;
  }
fim:
  free(a_pr); free(a_cre); free(a_cim);
  free(soma); free(soma2); free(scre); free(scim);
  return st;
}
