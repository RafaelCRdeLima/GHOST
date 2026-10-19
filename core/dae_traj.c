#include "dae_traj.h"
#include "dae_rng.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static dae_status dae_traj_ensemble_i(const dae_traj_cfg *cfg, uint64_t semente_base,
                                      dae_cheb *W, const dae_csr *H,
                                      dae_rho_acc *A, dae_traj_sink sink, void *user,
                                      dae_rho_acc *Aq, int32_t corte);

int32_t dae_traj_amostras(const dae_traj_cfg *cfg)
{
  if (!cfg || cfg->nt <= 0 || cfg->rho_stride <= 0) return 0;
  /* Inclui t = 0 e depois um a cada rho_stride passos. */
  return 1 + cfg->nt / cfg->rho_stride;
}

/* splitmix64 sobre (base, idx): fluxo proprio por trajetoria, derivavel do
 * indice e independente da ordem em que as threads executam. Um gerador
 * compartilhado faria o resultado depender do escalonador. */
uint64_t dae_traj_semente(uint64_t base, int32_t idx)
{
  uint64_t z = base + (uint64_t)0x9E3779B97F4A7C15ULL * (uint64_t)(idx + 1);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static dae_status acc_init(dae_rho_acc *A, int32_t n, int32_t n_amostras, int com_m2)
{
  size_t tot;
  if (!A || n <= 0 || n_amostras <= 0) return DAE_ERR_PARAM;
  A->n = n; A->n_amostras = n_amostras; A->proxima = 0; A->somadas = 0;
  A->n_traj_alvo = 0; A->m2 = NULL;
  tot = (size_t)n_amostras * (size_t)n * (size_t)n;
  A->re = (double *)calloc(tot, sizeof(double));
  A->im = (double *)calloc(tot, sizeof(double));
  if (!A->re || !A->im) { dae_rho_acc_free(A); return DAE_ERR_ALLOC; }
  if (com_m2) {
    A->m2 = (double *)calloc(tot, sizeof(double));
    if (!A->m2) { dae_rho_acc_free(A); return DAE_ERR_ALLOC; }
  }
  return DAE_OK;
}

dae_status dae_rho_acc_init(dae_rho_acc *A, int32_t n, int32_t n_amostras)
{ return acc_init(A, n, n_amostras, 0); }

dae_status dae_rho_acc_init_wk(dae_rho_acc *A, int32_t n, int32_t n_amostras)
{ return acc_init(A, n, n_amostras, 1); }

double dae_rho_acc_mod2_sem_vies(const dae_rho_acc *A, int32_t amostra,
                                 int32_t i, int32_t j, int32_t n_traj)
{
  size_t k;
  double q, u;
  if (!A || !A->m2 || n_traj < 2) return 0.0;
  k = (size_t)amostra * (size_t)A->n * (size_t)A->n +
      (size_t)i * (size_t)A->n + (size_t)j;
  q = A->re[k] * A->re[k] + A->im[k] * A->im[k];
  u = ((double)n_traj * q - A->m2[k]) / (double)(n_traj - 1);
  return u > 0.0 ? u : 0.0;
}

void dae_rho_acc_free(dae_rho_acc *A)
{
  if (!A) return;
  free(A->re); free(A->im); free(A->m2);
  A->re = NULL; A->im = NULL; A->m2 = NULL; A->n = 0; A->n_amostras = 0;
}

dae_status dae_rho_acc_somar(dae_rho_acc *A, int32_t idx,
                             const double *re, const double *im)
{
  int32_t s, i, j;
  if (!A || !A->re || !re || !im) return DAE_ERR_PARAM;
  /* A ORDEM E CONTRATO. Fora de ordem e erro, nao outro numero: soma de ponto
     flutuante nao e associativa, e a igualdade bit a bit entre os dois modos
     depende de as contribuicoes entrarem exatamente nesta sequencia. */
  if (idx != A->proxima) return DAE_ERR_PARAM;

  for (s = 0; s < A->n_amostras; ++s) {
    const double *pr = re + (size_t)s * (size_t)A->n;
    const double *pi = im + (size_t)s * (size_t)A->n;
    double *ar = A->re + (size_t)s * (size_t)A->n * (size_t)A->n;
    double *ai = A->im + (size_t)s * (size_t)A->n * (size_t)A->n;
    double *am = A->m2 ? A->m2 + (size_t)s * (size_t)A->n * (size_t)A->n : NULL;
    for (i = 0; i < A->n; ++i) {
      const double xr = pr[i], xi = pi[i];
      const double qi = xr * xr + xi * xi;
      for (j = 0; j < A->n; ++j) {
        /* psi_i psi_j^* */
        ar[(size_t)i * (size_t)A->n + (size_t)j] += xr * pr[j] + xi * pi[j];
        ai[(size_t)i * (size_t)A->n + (size_t)j] += xi * pr[j] - xr * pi[j];
        /* |X_ij|^2 = |psi_i|^2 |psi_j|^2, e e REAL */
        if (am) am[(size_t)i * (size_t)A->n + (size_t)j] +=
          qi * (pr[j] * pr[j] + pi[j] * pi[j]);
      }
    }
  }
  A->proxima = idx + 1;
  A->somadas += 1;
  return DAE_OK;
}

void dae_rho_acc_finalizar(dae_rho_acc *A)
{
  size_t tot, k;
  double inv;
  if (!A || !A->re || A->somadas <= 0) return;
  tot = (size_t)A->n_amostras * (size_t)A->n * (size_t)A->n;
  inv = 1.0 / (double)A->somadas;
  for (k = 0; k < tot; ++k) { A->re[k] *= inv; A->im[k] *= inv; }
  if (A->m2) for (k = 0; k < tot; ++k) A->m2[k] *= inv;
}

/* Colapso: escolhe j com probabilidade |psi_j|^2 e projeta. A busca e linear e
 * acumulativa — a MESMA travessia em qualquer plataforma, porque uma busca
 * binaria sobre soma parcial daria outro j para o mesmo sorteio se a soma
 * parcial fosse calculada em outra ordem. */
/* Tempo do proximo salto a partir de `t0`, dado um uniforme `u`.
 *
 * Com gamma constante e a forma fechada de sempre: -log(1-u)/gamma. Com
 * gamma(t) e a solucao de   integral_{t0}^{t} gamma(t') dt' = -log(1-u),
 * que com gamma constante POR PASSO se resolve andando passo a passo e
 * consumindo o "orcamento" -log(1-u). A taxa total continua independente do
 * estado, entao nao ha rejeicao — e o motivo de gamma(t) ser barato e
 * gamma_j (por sitio) nao ser. */
void dae_traj_cfg_default(dae_traj_cfg *cfg)
{
  if (!cfg) return;
  memset(cfg, 0, sizeof(*cfg));
  cfg->rho_stride = 1;
  cfg->saida = DAE_SAIDA_INDEFINIDA;
}

static double proximo_salto(const dae_traj_cfg *cfg, double t0, double u)
{
  double orc, t, gama;
  int32_t k;

  if (!cfg->gamma_tab) {
    if (!(cfg->gamma_deph > 0.0)) return (double)INFINITY;
    return t0 - log(1.0 - u) / cfg->gamma_deph;
  }
  orc = -log(1.0 - u);
  t = t0;
  k = (int32_t)(t / cfg->dt);
  if (k < 0) k = 0;
  while (k < cfg->nt) {
    const double t_fim = (double)(k + 1) * cfg->dt;
    gama = cfg->gamma_tab[k];
    if (gama > 0.0) {
      const double cabe = gama * (t_fim - t);
      if (cabe >= orc) return t + orc / gama;
      orc -= cabe;
    }
    t = t_fim;
    ++k;
  }
  return (double)INFINITY;            /* a taxa acabou antes do orcamento */
}

static void colapsar(double *re, double *im, int32_t n, double u, int32_t nf_par)
{
  double acc = 0.0;
  int32_t j, escolhido;

  if (nf_par <= 0) {                       /* caso generico: colapso em sitio */
    escolhido = n - 1;
    for (j = 0; j < n; ++j) {
      acc += re[j] * re[j] + im[j] * im[j];
      if (u < acc) { escolhido = j; break; }
    }
    for (j = 0; j < n; ++j) { re[j] = 0.0; im[j] = 0.0; }
    re[escolhido] = 1.0;
    return;
  }

  /* MERGULHADO: o bloco j sao as componentes j e nf_par+j juntas, e a soma dos
     dois modulos ao quadrado ja e |psi_j|^2 — a probabilidade de sabor, sem
     correcao nenhuma. A ordem da soma acumulada e a MESMA do caso generico,
     por bloco crescente, para o determinismo bit a bit nao depender disto. */
  escolhido = nf_par - 1;
  for (j = 0; j < nf_par; ++j) {
    acc += re[j] * re[j] + im[j] * im[j]
         + re[nf_par + j] * re[nf_par + j] + im[nf_par + j] * im[nf_par + j];
    if (u < acc) { escolhido = j; break; }
  }
  for (j = 0; j < n; ++j) { re[j] = 0.0; im[j] = 0.0; }
  /* volta para dentro do subespaco: (e_j, -i e_j)/sqrt(2) */
  re[escolhido] = 0.70710678118654752440;
  im[nf_par + escolhido] = -0.70710678118654752440;
}

dae_status dae_traj_ensemble_dupla(const dae_traj_cfg *cfg, uint64_t semente_base,
                                   dae_cheb *W, const dae_csr *H,
                                   dae_rho_acc *A, dae_rho_acc *Aq, int32_t corte)
{
  return dae_traj_ensemble_i(cfg, semente_base, W, H, A, NULL, NULL, Aq, corte);
}

dae_status dae_traj_ensemble(const dae_traj_cfg *cfg, uint64_t semente_base,
                             dae_cheb *W, const dae_csr *H,
                             dae_rho_acc *A, dae_traj_sink sink, void *user)
{
  return dae_traj_ensemble_i(cfg, semente_base, W, H, A, sink, user, NULL, 0);
}

static dae_status dae_traj_ensemble_i(const dae_traj_cfg *cfg, uint64_t semente_base,
                                      dae_cheb *W, const dae_csr *H,
                                      dae_rho_acc *A, dae_traj_sink sink, void *user,
                                      dae_rho_acc *Aq, int32_t corte)
{
  const int32_t n = H ? H->n : 0;
  int32_t na, idx, s;
  double *psre = NULL, *psim = NULL, *amre = NULL, *amim = NULL;
  dae_status st = DAE_OK;

  if (!cfg || !W || !H || n <= 0) return DAE_ERR_PARAM;
  if (cfg->saida == DAE_SAIDA_INDEFINIDA) return DAE_ERR_PARAM;
  if (cfg->n_traj <= 0 || cfg->gamma_deph < 0.0) return DAE_ERR_PARAM;
  if (cfg->sitio_inicial < 0 || cfg->sitio_inicial >= n) return DAE_ERR_PARAM;
  if (cfg->saida == DAE_SAIDA_ACUMULAR_RHO && !A) return DAE_ERR_PARAM;
  if (cfg->saida == DAE_SAIDA_ARQUIVAR_PSI && !sink) return DAE_ERR_PARAM;

  na = dae_traj_amostras(cfg);
  if (na <= 0) return DAE_ERR_PARAM;

  psre = (double *)calloc((size_t)n, sizeof(double));
  psim = (double *)calloc((size_t)n, sizeof(double));
  amre = (double *)calloc((size_t)na * (size_t)n, sizeof(double));
  amim = (double *)calloc((size_t)na * (size_t)n, sizeof(double));
  if (!psre || !psim || !amre || !amim) { st = DAE_ERR_ALLOC; goto fim; }

  for (idx = 0; idx < cfg->n_traj; ++idx) {
    dae_rng g;
    double t = 0.0, t_salto;
    int32_t passo;
    dae_rng_seed(&g, dae_traj_semente(semente_base, idx));

    if (cfg->init_re && cfg->init_im) {
      memcpy(psre, cfg->init_re, (size_t)n * sizeof(double));
      memcpy(psim, cfg->init_im, (size_t)n * sizeof(double));
    } else {
      memset(psre, 0, (size_t)n * sizeof(double));
      memset(psim, 0, (size_t)n * sizeof(double));
      psre[cfg->sitio_inicial] = 1.0;
    }

    /* Tempo do primeiro salto. gamma = 0 nunca salta: -log(u)/0 = +inf, e o
       teste de igualdade entre modos roda tambem nesse caso. */
    t_salto = proximo_salto(cfg, 0.0, dae_rng_uniform(&g));

    s = 0;
    memcpy(amre, psre, (size_t)n * sizeof(double));
    memcpy(amim, psim, (size_t)n * sizeof(double));
    s = 1;

    for (passo = 1; passo <= cfg->nt; ++passo) {
      const double t_fim = (double)passo * cfg->dt;
      while (t_salto < t_fim) {
        dae_cheb_info info;
        const double sub = t_salto - t;
        if (sub > 0.0) {
          st = dae_cheb_step(W, H, sub, psre, psim, &info);
          if (st != DAE_OK) goto fim;
        }
        colapsar(psre, psim, n, dae_rng_uniform(&g), cfg->nf_par);
        t = t_salto;
        t_salto = proximo_salto(cfg, t, dae_rng_uniform(&g));
      }
      if (t_fim > t) {
        dae_cheb_info info;
        st = dae_cheb_step(W, H, t_fim - t, psre, psim, &info);
        if (st != DAE_OK) goto fim;
        t = t_fim;
      }
      if (passo % cfg->rho_stride == 0 && s < na) {
        memcpy(amre + (size_t)s * (size_t)n, psre, (size_t)n * sizeof(double));
        memcpy(amim + (size_t)s * (size_t)n, psim, (size_t)n * sizeof(double));
        ++s;
      }
    }

    if (cfg->saida == DAE_SAIDA_ACUMULAR_RHO) {
      st = dae_rho_acc_somar(A, idx, amre, amim);
      if (st == DAE_OK && Aq && idx < corte) st = dae_rho_acc_somar(Aq, idx, amre, amim);
    } else {
      st = sink(user, idx, amre, amim, na, n);
    }
    if (st != DAE_OK) goto fim;
  }

fim:
  free(psre); free(psim); free(amre); free(amim);
  return st;
}
