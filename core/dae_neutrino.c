#include "dae_neutrino.h"

#include "dae_cheb.h"
#include "dae_rng.h"
#include "dae_ham.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void dae_nu_params_default(dae_nu_params *P)
{
  if (!P) return;
  memset(P, 0, sizeof(*P));
  P->n_flavors = 3;
  /* NuFIT-ish, ordenamento normal. São DEFAULTS, não medidas: quem publica
     número põe os dele no spec, que é onde o valor fica registrado. */
  P->theta12 = 33.4 * 3.14159265358979323846 / 180.0;
  P->theta13 =  8.6 * 3.14159265358979323846 / 180.0;
  P->theta23 = 49.0 * 3.14159265358979323846 / 180.0;
  P->delta_cp = 0.0;
  P->dm21_sq = 7.42e-5;
  P->dm31_sq = 2.51e-3;
  P->energy_eV = 1.0e6;         /* 1 MeV */
  P->antineutrino = 0;
}

dae_status dae_nu_pmns(const dae_nu_params *P, double *ure, double *uim)
{
  int32_t nf, i;
  double c12, s12, c13, s13, c23, s23, cd, sd, dcp;

  if (!P || !ure || !uim) return DAE_ERR_PARAM;
  nf = P->n_flavors;
  if (nf < 2 || nf > DAE_NU_MAX_SABORES) return DAE_ERR_PARAM;

  for (i = 0; i < nf * nf; ++i) { ure[i] = 0.0; uim[i] = 0.0; }

  if (nf == 2) {
    /* Dois sabores: uma rotação real. Não há fase CP em 2x2 — ela seria
       removível por redefinição de fase, e é por isso que delta_cp só ganha
       significado com três gerações. */
    c12 = cos(P->theta12); s12 = sin(P->theta12);
    ure[0] =  c12; ure[1] = s12;
    ure[2] = -s12; ure[3] = c12;
    return DAE_OK;
  }

  /* Antineutrino: delta -> -delta. O potencial de matéria troca de sinal em
     dae_nu_hamiltoniano; os dois juntos é que dão a assimetria CP observável. */
  dcp = P->antineutrino ? -P->delta_cp : P->delta_cp;

  c12 = cos(P->theta12); s12 = sin(P->theta12);
  c13 = cos(P->theta13); s13 = sin(P->theta13);
  c23 = cos(P->theta23); s23 = sin(P->theta23);
  cd  = cos(dcp);        sd  = sin(dcp);         /* e^{+i delta} = cd + i sd */

  /* Convenção do PDG: U_e3 = s13 e^{-i delta}; os demais trazem e^{+i delta}. */
  ure[0 * 3 + 0] = c12 * c13;
  ure[0 * 3 + 1] = s12 * c13;
  ure[0 * 3 + 2] =  s13 * cd;  uim[0 * 3 + 2] = -s13 * sd;

  ure[1 * 3 + 0] = -s12 * c23 - c12 * s23 * s13 * cd;
  uim[1 * 3 + 0] =            - c12 * s23 * s13 * sd;
  ure[1 * 3 + 1] =  c12 * c23 - s12 * s23 * s13 * cd;
  uim[1 * 3 + 1] =            - s12 * s23 * s13 * sd;
  ure[1 * 3 + 2] =  s23 * c13;

  ure[2 * 3 + 0] =  s12 * s23 - c12 * c23 * s13 * cd;
  uim[2 * 3 + 0] =            - c12 * c23 * s13 * sd;
  ure[2 * 3 + 1] = -c12 * s23 - s12 * c23 * s13 * cd;
  uim[2 * 3 + 1] =            - s12 * c23 * s13 * sd;
  ure[2 * 3 + 2] =  c23 * c13;

  return DAE_OK;
}

dae_status dae_nu_hamiltoniano(const dae_nu_params *P, double rho_g_cm3, double ye,
                               double *hre, double *him)
{
  double ure[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  double uim[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  double lam[DAE_NU_MAX_SABORES];
  double v, tr;
  int32_t nf, a, b, i;
  dae_status st;

  if (!P || !hre || !him) return DAE_ERR_PARAM;
  nf = P->n_flavors;
  if (nf < 2 || nf > DAE_NU_MAX_SABORES) return DAE_ERR_PARAM;
  if (!(P->energy_eV > 0.0)) return DAE_ERR_PARAM;
  if (rho_g_cm3 < 0.0 || ye < 0.0) return DAE_ERR_PARAM;
  /* NaN em qualquer entrada envenena tudo em silêncio; recusa aqui. */
  if (P->energy_eV != P->energy_eV || rho_g_cm3 != rho_g_cm3 || ye != ye)
    return DAE_ERR_PARAM;

  st = dae_nu_pmns(P, ure, uim);
  if (st != DAE_OK) return st;

  /* Autovalores de massa, em eV: Dm^2/2E. O primeiro é zero por convenção —
     só diferenças são observáveis. */
  lam[0] = 0.0;
  lam[1] = P->dm21_sq / (2.0 * P->energy_eV);
  if (nf == 3) lam[2] = P->dm31_sq / (2.0 * P->energy_eV);

  /* H_ab = sum_i U_ai lam_i conj(U_bi) */
  for (a = 0; a < nf; ++a) {
    for (b = 0; b < nf; ++b) {
      double sr = 0.0, si = 0.0;
      for (i = 0; i < nf; ++i) {
        const double ar = ure[a * nf + i], ai = uim[a * nf + i];
        const double br = ure[b * nf + i], bi = -uim[b * nf + i];   /* conj */
        sr += lam[i] * (ar * br - ai * bi);
        si += lam[i] * (ar * bi + ai * br);
      }
      hre[a * nf + b] = sr;
      him[a * nf + b] = si;
    }
  }

  /* Wolfenstein: o potencial de corrente carregada afeta SÓ nu_e, e troca de
     sinal para antineutrino. É o que faz a ressonância MSW acontecer em nu ou
     em anti-nu conforme a hierarquia. */
  v = DAE_NU_V_POR_RHO * ye * rho_g_cm3;
  if (P->antineutrino) v = -v;
  hre[0] += v;

  /* Remove o traço: fase global não é observável, e centrar o espectro reduz
     a ordem de Chebyshev. */
  tr = 0.0;
  for (a = 0; a < nf; ++a) tr += hre[a * nf + a];
  tr /= (double)nf;
  for (a = 0; a < nf; ++a) hre[a * nf + a] -= tr;

  return DAE_OK;
}

int dae_nu_precisa_mergulho(const double *him, int32_t nf)
{
  int32_t i;
  if (!him || nf <= 0) return 0;
  for (i = 0; i < nf * nf; ++i) if (him[i] != 0.0) return 1;
  return 0;
}

dae_status dae_nu_adjacencia(const double *hre, const double *him, int32_t nf,
                             dae_csr *A, int *mergulhado)
{
  int32_t *ei = NULL, *ej = NULL;
  double *w = NULL;
  int32_t n, ne = 0, cap, a, b;
  int merg;
  dae_status st;

  if (!hre || !him || !A || !mergulhado || nf < 2 || nf > DAE_NU_MAX_SABORES)
    return DAE_ERR_PARAM;

  merg = dae_nu_precisa_mergulho(him, nf);
  n = merg ? 2 * nf : nf;
  cap = n * (n + 1) / 2;                 /* triângulo superior com diagonal */

  ei = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  ej = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  w  = (double  *)malloc((size_t)cap * sizeof(double));
  if (!ei || !ej || !w) { free(ei); free(ej); free(w); return DAE_ERR_ALLOC; }

  for (a = 0; a < n; ++a) {
    for (b = a; b < n; ++b) {
      double v;
      if (!merg) {
        v = hre[a * nf + b];
      } else if (a < nf && b < nf) {
        v = hre[a * nf + b];                  /* bloco X  */
      } else if (a < nf) {
        v = -him[a * nf + (b - nf)];          /* bloco -Y */
      } else {
        v = hre[(a - nf) * nf + (b - nf)];    /* bloco X  */
      }
      /* O pipeline monta H = -gamma*A com gamma = 1, então A carrega -H. */
      v = -v;
      if (v != 0.0) { ei[ne] = a; ej[ne] = b; w[ne] = v; ++ne; }
    }
  }

  st = dae_csr_from_edges(A, n, ei, ej, w, ne, NULL);
  free(ei); free(ej); free(w);
  if (st != DAE_OK) return st;

  *mergulhado = merg;
  return DAE_OK;
}

void dae_nu_inicial(int32_t sabor, int32_t nf, int mergulhado,
                    double *re, double *im)
{
  const int32_t n = mergulhado ? 2 * nf : nf;
  int32_t i;
  if (!re || !im || nf <= 0 || sabor < 0 || sabor >= nf) return;
  for (i = 0; i < n; ++i) { re[i] = 0.0; im[i] = 0.0; }
  if (!mergulhado) {
    re[sabor] = 1.0;
  } else {
    /* phi = (psi, -i psi)/sqrt(2), com psi = e_sabor */
    const double c = 1.0 / sqrt(2.0);
    re[sabor] = c;
    im[nf + sabor] = -c;
  }
}

void dae_nu_projeta(const double *re, const double *im, int32_t nf, int mergulhado,
                    double *pre, double *pim)
{
  int32_t a;
  if (!re || !im || !pre || !pim || nf <= 0) return;
  if (!mergulhado) {
    for (a = 0; a < nf; ++a) { pre[a] = re[a]; pim[a] = im[a]; }
    return;
  }
  /* psi_a = (phi_a + i phi_{n+a})/sqrt(2). É a projeção ORTOGONAL no subespaço
     invariante, então também corrige a deriva numérica em vez de só desfazer
     o mergulho. */
  for (a = 0; a < nf; ++a) {
    const double c = 1.0 / sqrt(2.0);
    pre[a] = (re[a] - im[nf + a]) * c;
    pim[a] = (im[a] + re[nf + a]) * c;
  }
}

/* ------------------------------------------------------- perfil de materia -- */

void dae_nu_perfil_default(dae_nu_perfil *M)
{
  if (!M) return;
  memset(M, 0, sizeof(*M));
  M->kind = DAE_NU_PERFIL_CONSTANTE;
  M->ye = 0.5;
  M->rho = 0.0;
  M->rho_ref = 0.0;
  M->r_ref = 1.0;
  M->expoente = 3.0;
  M->r_start = 0.0;
  M->ruido_rms = 0.0;
  M->ruido_p = 5.0 / 3.0;
  M->ruido_n_min = 1;
  M->ruido_n_max = 200;
  M->ruido_semente = 2026;
}

void dae_nu_ruido_prepara(dae_nu_perfil *M, double *amp, double *fase)
{
  dae_rng g;
  int32_t n, k = 0;
  if (!M || !amp || !fase) return;
  dae_rng_seed(&g, M->ruido_semente);
  for (n = M->ruido_n_min; n <= M->ruido_n_max; ++n, ++k) {
    amp[k] = pow((double)n, -0.5 * M->ruido_p);
    fase[k] = 2.0 * 3.14159265358979323846 * dae_rng_uniform(&g);
  }
  M->ruido_amp = amp;
  M->ruido_fase = fase;
}

/* O fator de flutuacao em `r`. 1.0 quando o ruido esta desligado.
 *
 * As fases vem de dae_rng semeado UMA vez por perfil: o campo tem de ser o
 * mesmo em toda chamada, senao cada fatia da propagacao veria um meio
 * diferente e o resultado nao seria de perfil nenhum. Recalcular a cada
 * chamada custa O(n_max) e e o preco de nao guardar estado — CONVENTIONS.md
 * parte 1 proibe estado global mutavel, e um cache aqui seria exatamente isso. */
static double fator_ruido(const dae_nu_perfil *M, double r_km)
{
  dae_rng g;
  double soma = 0.0, norma = 0.0, L;
  int32_t n;

  if (!(M->ruido_rms > 0.0) || M->ruido_n_max < M->ruido_n_min) return 1.0;
  L = (M->r_tab && M->n_tab > 1) ? (M->r_tab[M->n_tab - 1] - M->r_tab[0])
                                 : (M->r_ref > 0.0 ? M->r_ref : 1.0);
  if (!(L > 0.0)) return 1.0;

  if (M->ruido_amp && M->ruido_fase) {
    int32_t j = 0;
    for (n = M->ruido_n_min; n <= M->ruido_n_max; ++n, ++j) {
      const double a = M->ruido_amp[j];
      const double k = 2.0 * 3.14159265358979323846 * (double)n / L;
      soma += a * cos(k * r_km + M->ruido_fase[j]);
      norma += 0.5 * a * a;
    }
  } else {
    dae_rng_seed(&g, M->ruido_semente);
    for (n = M->ruido_n_min; n <= M->ruido_n_max; ++n) {
      const double a = pow((double)n, -0.5 * M->ruido_p);
      const double fase = 2.0 * 3.14159265358979323846 * dae_rng_uniform(&g);
      const double k = 2.0 * 3.14159265358979323846 * (double)n / L;
      soma += a * cos(k * r_km + fase);
      norma += 0.5 * a * a;               /* <cos^2> = 1/2 */
    }
  }
  if (!(norma > 0.0)) return 1.0;
  /* normaliza para que o contraste tenha o rms pedido */
  return 1.0 + M->ruido_rms * soma / sqrt(norma);
}

/* Busca binaria: o maior i com tab[i] <= x. Devolve -1 antes do inicio. */
static int32_t acha(const double *tab, int32_t n, double x)
{
  int32_t lo = 0, hi = n - 1;
  if (n <= 0 || x < tab[0]) return -1;
  if (x >= tab[n - 1]) return n - 1;
  while (hi - lo > 1) {
    const int32_t m = lo + (hi - lo) / 2;
    if (tab[m] <= x) lo = m; else hi = m;
  }
  return lo;
}

/* Interpola `y` em log-log sobre a grade `x`. Fora da tabela, SEGURA no valor
   da borda em vez de extrapolar: extrapolar lei de potencia num perfil
   estelar erra por ordens de grandeza, e erra em silencio. */
static double interp_loglog(const double *x, const double *y, int32_t n, double xv)
{
  int32_t i;
  double t;
  if (n <= 0) return 0.0;
  if (n == 1) return y[0];
  i = acha(x, n, xv);
  if (i < 0) return y[0];
  if (i >= n - 1) return y[n - 1];
  if (!(x[i] > 0.0) || !(x[i + 1] > 0.0) || !(y[i] > 0.0) || !(y[i + 1] > 0.0)) {
    /* algum ponto nao-positivo: cai para linear, que ao menos e definido */
    t = (xv - x[i]) / (x[i + 1] - x[i]);
    return y[i] + t * (y[i + 1] - y[i]);
  }
  t = (log(xv) - log(x[i])) / (log(x[i + 1]) - log(x[i]));
  return exp(log(y[i]) + t * (log(y[i + 1]) - log(y[i])));
}

double dae_nu_densidade(const dae_nu_perfil *M, double r_km)
{
  double base;
  if (!M) return 0.0;
  if (M->kind == DAE_NU_PERFIL_CONSTANTE) return M->rho * fator_ruido(M, r_km);
  if (M->kind == DAE_NU_PERFIL_TABELA) {
    if (!M->r_tab || !M->rho_tab || M->n_tab <= 0) return 0.0;
    base = interp_loglog(M->r_tab, M->rho_tab, M->n_tab, r_km);
    return base * fator_ruido(M, r_km);
  }
  /* rho(r) = rho_ref * (r_ref / r)^expoente. r <= 0 nao tem sentido fisico e
     daria inf/NaN silencioso; devolve a densidade de referencia. */
  if (!(r_km > 0.0) || !(M->r_ref > 0.0)) return M->rho_ref;
  base = M->rho_ref * pow(M->r_ref / r_km, M->expoente);
  return base * fator_ruido(M, r_km);
}

double dae_nu_ye(const dae_nu_perfil *M, double r_km)
{
  if (!M) return 0.0;
  if (M->kind != DAE_NU_PERFIL_TABELA || !M->ye_tab || !M->r_tab || M->n_tab <= 0)
    return M->ye;
  return interp_loglog(M->r_tab, M->ye_tab, M->n_tab, r_km);
}

/* Um passo de comprimento (r1-r0) km com H congelado na densidade do ponto
 * MEDIO — o ponto medio, e nao a borda, porque com H linear na fatia ele
 * cancela o erro de primeira ordem de graca. */
static dae_status passo_fatia(const dae_nu_params *P, const dae_nu_perfil *M,
                              double r0, double r1, int merg,
                              double *re, double *im)
{
  double hre[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  double him[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  dae_csr A, H;
  dae_cheb W;
  dae_cheb_info info;
  int merg_local = 0;
  dae_status st;

  { const double rm = 0.5 * (r0 + r1);
    st = dae_nu_hamiltoniano(P, dae_nu_densidade(M, rm), dae_nu_ye(M, rm), hre, him); }
  if (st != DAE_OK) return st;
  st = dae_nu_adjacencia(hre, him, P->n_flavors, &A, &merg_local);
  if (st != DAE_OK) return st;
  /* O tamanho do espaco tem de ser o MESMO em toda fatia, senao o vetor de
     estado mudaria de dimensao no meio do caminho. Im(H) nao depende de rho —
     a materia entra real, na diagonal — entao merg e constante ao longo do
     perfil e este guarda nunca deveria disparar. Fica porque se disparar, o
     estado ja esta corrompido e continuar seria pior. */
  if (merg_local != merg) { dae_csr_free(&A); return DAE_ERR_PARAM; }

  st = dae_hamiltonian(&H, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) { dae_csr_free(&A); return st; }
  st = dae_cheb_init(&W, &H, 0);
  if (st != DAE_OK) { dae_csr_free(&A); dae_csr_free(&H); return st; }

  st = dae_cheb_step(&W, &H, (r1 - r0) * DAE_NU_T_POR_KM, re, im, &info);

  dae_cheb_free(&W);
  dae_csr_free(&H);
  dae_csr_free(&A);
  return st;
}

/* Escala espectral de H no raio `r_km`, por Gershgorin. E o `a` que o
 * propagador de Chebyshev usa, e portanto o que fixa alpha = a*dt.
 *
 * Recebe o RAIO, e nao a densidade, para avaliar rho e Y_e pelo MESMO caminho
 * que passo_fatia usa. Usar aqui o `M->ye` escalar seria errado com
 * `ye_table`: os dois divergem, e um escalar menor que a tabela faz o piso
 * subestimar alpha ate o Chebyshev recusar a fatia. Ver t105. */
static double escala_em(const dae_nu_params *P, const dae_nu_perfil *M, double r_km)
{
  double hre[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  double him[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  dae_csr A;
  int merg_local = 0;
  double lo, hi;
  if (dae_nu_hamiltoniano(P, dae_nu_densidade(M, r_km), dae_nu_ye(M, r_km),
                          hre, him) != DAE_OK) return 0.0;
  if (dae_nu_adjacencia(hre, him, P->n_flavors, &A, &merg_local) != DAE_OK) return 0.0;
  dae_csr_gershgorin(&A, &lo, &hi);
  dae_csr_free(&A);
  return 0.5 * (hi - lo);
}

dae_status dae_nu_propaga(const dae_nu_params *P, const dae_nu_perfil *M,
                          double r0_km, double r1_km, double tol,
                          int merg, double *re, double *im, int32_t *subdiv)
{
  /* Teto de subdivisoes: 2^20 sub-passos numa fatia so. Precisa ser folgado
     porque o PISO de factibilidade abaixo ja consome niveis sozinho em regiao
     densa — perto da protoestrela ele pede ~2^17 antes de o criterio de erro
     dizer qualquer coisa. Se nem assim converge, o pedido esta fora do que
     este metodo resolve, e RECUSAR e melhor que devolver um numero. */
  const int32_t TETO = 20;
  const int32_t n = merg ? 2 * P->n_flavors : P->n_flavors;
  double a_re[2 * DAE_NU_MAX_SABORES], a_im[2 * DAE_NU_MAX_SABORES];
  double b_re[2 * DAE_NU_MAX_SABORES], b_im[2 * DAE_NU_MAX_SABORES];
  int32_t nivel, nivel_min = 0, i;
  dae_status st;

  if (!P || !M || !re || !im) return DAE_ERR_PARAM;
  if (subdiv) *subdiv = 1;
  if (!(r1_km > r0_km)) return DAE_OK;              /* fatia vazia */

  /* PISO DE FACTIBILIDADE, antes do criterio de erro. Sao coisas diferentes:
     o criterio de erro pergunta "a resposta esta convergida?", este pergunta
     "o propagador consegue dar UM passo desse tamanho?".
     Perto de uma protoestrela V e ~1e7 vezes maior que Dm2/2E, e um passo de
     30 km daria alpha = a*dt ~ 6e7 — o Chebyshev precisaria de K da mesma
     ordem e devolve DAE_ERR_TOOBIG. Fisicamente esse regime e inofensivo (com
     V dominante o nu_e E autoestado de materia, e a fase rapida nao muda
     observavel nenhum), mas o propagador nao sabe disso: tem de andar em
     passos que ele consiga dar. */
  {
    const double ALPHA_MAX = 200.0;
    const double a0 = escala_em(P, M, r0_km);
    const double a1 = escala_em(P, M, r1_km);
    const double a = (a0 > a1) ? a0 : a1;
    const double alpha = a * (r1_km - r0_km) * DAE_NU_T_POR_KM;
    nivel_min = 0;
    if (alpha > ALPHA_MAX) {
      double razao = alpha / ALPHA_MAX;
      while (razao > 1.0 && nivel_min < TETO) { razao *= 0.5; ++nivel_min; }
    }
  }

  for (nivel = nivel_min; nivel <= TETO; ++nivel) {
    const int32_t k = (int32_t)1 << nivel;          /* 1, 2, 4, ... sub-passos */
    const double h = (r1_km - r0_km) / (double)k;
    double dif = 0.0;

    /* `a` = k sub-passos; `b` = 2k. Os dois partem do MESMO estado de entrada. */
    for (i = 0; i < n; ++i) { a_re[i] = re[i]; a_im[i] = im[i]; }
    for (i = 0; i < k; ++i) {
      st = passo_fatia(P, M, r0_km + (double)i * h, r0_km + (double)(i + 1) * h,
                       merg, a_re, a_im);
      if (st != DAE_OK) return st;
    }
    for (i = 0; i < n; ++i) { b_re[i] = re[i]; b_im[i] = im[i]; }
    for (i = 0; i < 2 * k; ++i) {
      const double h2 = h * 0.5;
      st = passo_fatia(P, M, r0_km + (double)i * h2, r0_km + (double)(i + 1) * h2,
                       merg, b_re, b_im);
      if (st != DAE_OK) return st;
    }

    for (i = 0; i < n; ++i) {
      const double dr = a_re[i] - b_re[i], di = a_im[i] - b_im[i];
      const double d = dr * dr + di * di;
      if (d > dif) dif = d;
    }
    dif = sqrt(dif);

    if (dif <= tol) {
      /* Fica com a resposta MAIS FINA: ela custou o mesmo e e melhor. */
      for (i = 0; i < n; ++i) { re[i] = b_re[i]; im[i] = b_im[i]; }
      if (subdiv) *subdiv = 2 * k;
      return DAE_OK;
    }
  }
  return DAE_ERR_TOOBIG;
}
