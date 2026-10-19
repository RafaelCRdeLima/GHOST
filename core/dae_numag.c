#include "dae_numag.h"

#include "dae_cheb.h"
#include "dae_ham.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void dae_numag_params_default(dae_numag_params *P)
{
  if (!P) return;
  memset(P, 0, sizeof(*P));
  dae_nu_params_default(&P->sabor);
  P->sabor.energy_eV = 10.0e6;
  P->natureza = DAE_NUMAG_DIRAC;
  P->n_ativos = 1;
}

static double norma3(const double v[3])
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void vetorial(const double a[3], const double b[3], double c[3])
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

dae_status dae_numag_base_transversa(const double k[3], double e1[3], double e2[3])
{
  static const double Z[3] = { 0.0, 0.0, 1.0 };
  static const double X[3] = { 1.0, 0.0, 0.0 };
  double kh[3], nk, n1;
  int i;

  if (!k || !e1 || !e2) return DAE_ERR_PARAM;
  nk = norma3(k);
  if (!(nk > 0.0)) return DAE_ERR_PARAM;          /* também recusa NaN */
  for (i = 0; i < 3; ++i) kh[i] = k[i] / nk;

  /* z x k, a menos que k esteja perto de z: aí x x k. O limiar só precisa
     evitar dividir por quase zero; a base é qualquer uma, desde que FIXA para
     uma dada direção, e é. */
  vetorial(Z, kh, e1);
  n1 = norma3(e1);
  if (n1 < 1e-3) { vetorial(X, kh, e1); n1 = norma3(e1); }
  for (i = 0; i < 3; ++i) e1[i] /= n1;
  vetorial(kh, e1, e2);
  return DAE_OK;
}

void dae_numag_projeta_campo(const double b[3], const double e1[3],
                             const double e2[3], double b_t[2])
{
  if (!b || !e1 || !e2 || !b_t) return;
  /* As componentes de B em e1 e e2 SÃO a parte transversal: a componente ao
     longo de k não tem projeção nelas. */
  b_t[0] = b[0] * e1[0] + b[1] * e1[1] + b[2] * e1[2];
  b_t[1] = b[0] * e2[0] + b[1] * e2[1] + b[2] * e2[2];
}

dae_status dae_numag_acoplamento(const dae_numag_meio *m,
                                 double *b_perp_gauss, double *phi_rad)
{
  double b1, b2;

  if (!m || !b_perp_gauss || !phi_rad) return DAE_ERR_PARAM;
  b1 = m->b_t[0];
  b2 = m->b_t[1];
  if (b1 != b1 || b2 != b2) return DAE_ERR_PARAM;
  *b_perp_gauss = sqrt(b1 * b1 + b2 * b2);
  *phi_rad = atan2(b2, b1);
  return DAE_OK;
}

double dae_numag_desajuste(const dae_numag_meio *m)
{
  if (!m) return 0.0;
  /* sqrt(2) G_F (n_e - n_n/2 + 2 (n_nue - n_nuebar)), com n_e = Y_e n_b e
     n_n = (1 - Y_e) n_b em matéria neutra de e, p, n. DAE_NU_V_POR_RHO é
     sqrt(2) G_F n_b por rho. */
  return DAE_NU_V_POR_RHO * m->rho * ((3.0 * m->ye - 1.0) * 0.5 + 2.0 * m->ynu);
}

dae_status dae_numag_hamiltoniano(const dae_numag_params *P,
                                  const dae_numag_meio *m,
                                  double *hre, double *him)
{
  double bperp, phi, omega, delta;
  dae_status st;

  if (!P || !m || !hre || !him) return DAE_ERR_PARAM;
  if (P->natureza != DAE_NUMAG_DIRAC || P->sabor.antineutrino || P->n_ativos != 1)
    return DAE_ERR_UNIMPLEMENTED;
  if (!(m->rho >= 0.0) || !(m->ye >= 0.0)) return DAE_ERR_PARAM;   /* e NaN */
  if (!(fabs(m->ynu) <= 1.0)) return DAE_ERR_PARAM;   /* |Y_nue| <= 1, e NaN */
  /* Momento diagonal de Dirac é REAL: uma parte imaginária aqui é erro de
     quem montou P, não fase física. */
  if (P->mu_im[0] != 0.0 || P->mu_re[0] != P->mu_re[0]) return DAE_ERR_PARAM;

  st = dae_numag_acoplamento(m, &bperp, &phi);
  if (st != DAE_OK) return st;

  omega = P->mu_re[0] * DAE_NUMAG_MU_B_EV_POR_T * bperp * DAE_NUMAG_T_POR_GAUSS;
  delta = dae_numag_desajuste(m);

  /* Traço do ESPAÇO INTEIRO removido: diag(Delta/2, -Delta/2). Removê-lo por
     bloco daria diag(0, 0) e apagaria a ressonância (CONVENTIONS.md 13.4). */
  hre[0] = 0.5 * delta;  him[0] = 0.0;
  hre[3] = -0.5 * delta; him[3] = 0.0;
  /* H_LR = Omega e^{+i phi}; H_RL é o conjugado. */
  hre[1] = omega * cos(phi); him[1] = omega * sin(phi);
  hre[2] = hre[1];           him[2] = -him[1];
  return DAE_OK;
}

dae_status dae_numag_adjacencia(const double *hre, const double *him, int32_t n,
                                dae_csr *A)
{
  enum { NMAX = 2 * DAE_NUMAG_MAX_ESTADOS, CAP = NMAX * (NMAX + 1) / 2 };
  int32_t ei[CAP], ej[CAP];
  double w[CAP];
  int32_t ne = 0, a, b;

  if (!hre || !him || !A || n < 1 || n > DAE_NUMAG_MAX_ESTADOS) return DAE_ERR_PARAM;

  /* H~ = [X -Y; Y X] com H = X + iY; triângulo superior, A = -H~. É o mesmo
     mergulho de dae_nu_adjacencia, mas SEMPRE ligado. */
  for (a = 0; a < 2 * n; ++a) {
    for (b = a; b < 2 * n; ++b) {
      double v;
      if (a < n && b < n)  v = hre[a * n + b];
      else if (a < n)      v = -him[a * n + (b - n)];
      else                 v = hre[(a - n) * n + (b - n)];
      v = -v;
      if (v != 0.0) { ei[ne] = a; ej[ne] = b; w[ne] = v; ++ne; }
    }
  }
  /* ne = 0 (H nulo: mu = 0 em Y_e = 1/3) é aceito e mantém a dimensão 2n. */
  return dae_csr_from_edges(A, 2 * n, ei, ej, w, ne, NULL);
}

/* Um passo de (s1 - s0) km com H congelado no PONTO MÉDIO, como em
   dae_nu_propaga: com H linear na fatia, o ponto médio cancela o erro de
   primeira ordem. */
static dae_status numag_passo_fatia(const dae_numag_params *P, dae_numag_meio_fn meio,
                              void *user, double s0, double s1,
                              double *re, double *im)
{
  double hre[DAE_NUMAG_MAX_ESTADOS * DAE_NUMAG_MAX_ESTADOS];
  double him[DAE_NUMAG_MAX_ESTADOS * DAE_NUMAG_MAX_ESTADOS];
  dae_numag_meio m;
  dae_csr A, H;
  dae_cheb W;
  dae_cheb_info info;
  dae_status st;

  memset(&m, 0, sizeof(m));
  memset(&m, 0, sizeof(m));
  meio(user, 0.5 * (s0 + s1), &m);
  st = dae_numag_hamiltoniano(P, &m, hre, him);
  if (st != DAE_OK) return st;
  st = dae_numag_adjacencia(hre, him, 2 * P->n_ativos, &A);
  if (st != DAE_OK) return st;
  /* H = 0 (mu = 0 em Y_e = 1/3): exp(-i 0 t) = 1. Sem isto, o Chebyshev
     força a escala espectral de uma matriz nula a 1, o piso (que viu escala
     0) não subdivide, e um passo de alpha ~ 1e6 é recusado como TOOBIG. */
  if (A.nnz == 0) { dae_csr_free(&A); return DAE_OK; }
  st = dae_hamiltonian(&H, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) { dae_csr_free(&A); return st; }
  st = dae_cheb_init(&W, &H, 0);
  if (st != DAE_OK) { dae_csr_free(&A); dae_csr_free(&H); return st; }
  st = dae_cheb_step(&W, &H, (s1 - s0) * DAE_NU_T_POR_KM, re, im, &info);
  dae_cheb_free(&W);
  dae_csr_free(&H);
  dae_csr_free(&A);
  return st;
}

/* Raio espectral de H no ponto s, por Gershgorin sobre a adjacência — o `a`
   que fixa alpha = a dt. Lê o meio pela MESMA função da propagação. */
static double numag_escala_em(const dae_numag_params *P, dae_numag_meio_fn meio,
                        void *user, double s)
{
  double hre[DAE_NUMAG_MAX_ESTADOS * DAE_NUMAG_MAX_ESTADOS];
  double him[DAE_NUMAG_MAX_ESTADOS * DAE_NUMAG_MAX_ESTADOS];
  dae_numag_meio m;
  dae_csr A;
  double lo, hi;

  memset(&m, 0, sizeof(m));
  memset(&m, 0, sizeof(m));
  meio(user, s, &m);
  if (dae_numag_hamiltoniano(P, &m, hre, him) != DAE_OK) return 0.0;
  if (dae_numag_adjacencia(hre, him, 2 * P->n_ativos, &A) != DAE_OK) return 0.0;
  dae_csr_gershgorin(&A, &lo, &hi);
  dae_csr_free(&A);
  return 0.5 * (hi - lo);
}

dae_status dae_numag_propaga(const dae_numag_params *P, dae_numag_meio_fn meio,
                             void *user, double s0_km, double s1_km, double tol,
                             double *re, double *im, int32_t *subdiv)
{
  /* Os mesmos tetos de dae_nu_propaga, pelas mesmas razões. */
  const int32_t TETO = 20;
  const double ALPHA_MAX = 200.0;
  double a_re[2 * DAE_NUMAG_MAX_ESTADOS], a_im[2 * DAE_NUMAG_MAX_ESTADOS];
  double b_re[2 * DAE_NUMAG_MAX_ESTADOS], b_im[2 * DAE_NUMAG_MAX_ESTADOS];
  int32_t n, nivel, nivel_min = 0, i;
  dae_status st;

  if (!P || !meio || !re || !im || !(tol > 0.0)) return DAE_ERR_PARAM;
  if (P->n_ativos < 1 || P->n_ativos > DAE_NUMAG_MAX_ATIVOS) return DAE_ERR_PARAM;
  n = 4 * P->n_ativos;                      /* 2 helicidades, mergulhado */
  if (subdiv) *subdiv = 1;
  if (!(s1_km > s0_km)) return DAE_OK;

  {
    const double a0 = numag_escala_em(P, meio, user, s0_km);
    const double a1 = numag_escala_em(P, meio, user, s1_km);
    const double a = (a0 > a1) ? a0 : a1;
    const double alpha = a * (s1_km - s0_km) * DAE_NU_T_POR_KM;
    if (alpha > ALPHA_MAX) {
      double razao = alpha / ALPHA_MAX;
      while (razao > 1.0 && nivel_min < TETO) { razao *= 0.5; ++nivel_min; }
    }
  }

  for (nivel = nivel_min; nivel <= TETO; ++nivel) {
    const int32_t k = (int32_t)1 << nivel;
    const double h = (s1_km - s0_km) / (double)k;
    double dif = 0.0;

    for (i = 0; i < n; ++i) { a_re[i] = re[i]; a_im[i] = im[i]; }
    for (i = 0; i < k; ++i) {
      st = numag_passo_fatia(P, meio, user, s0_km + (double)i * h,
                       s0_km + (double)(i + 1) * h, a_re, a_im);
      if (st != DAE_OK) return st;
    }
    for (i = 0; i < n; ++i) { b_re[i] = re[i]; b_im[i] = im[i]; }
    for (i = 0; i < 2 * k; ++i) {
      const double h2 = 0.5 * h;
      st = numag_passo_fatia(P, meio, user, s0_km + (double)i * h2,
                       s0_km + (double)(i + 1) * h2, b_re, b_im);
      if (st != DAE_OK) return st;
    }
    for (i = 0; i < n; ++i) {
      const double dr = a_re[i] - b_re[i], di = a_im[i] - b_im[i];
      const double d = dr * dr + di * di;
      if (d > dif) dif = d;
    }
    if (sqrt(dif) <= tol) {
      for (i = 0; i < n; ++i) { re[i] = b_re[i]; im[i] = b_im[i]; }
      if (subdiv) *subdiv = 2 * k;
      return DAE_OK;
    }
  }
  return DAE_ERR_TOOBIG;
}

/* ------------------------------------------------- trajetória tabelada -- */

/* O maior i com s[i] <= x; -1 antes do início. Busca binária. */
static int32_t numag_acha(const double *s, int32_t n, double x)
{
  int32_t lo = 0, hi = n - 1;
  if (n <= 0 || x < s[0]) return -1;
  if (x >= s[n - 1]) return n - 1;
  while (hi - lo > 1) {
    const int32_t m = lo + (hi - lo) / 2;
    if (s[m] <= x) lo = m; else hi = m;
  }
  return lo;
}

/* Exata nos dois nós: y[i] + 1 * (y[i+1] - y[i]) não devolve y[i+1] bit a
   bit em geral. */
static double numag_lin(const double *y, int32_t i, double t)
{
  if (t == 1.0) return y[i + 1];
  return y[i] + t * (y[i + 1] - y[i]);
}

void dae_numag_meio_tabela(void *user, double s_km, dae_numag_meio *m)
{
  const dae_numag_trajetoria *T = (const dae_numag_trajetoria *)user;
  int32_t i;
  double t;

  if (!T || !m || T->n <= 0) return;
  if (T->n == 1 || s_km <= T->s_km[0]) { i = 0; t = 0.0; }
  else {
    i = numag_acha(T->s_km, T->n, s_km);
    if (i >= T->n - 1) { i = T->n - 2; t = 1.0; }      /* segura na borda */
    else t = (s_km - T->s_km[i]) / (T->s_km[i + 1] - T->s_km[i]);
  }
  if (T->n == 1) {
    m->rho = T->rho[0]; m->ye = T->ye[0];
    m->b_t[0] = T->b1[0]; m->b_t[1] = T->b2[0];
    m->razao_energia = T->razao_energia ? T->razao_energia[0] : 1.0;
    m->ynu = T->ynu ? T->ynu[0] : 0.0;
    return;
  }
  /* Nos nós, o valor TABELADO, exato: exp(log(rho)) não devolve rho bit a
     bit, e o que o spec diz num nó é o que vale ali. */
  if (t == 0.0)      m->rho = T->rho[i];
  else if (t == 1.0) m->rho = T->rho[i + 1];
  else if (T->rho[i] > 0.0 && T->rho[i + 1] > 0.0)
    m->rho = T->rho[i] * pow(T->rho[i + 1] / T->rho[i], t);
  else
    m->rho = numag_lin(T->rho, i, t);
  m->ye = numag_lin(T->ye, i, t);
  m->b_t[0] = numag_lin(T->b1, i, t);
  m->b_t[1] = numag_lin(T->b2, i, t);
  m->razao_energia = T->razao_energia ? numag_lin(T->razao_energia, i, t) : 1.0;
  m->ynu = T->ynu ? numag_lin(T->ynu, i, t) : 0.0;
}

dae_status dae_numag_propaga_tabela(const dae_numag_params *P,
                                    const dae_numag_trajetoria *T,
                                    double s0_km, double s1_km, double tol,
                                    double *re, double *im, int32_t *subdiv)
{
  dae_numag_trajetoria Tl;
  double a = s0_km;
  int32_t i, maximo = 0, sub;
  dae_status st;

  if (!T || T->n < 2) return DAE_ERR_PARAM;
  Tl = *T;
  if (subdiv) *subdiv = 1;
  i = numag_acha(T->s_km, T->n, s0_km);
  if (i < 0) i = 0;
  for (; i + 1 < T->n && T->s_km[i + 1] < s1_km; ++i) {
    if (T->s_km[i + 1] <= a) continue;
    sub = 0;
    st = dae_numag_propaga(P, dae_numag_meio_tabela, &Tl, a, T->s_km[i + 1], tol, re, im, &sub);
    if (st != DAE_OK) return st;
    if (sub > maximo) maximo = sub;
    a = T->s_km[i + 1];
  }
  sub = 0;
  st = dae_numag_propaga(P, dae_numag_meio_tabela, &Tl, a, s1_km, tol, re, im, &sub);
  if (st != DAE_OK) return st;
  if (sub > maximo) maximo = sub;
  if (subdiv) *subdiv = maximo;
  return DAE_OK;
}
