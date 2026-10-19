/* dae_numag_janelas.c — integrador por janelas de ressonância.
 * O desenho e o critério estão no cabeçalho de dae_numag.h, seção
 * "integrador por janelas". */
#include "dae_numag.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void dae_numag_janelas_default(dae_numag_janelas_cfg *c)
{
  if (!c) return;
  c->K = 1000.0;
  /* O erro do trecho adiabático segue eps~ nas bordas das janelas (é a
     admixtura que a força bruta acompanha e a fórmula adiabática ignora):
     eps_max controla a precisão. Medido contra força bruta em t110. */
  c->eps_max = 1e-6;
  c->eps_tol = 1e-5;
  c->tol = 1e-10;
  c->tol_fase = 1e-7;
  c->sub = 8;
}

/* O meio e as suas derivadas em s, ANALÍTICOS no intervalo i da tabela — a
   mesma interpolação de dae_numag_meio_tabela (rho geométrico, o resto
   linear), derivada à mão. Dentro de um intervalo tudo é suave; nos nós,
   phi' e Delta_ef podem saltar. */
typedef struct {
  double delta, omega, domega, phi, dphi, ddphi;
  double def, ddef, E, theta, dtheta, eps;
  int    nulo;
} jan_pt;

static void jan_ponto(const dae_numag_trajetoria *T, int32_t i, double s,
                      double cmu, jan_pt *p)
{
  const double TK = DAE_NU_T_POR_KM;
  const double h = T->s_km[i + 1] - T->s_km[i];
  double t = (s - T->s_km[i]) / h;
  double rho, drho, ye, dye, yn, dyn, x, b1, b2, db1, db2, bb, bdb, N, ddelta;

  if (t < 0.0) t = 0.0;
  if (t > 1.0) t = 1.0;
  if (T->rho[i] > 0.0 && T->rho[i + 1] > 0.0) {
    const double lr = log(T->rho[i + 1] / T->rho[i]);
    rho = (t == 0.0) ? T->rho[i] : (t == 1.0) ? T->rho[i + 1]
                                              : T->rho[i] * pow(T->rho[i + 1] / T->rho[i], t);
    drho = rho * lr / h;
  } else {
    rho = (t == 1.0) ? T->rho[i + 1] : T->rho[i] + t * (T->rho[i + 1] - T->rho[i]);
    drho = (T->rho[i + 1] - T->rho[i]) / h;
  }
  ye  = (t == 1.0) ? T->ye[i + 1] : T->ye[i] + t * (T->ye[i + 1] - T->ye[i]);
  dye = (T->ye[i + 1] - T->ye[i]) / h;
  b1  = (t == 1.0) ? T->b1[i + 1] : T->b1[i] + t * (T->b1[i + 1] - T->b1[i]);
  b2  = (t == 1.0) ? T->b2[i + 1] : T->b2[i] + t * (T->b2[i + 1] - T->b2[i]);
  db1 = (T->b1[i + 1] - T->b1[i]) / h;
  db2 = (T->b2[i + 1] - T->b2[i]) / h;

  if (T->ynu) {
    yn  = (t == 1.0) ? T->ynu[i + 1] : T->ynu[i] + t * (T->ynu[i + 1] - T->ynu[i]);
    dyn = (T->ynu[i + 1] - T->ynu[i]) / h;
  } else {
    yn = 0.0; dyn = 0.0;
  }

  /* o mesmo Delta de dae_numag_desajuste: x = (3 Y_e - 1)/2 + 2 Y_nue.
     A ordem das operações na derivada é deliberada: com Y_nue = 0, os
     fatores 0,5 são exatos e o resultado coincide BIT A BIT com o da fórmula
     sem fundo (escrever rho * 1,5 dY_e em vez de 0,5 (3 rho dY_e) mudaria o
     arredondamento). */
  x = (3.0 * ye - 1.0) * 0.5 + 2.0 * yn;
  p->delta = DAE_NU_V_POR_RHO * rho * x;
  ddelta = DAE_NU_V_POR_RHO * 0.5 * (drho * (3.0 * ye - 1.0) + 3.0 * rho * dye)
           + DAE_NU_V_POR_RHO * 2.0 * (drho * yn + rho * dyn);

  bb = sqrt(b1 * b1 + b2 * b2);
  p->nulo = !(bb > 0.0);
  bdb = b1 * db1 + b2 * db2;
  /* Com b linear no intervalo, N = b1 b2' - b2 b1' é CONSTANTE nele: phi' =
     N/|b|^2, e phi'' = -2 N (b.b')/|b|^4. */
  N = b1 * db2 - b2 * db1;
  p->omega = cmu * bb;
  p->domega = p->nulo ? 0.0 : cmu * bdb / bb;
  p->phi = atan2(b2, b1);
  p->dphi = p->nulo ? 0.0 : N / (bb * bb);
  p->ddphi = p->nulo ? 0.0 : -2.0 * N * bdb / (bb * bb * bb * bb);

  p->def = p->delta + p->dphi / TK;               /* eV    */
  p->ddef = ddelta + p->ddphi / TK;               /* eV/km */
  p->E = sqrt(p->def * p->def + 4.0 * p->omega * p->omega);
  p->theta = 0.5 * atan2(2.0 * p->omega, p->def);
  if (p->E > 0.0) {
    p->dtheta = (p->def * p->domega - p->omega * p->ddef) / (p->E * p->E);   /* 1/km */
    p->eps = p->dtheta / (TK * p->E);                                        /* adim. */
  } else {
    p->dtheta = 0.0;
    p->eps = HUGE_VAL;
  }
}

/* O intervalo i com s_i <= s < s_{i+1}, preso a [0, n-2]. */
static int32_t jan_intervalo(const dae_numag_trajetoria *T, double s)
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

/* Gauss-Legendre de 8 pontos de E~ em [a, b], no intervalo i. */
static double jan_gl8(const dae_numag_trajetoria *T, int32_t i, double cmu,
                      double a, double b)
{
  static const double X[4] = { 0.1834346424956498, 0.5255324099163290,
                               0.7966664774136267, 0.9602898564975363 };
  static const double W[4] = { 0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763 };
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  double soma = 0.0;
  int k;
  jan_pt p;
  for (k = 0; k < 4; ++k) {
    jan_ponto(T, i, c - h * X[k], cmu, &p); soma += W[k] * p.E;
    jan_ponto(T, i, c + h * X[k], cmu, &p); soma += W[k] * p.E;
  }
  return soma * h;
}

/* int_a^b E~ ds, adaptativo até que a FASE 0,5 T int E~ ds mude menos que
   `tol` rad entre um nível e o seguinte. O piso relativo 64 eps evita pedir
   à quadratura mais dígitos do que o double carrega numa fase de 1e7 rad. */
static double jan_integra(const dae_numag_trajetoria *T, int32_t i, double cmu,
                          double a, double b, double inteiro, double tol, int prof)
{
  const double TK = DAE_NU_T_POR_KM;
  const double m = 0.5 * (a + b);
  const double l = jan_gl8(T, i, cmu, a, m), r = jan_gl8(T, i, cmu, m, b);
  const double piso = 64.0 * DBL_EPSILON * fabs(l + r) * 0.5 * TK;
  if (prof >= 40 || 0.5 * TK * fabs(l + r - inteiro) <= (tol > piso ? tol : piso))
    return l + r;
  return jan_integra(T, i, cmu, a, m, l, 0.5 * tol, prof + 1)
       + jan_integra(T, i, cmu, m, b, r, 0.5 * tol, prof + 1);
}

/* Trecho adiabático [a, b] dentro do intervalo i, sobre psi FÍSICO (2
   complexos). Laboratório -> girante -> autoestados de H~ -> fases ->
   girante -> laboratório. */
static void jan_adiabatico(const dae_numag_trajetoria *T, int32_t i, double cmu,
                           double a, double b, double tol_fase,
                           double *pr, double *pi)
{
  const double TK = DAE_NU_T_POR_KM, PI = 3.14159265358979323846;
  jan_pt pa, pb;
  double fa, fb, dfi, fase, ca, sa, cb, sb, c, s;
  double t0r, t0i, t1r, t1i, cpr, cpi, cmr, cmi, xr, xi;

  jan_ponto(T, i, a, cmu, &pa);
  jan_ponto(T, i, b, cmu, &pb);
  /* phi contínuo no trecho: a interpolação linear de b não cruza a origem
     num trecho adiabático (B_perp = 0 é janela), então a variação cabe em
     (-pi, pi]. */
  dfi = pb.phi - pa.phi;
  while (dfi > PI) dfi -= 2.0 * PI;
  while (dfi <= -PI) dfi += 2.0 * PI;
  fa = pa.phi; fb = fa + dfi;
  fase = 0.5 * TK * jan_integra(T, i, cmu, a, b, jan_gl8(T, i, cmu, a, b), tol_fase, 0);

  /* psi~ = R(a)^dag psi, R = diag(e^{i phi/2}, e^{-i phi/2}) */
  c = cos(0.5 * fa); s = sin(0.5 * fa);
  t0r = c * pr[0] + s * pi[0]; t0i = c * pi[0] - s * pr[0];      /* e^{-i fa/2} psi0 */
  t1r = c * pr[1] - s * pi[1]; t1i = c * pi[1] + s * pr[1];      /* e^{+i fa/2} psi1 */
  /* c+- = V(a)^T psi~ */
  ca = cos(pa.theta); sa = sin(pa.theta);
  cpr = ca * t0r + sa * t1r; cpi = ca * t0i + sa * t1i;
  cmr = -sa * t0r + ca * t1r; cmi = -sa * t0i + ca * t1i;
  /* c+ e^{-i fase}, c- e^{+i fase} */
  c = cos(fase); s = sin(fase);
  xr = c * cpr + s * cpi; xi = c * cpi - s * cpr; cpr = xr; cpi = xi;
  xr = c * cmr - s * cmi; xi = c * cmi + s * cmr; cmr = xr; cmi = xi;
  /* psi~ = V(b) c */
  cb = cos(pb.theta); sb = sin(pb.theta);
  t0r = cb * cpr - sb * cmr; t0i = cb * cpi - sb * cmi;
  t1r = sb * cpr + cb * cmr; t1i = sb * cpi + cb * cmi;
  /* psi = R(b) psi~ */
  c = cos(0.5 * fb); s = sin(0.5 * fb);
  pr[0] = c * t0r - s * t0i; pi[0] = c * t0i + s * t0r;          /* e^{+i fb/2} */
  pr[1] = c * t1r + s * t1i; pi[1] = c * t1i - s * t1r;          /* e^{-i fb/2} */
}

static void jan_de_mergulho(const double *re, const double *im, double *pr, double *pi)
{
  const double c = 1.0 / sqrt(2.0);
  pr[0] = (re[0] - im[2]) * c; pi[0] = (im[0] + re[2]) * c;
  pr[1] = (re[1] - im[3]) * c; pi[1] = (im[1] + re[3]) * c;
}

static void jan_para_mergulho(const double *pr, const double *pi, double *re, double *im)
{
  const double c = 1.0 / sqrt(2.0);
  int a;
  for (a = 0; a < 2; ++a) {
    re[a] = pr[a] * c; im[a] = pi[a] * c;
    re[2 + a] = pi[a] * c; im[2 + a] = -pr[a] * c;
  }
}

/* ------------------------------------------------------- classificação -- */

/* Fim do bloco de força bruta que começa na célula k e não passa de d: células
   consecutivas juntas até a fase acumulada estimada (E~max . T . comprimento)
   chegar a ALPHA_BLOCO, e NUNCA atravessando um nó da tabela. Chamar
   dae_numag_propaga célula a célula refaria a matriz e o Chebyshev a cada uma;
   atravessar um nó põe a quina do meio interpolado dentro de um passo, a regra
   do ponto médio cai para 1a ordem e o controle de passo subdivide até o teto
   (medido: 2^21 subpassos). */
#define JAN_ALPHA_BLOCO 200.0


typedef struct {
  int32_t n;
  double  *a, *b;
  int32_t *iv;
  int     *jan;
  jan_pt  *L, *M, *R;
} jan_celulas;

static void jan_celulas_free(jan_celulas *C)
{
  free(C->a); free(C->b); free(C->iv); free(C->jan);
  free(C->L); free(C->M); free(C->R);
  memset(C, 0, sizeof(*C));
}

static dae_status jan_classifica(const dae_numag_params *P, const dae_numag_trajetoria *T,
                                 double s0, double s1, const dae_numag_janelas_cfg *cfg,
                                 double cmu, jan_celulas *C, dae_numag_janelas_info *info)
{
  const int32_t i0 = jan_intervalo(T, s0), i1 = jan_intervalo(T, s1);
  const int32_t cap = (i1 - i0 + 1) * (cfg->sub + 1) + 2;
  int32_t i, j, k = 0, c;
  /* Célula VIRTUAL à esquerda de s0, na grade global de subcélulas (a mesma
     que uma chamada única começando antes de s0 usaria): a classificação dos
     zeros e a continuidade dos intervalos nulos não podem depender de onde
     o chamador fatia a trajetória (nt do CLI). No início da
     tabela não há lado esquerdo. */
  jan_pt vL, vM, vR;
  int tem_v = 0;

  (void)P;
  memset(C, 0, sizeof(*C));
  C->a = (double *)malloc((size_t)cap * sizeof(double));
  C->b = (double *)malloc((size_t)cap * sizeof(double));
  C->iv = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  C->jan = (int *)malloc((size_t)cap * sizeof(int));
  C->L = (jan_pt *)malloc((size_t)cap * sizeof(jan_pt));
  C->M = (jan_pt *)malloc((size_t)cap * sizeof(jan_pt));
  C->R = (jan_pt *)malloc((size_t)cap * sizeof(jan_pt));
  if (!C->a || !C->b || !C->iv || !C->jan || !C->L || !C->M || !C->R) {
    jan_celulas_free(C);
    return DAE_ERR_ALLOC;
  }

  for (i = i0; i <= i1; ++i) {
    const double si = T->s_km[i], h = T->s_km[i + 1] - si;
    const double lo = (s0 > si) ? s0 : si;
    const double hi = (s1 < si + h) ? s1 : si + h;
    double x = lo;
    if (!(hi > lo)) continue;
    for (j = 1; j <= cfg->sub; ++j) {
      const double y = (j == cfg->sub) ? hi : si + h * (double)j / (double)cfg->sub;
      if (y <= x) continue;
      if (y > hi) break;
      C->a[k] = x; C->b[k] = y; C->iv[k] = i; ++k;
      x = y;
    }
    if (x < hi) { C->a[k] = x; C->b[k] = hi; C->iv[k] = i; ++k; }
  }
  C->n = k;

  if (s0 > T->s_km[0]) {
    /* intervalo à esquerda: num nó, o anterior (Delta_ef salta nos nós pelo
       termo phi'/T, como em dae_numag_angulos com lado_esquerdo) */
    int32_t il = jan_intervalo(T, s0);
    double hl, g;
    if (il > 0 && s0 <= T->s_km[il]) --il;
    hl = T->s_km[il + 1] - T->s_km[il];
    j = (int32_t)ceil((s0 - T->s_km[il]) / hl * (double)cfg->sub) - 1;
    if (j < 0) j = 0;
    g = T->s_km[il] + hl * (double)j / (double)cfg->sub;
    while (j > 0 && !(g < s0)) { --j; g = T->s_km[il] + hl * (double)j / (double)cfg->sub; }
    if (g < s0) {
      jan_ponto(T, il, g, cmu, &vL);
      jan_ponto(T, il, 0.5 * (g + s0), cmu, &vM);
      jan_ponto(T, il, s0, cmu, &vR);
      tem_v = 1;
    }
  }

  for (c = 0; c < C->n; ++c) {
    jan_pt *L = &C->L[c], *M = &C->M[c], *R = &C->R[c];
    /* a célula anterior na TRAJETÓRIA: a do trecho, ou a virtual em c = 0 */
    const jan_pt *aL = (c > 0) ? &C->L[c - 1] : (tem_v ? &vL : NULL);
    const jan_pt *aM = (c > 0) ? &C->M[c - 1] : (tem_v ? &vM : NULL);
    const jan_pt *aR = (c > 0) ? &C->R[c - 1] : (tem_v ? &vR : NULL);
    const double lim = cfg->K * 2.0;
    jan_ponto(T, C->iv[c], C->a[c], cmu, L);
    jan_ponto(T, C->iv[c], 0.5 * (C->a[c] + C->b[c]), cmu, M);
    jan_ponto(T, C->iv[c], C->b[c], cmu, R);
    C->jan[c] =
      L->nulo || M->nulo || R->nulo ||
      !(L->E > 0.0 && M->E > 0.0 && R->E > 0.0) ||
      fabs(L->def) < lim * L->omega || fabs(M->def) < lim * M->omega ||
      fabs(R->def) < lim * R->omega ||
      fabs(L->eps) > cfg->eps_max || fabs(M->eps) > cfg->eps_max ||
      fabs(R->eps) > cfg->eps_max;
    /* Zeros de Delta_ef. Um zero exatamente num ponto de amostragem é da
       célula à DIREITA dele, para não ser contado duas vezes — inclusive entre
       fatias de uma execução —, e é localizado NO ponto, sem bisseção (que,
       com L->def = 0, levaria o extremo esquerdo até o direito numa raiz
       ascendente). Classificação pelos sinais dos dois lados: troca de sinal
       é cruzamento; mesmo sinal é contato; zero nos três pontos da célula é
       intervalo nulo. O lado esquerdo vem da célula anterior na trajetória
       — em c = 0, da célula virtual antes de s0 —, de modo que uma chamada
       e várias fatias contíguas dão a mesma classificação. Só no início da
       TABELA não há lado esquerdo: um zero ali com R != 0 conta como
       cruzamento (o raio sai da ressonância). Um zero exatamente em s1 não é
       atravessado e não conta: a fatia seguinte o classifica. */
    if (L->def == 0.0) {
      C->jan[c] = 1;
      if (M->def == 0.0 && R->def == 0.0) {
        if (info) {
          if (!(aL && aL->def == 0.0 && aM->def == 0.0 && aR->def == 0.0))
            ++info->n_intervalos_nulos;
          info->comprimento_nulo_km += C->b[c] - C->a[c];
        }
      } else if (info) {
        /* os lados: o valor não nulo mais próximo do zero em cada célula */
        const double dir = (M->def != 0.0) ? M->def : R->def;
        const double esq = !aL ? 0.0 : (aM->def != 0.0) ? aM->def : aL->def;
        if (aL && aL->def == 0.0 && aM->def == 0.0) {
          /* fim de um intervalo nulo: não é um cruzamento isolado novo */
        } else if (!aL || (esq < 0.0) != (dir < 0.0)) {
          ++info->n_cruzamentos; info->cruzamento_s_km = C->a[c];
        } else {
          ++info->n_contatos;
        }
      }
    } else if ((L->def < 0.0 && R->def > 0.0) || (L->def > 0.0 && R->def < 0.0)) {
      double x0 = C->a[c], x1 = C->b[c];
      int it;
      jan_pt q;
      for (it = 0; it < 80; ++it) {
        const double xm = 0.5 * (x0 + x1);
        jan_ponto(T, C->iv[c], xm, cmu, &q);
        if (q.def == 0.0) { x0 = x1 = xm; break; }
        if ((q.def < 0.0) == (L->def < 0.0)) x0 = xm; else x1 = xm;
      }
      C->jan[c] = 1;
      if (info) { ++info->n_cruzamentos; info->cruzamento_s_km = 0.5 * (x0 + x1); }
    }
  }

  /* Trechos adiabáticos: a cota de primeira ordem decide se ficam. */
  c = 0;
  while (c < C->n) {
    int32_t d = c;
    double cota, saltos = 0.0, emax = 0.0;
    if (C->jan[c]) { ++c; continue; }
    while (d + 1 < C->n && !C->jan[d + 1]) ++d;
    cota = fabs(C->L[c].eps) + fabs(C->R[d].eps);
    for (k = c; k <= d; ++k) {
      cota += fabs(C->M[k].eps - C->L[k].eps) + fabs(C->R[k].eps - C->M[k].eps);
      if (fabs(C->L[k].eps) > emax) emax = fabs(C->L[k].eps);
      if (fabs(C->M[k].eps) > emax) emax = fabs(C->M[k].eps);
      if (fabs(C->R[k].eps) > emax) emax = fabs(C->R[k].eps);
      if (k < d && C->iv[k + 1] != C->iv[k]) {
        /* nó da tabela: theta~ e eps~ podem saltar */
        const double dt = fabs(C->L[k + 1].theta - C->R[k].theta);
        saltos += dt;
        cota += dt + fabs(C->L[k + 1].eps - C->R[k].eps);
      }
    }
    if (cota > cfg->eps_tol) {
      for (k = c; k <= d; ++k) C->jan[k] = 1;     /* não confia: força bruta */
    } else if (info) {
      info->cota_adiabatica += cota;
      info->saltos_tabela += saltos;
      if (emax > info->eps_fora_max) info->eps_fora_max = emax;
    }
    c = d + 1;
  }
  return DAE_OK;
}

static int32_t jan_fim_bloco(const jan_celulas *C, int32_t k, int32_t d)
{
  double alpha = 0.0;
  int32_t e = k;
  for (;;) {
    double E = C->L[e].E;
    if (C->M[e].E > E) E = C->M[e].E;
    if (C->R[e].E > E) E = C->R[e].E;
    alpha += 0.5 * E * DAE_NU_T_POR_KM * (C->b[e] - C->a[e]);
    if (e == d || alpha >= JAN_ALPHA_BLOCO || C->iv[e + 1] != C->iv[k]) return e;
    ++e;
  }
}

static dae_status jan_valida(const dae_numag_params *P, const dae_numag_trajetoria *T,
                             double s0, double s1, const dae_numag_janelas_cfg *cfg)
{
  if (!P || !T || !cfg || T->n < 2) return DAE_ERR_PARAM;
  if (P->natureza != DAE_NUMAG_DIRAC || P->sabor.antineutrino || P->n_ativos != 1)
    return DAE_ERR_UNIMPLEMENTED;
  if (!(cfg->K > 0.0) || !(cfg->eps_max > 0.0) || !(cfg->eps_tol > 0.0) ||
      !(cfg->tol > 0.0) || !(cfg->tol_fase > 0.0) || cfg->sub < 1)
    return DAE_ERR_PARAM;
  /* As derivadas analíticas só valem dentro da tabela. Arredondamento de
     s0 + k (t1/nt) não é pedido fora dela: até 1e-12 do comprimento vira a
     borda (ver jan_borda). */
  {
    const double L = T->s_km[T->n - 1] - T->s_km[0], f = 1e-12 * L;
    if (s0 < T->s_km[0] - f || s1 > T->s_km[T->n - 1] + f) return DAE_ERR_PARAM;
  }
  return DAE_OK;
}

static void jan_borda(const dae_numag_trajetoria *T, double *s0, double *s1)
{
  if (*s0 < T->s_km[0]) *s0 = T->s_km[0];
  if (*s1 > T->s_km[T->n - 1]) *s1 = T->s_km[T->n - 1];
}

dae_status dae_numag_propaga_janelas(const dae_numag_params *P,
                                     const dae_numag_trajetoria *T,
                                     double s0, double s1,
                                     const dae_numag_janelas_cfg *cfg,
                                     double *re, double *im,
                                     dae_numag_janelas_info *info)
{
  const double cmu = P ? P->mu_re[0] * DAE_NUMAG_MU_B_EV_POR_T * DAE_NUMAG_T_POR_GAUSS : 0.0;
  dae_numag_trajetoria Tl;
  jan_celulas C;
  dae_status st;
  int32_t c;
  double pr[2], pim[2];

  st = jan_valida(P, T, s0, s1, cfg);
  if (st != DAE_OK) return st;
  if (!re || !im) return DAE_ERR_PARAM;
  jan_borda(T, &s0, &s1);
  if (info) info->terminou_em_janela = 0;
  if (!(s1 > s0)) return DAE_OK;
  Tl = *T;                         /* o callback recebe void* não-const */
  st = jan_classifica(P, T, s0, s1, cfg, cmu, &C, info);
  if (st != DAE_OK) return st;
  if (info && C.n > 0) info->terminou_em_janela = C.jan[C.n - 1];

  c = 0;
  while (c < C.n && st == DAE_OK) {
    int32_t d = c;
    while (d + 1 < C.n && C.jan[d + 1] == C.jan[c]) ++d;
    if (C.jan[c]) {
      int32_t k;
      if (info) { ++info->n_janelas; info->comprimento_janela_km += C.b[d] - C.a[c]; }
      for (k = c; k <= d && st == DAE_OK; ) {
        int32_t sub = 0;
        const int32_t e = jan_fim_bloco(&C, k, d);
        st = dae_numag_propaga(P, dae_numag_meio_tabela, &Tl, C.a[k], C.b[e],
                               cfg->tol, re, im, &sub);
        if (info && sub > info->subdiv_max) info->subdiv_max = sub;
        k = e + 1;
      }
    } else {
      int32_t k = c;
      if (info) ++info->n_trechos;
      jan_de_mergulho(re, im, pr, pim);
      while (k <= d) {
        int32_t e = k;
        while (e + 1 <= d && C.iv[e + 1] == C.iv[k]) ++e;
        jan_adiabatico(T, C.iv[k], cmu, C.a[k], C.b[e], cfg->tol_fase, pr, pim);
        k = e + 1;
      }
      jan_para_mergulho(pr, pim, re, im);
    }
    c = d + 1;
  }
  jan_celulas_free(&C);
  return st;
}

dae_status dae_numag_angulos(const dae_numag_params *P, const dae_numag_trajetoria *T,
                             double s, int lado_esquerdo, double *theta, double *phi)
{
  int32_t i;
  jan_pt p;
  if (!P || !T || T->n < 2 || !theta || !phi) return DAE_ERR_PARAM;
  i = jan_intervalo(T, s);
  if (lado_esquerdo && i > 0 && s == T->s_km[i]) --i;
  jan_ponto(T, i, s, P->mu_re[0] * DAE_NUMAG_MU_B_EV_POR_T * DAE_NUMAG_T_POR_GAUSS, &p);
  *theta = p.theta;
  *phi = p.phi;
  return DAE_OK;
}

/* Aplica a força bruta das células c..d às duas colunas de U (laboratório),
   no lugar. */
static dae_status jan_U_janela(const dae_numag_params *P, dae_numag_trajetoria *Tl,
                               const jan_celulas *C, int32_t c, int32_t d,
                               const dae_numag_janelas_cfg *cfg, double ur[4], double ui[4],
                               dae_numag_janelas_info *info)
{
  int col;
  dae_status st = DAE_OK;
  for (col = 0; col < 2 && st == DAE_OK; ++col) {
    double re[4], im[4], pr[2], pim[2];
    int32_t k;
    pr[0] = ur[col]; pim[0] = ui[col]; pr[1] = ur[2 + col]; pim[1] = ui[2 + col];
    jan_para_mergulho(pr, pim, re, im);
    for (k = c; k <= d && st == DAE_OK; ) {
      int32_t sub = 0;
      const int32_t e = jan_fim_bloco(C, k, d);
      st = dae_numag_propaga(P, dae_numag_meio_tabela, Tl, C->a[k], C->b[e],
                             cfg->tol, re, im, &sub);
      if (info && sub > info->subdiv_max) info->subdiv_max = sub;
      k = e + 1;
    }
    jan_de_mergulho(re, im, pr, pim);
    ur[col] = pr[0]; ui[col] = pim[0]; ur[2 + col] = pr[1]; ui[2 + col] = pim[1];
  }
  return st;
}

/* prob[m*2+n] = |<m~(b)| R(b)^dag U R(a) |n~(a)>|^2 */
static void jan_prob(const double ur[4], const double ui[4], double ta, double fa,
                     double tb, double fb, double prob[4])
{
  int m, n;
  for (m = 0; m < 2; ++m) for (n = 0; n < 2; ++n) {
    const double na0 = (n == 0) ? cos(ta) : -sin(ta), na1 = (n == 0) ? sin(ta) : cos(ta);
    const double mb0 = (m == 0) ? cos(tb) : -sin(tb), mb1 = (m == 0) ? sin(tb) : cos(tb);
    const double v0r = cos(0.5 * fa) * na0, v0i = sin(0.5 * fa) * na0;
    const double v1r = cos(0.5 * fa) * na1, v1i = -sin(0.5 * fa) * na1;
    const double w0r = ur[0] * v0r - ui[0] * v0i + ur[1] * v1r - ui[1] * v1i;
    const double w0i = ur[0] * v0i + ui[0] * v0r + ur[1] * v1i + ui[1] * v1r;
    const double w1r = ur[2] * v0r - ui[2] * v0i + ur[3] * v1r - ui[3] * v1i;
    const double w1i = ur[2] * v0i + ui[2] * v0r + ur[3] * v1i + ui[3] * v1r;
    const double cb = cos(0.5 * fb), sb = sin(0.5 * fb);
    const double z0r = cb * w0r + sb * w0i, z0i = cb * w0i - sb * w0r;
    const double z1r = cb * w1r - sb * w1i, z1i = cb * w1i + sb * w1r;
    const double ar = mb0 * z0r + mb1 * z1r, ai = mb0 * z0i + mb1 * z1i;
    prob[m * 2 + n] = ar * ar + ai * ai;
  }
}

static void jan_aplica(const double prob[4], const double in[2], double out[2])
{
  const double p0 = in[0], p1 = in[1];
  out[0] = prob[0] * p0 + prob[1] * p1;
  out[1] = prob[2] * p0 + prob[3] * p1;
}

static void jan_mistura(double *pop, double dt)
{
  if (dt != 0.0) {
    const double c2 = cos(dt) * cos(dt), s2 = 1.0 - c2;
    const double p0 = pop[0], p1 = pop[1];
    pop[0] = c2 * p0 + s2 * p1;
    pop[1] = s2 * p0 + c2 * p1;
  }
}

dae_status dae_numag_propaga_janelas_inc(const dae_numag_params *P,
                                         const dae_numag_trajetoria *T,
                                         double s0, double s1,
                                         const dae_numag_janelas_cfg *cfg,
                                         double *pop,
                                         dae_numag_janela_aberta *ab,
                                         dae_numag_janelas_info *info)
{
  const double cmu = P ? P->mu_re[0] * DAE_NUMAG_MU_B_EV_POR_T * DAE_NUMAG_T_POR_GAUSS : 0.0;
  dae_numag_trajetoria Tl;
  dae_numag_janela_aberta local;
  jan_celulas C;
  dae_status st;
  int32_t c;
  double theta_atual;

  st = jan_valida(P, T, s0, s1, cfg);
  if (st != DAE_OK) return st;
  if (!pop) return DAE_ERR_PARAM;
  jan_borda(T, &s0, &s1);
  if (info) info->terminou_em_janela = 0;
  if (!ab) { memset(&local, 0, sizeof local); ab = &local; }
  if (!(s1 > s0)) return DAE_OK;
  Tl = *T;
  st = jan_classifica(P, T, s0, s1, cfg, cmu, &C, info);
  if (st != DAE_OK) return st;
  if (info && C.n > 0) info->terminou_em_janela = C.jan[C.n - 1];

  /* Janela aberta numa chamada anterior e que NÃO continua aqui: fecha em s0. */
  if (ab->aberta && C.n > 0 && !C.jan[0]) {
    double prob[4];
    jan_prob(ab->ur, ab->ui, ab->theta_a, ab->phi_a, C.L[0].theta, C.L[0].phi, prob);
    jan_aplica(prob, ab->pop_a, pop);
    ab->aberta = 0;
  }

  theta_atual = C.n ? C.L[0].theta : 0.0;
  c = 0;
  while (c < C.n && st == DAE_OK) {
    int32_t d = c;
    const int continua = (c == 0 && ab->aberta && C.jan[0]);
    if (!continua) {
      jan_mistura(pop, C.L[c].theta - theta_atual);    /* nó: base muda de vez */
      theta_atual = C.L[c].theta;
    }
    while (d + 1 < C.n && C.jan[d + 1] == C.jan[c]) ++d;
    if (C.jan[c]) {
      double prob[4];
      if (info) { ++info->n_janelas; info->comprimento_janela_km += C.b[d] - C.a[c]; }
      if (!continua) {
        /* janela nova: U parcial = identidade, base e populações do início */
        ab->ur[0] = 1.0; ab->ui[0] = 0.0; ab->ur[1] = 0.0; ab->ui[1] = 0.0;
        ab->ur[2] = 0.0; ab->ui[2] = 0.0; ab->ur[3] = 1.0; ab->ui[3] = 0.0;
        ab->theta_a = C.L[c].theta; ab->phi_a = C.L[c].phi;
        ab->pop_a[0] = pop[0]; ab->pop_a[1] = pop[1];
        ab->aberta = 1;
      }
      st = jan_U_janela(P, &Tl, &C, c, d, cfg, ab->ur, ab->ui, info);
      if (st != DAE_OK) break;
      jan_prob(ab->ur, ab->ui, ab->theta_a, ab->phi_a, C.R[d].theta, C.R[d].phi, prob);
      jan_aplica(prob, ab->pop_a, pop);
      /* Termina dentro da trajetória pedida: fecha. Termina em s1: fica
         aberta (se houver quem continue) e `pop` é a provisória em s1. */
      if (d < C.n - 1 || ab == &local) ab->aberta = 0;
      theta_atual = C.R[d].theta;
    } else {
      int32_t k;
      if (info) ++info->n_trechos;
      for (k = c; k <= d; ++k) {
        jan_mistura(pop, C.L[k].theta - theta_atual);
        theta_atual = C.R[k].theta;   /* dentro da célula a base gira de leve
                                         e as populações ficam */
      }
    }
    c = d + 1;
  }
  jan_celulas_free(&C);
  return st;
}
