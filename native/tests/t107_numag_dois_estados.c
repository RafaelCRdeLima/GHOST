/* t107_numag_dois_estados.c — ACEITAÇÃO 107: o módulo dae_numag, um sabor
 * Dirac (nu_eL, nu_eR), contra as respostas exatas.
 *
 * O que cada bloco prova, e a companheira que impede o teste de ser vazio:
 *
 *  1. CAMADA CONSTANTE contra P = 4 Omega^2/(Delta^2 + 4 Omega^2) sin^2(...),
 *     Delta/(2 Omega) = 0, 1, 5, em unidades físicas (rho, Y_e, mu_nu, B, km).
 *     Companheira: Delta = 0 atinge P = 1.
 *  2. RAMPA LINEAR contra Landau-Zener, P_{L->R} = 1 - exp(-gamma), para
 *     gamma = 0,1, 1 e 3, com extremos a 4000 Omega da ressonância.
 *     Companheira: extremos a 2 Omega divergem de LZ de forma mensurável.
 *  3. mu = 0 não converte (exatamente); mu != 0 converte.
 *  4. B paralelo a k não converte; B perpendicular converte.
 *  5. Hermiticidade, traço nulo no espaço inteiro, norma.
 *  6. POSIÇÃO DA RESSONÂNCIA: Y_e = 1/3 converte tudo, Y_e = 0,1 dá a
 *     supressão P_max = 1,91e-11. Um traço removido POR BLOCO
 *     zeraria Delta e levaria a amplitude a 1 em Y_e = 0,1; o mutante que faz
 *     isso falha aqui (e nos blocos 1, 2, 5 e 7).
 *  7. CAMPO GIRANTE e o SINAL de Delta_ef (CONVENTIONS.md 13.3). Com
 *     phi_B = phi' s, a resposta é a de camada constante com
 *     Delta_ef = Delta + phi'/DAE_NU_T_POR_KM. Exige ressonância completa em
 *     Delta = -phi'/T, onde Delta NÃO se anula, e supressão em Delta = +phi'/T.
 *     Se o sinal estivesse trocado, as duas exigências falhariam juntas.
 *     Companheira: phi_B constante e não nulo dá o mesmo P que phi_B = 0.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

/* rho fixo, alto: Delta = V rho (3 Y_e - 1)/2 cobre +-3,8 eV com Y_e em [0, 1]. */
#define RHO 1.0e14
#define VR  (DAE_NU_V_POR_RHO * RHO)

typedef struct {
  double ye0, dye_ds, s_ref;   /* Y_e(s) = ye0 + dye_ds (s - s_ref)          */
  double b_perp, b_par;        /* gauss                                       */
  double phi0, dphi_ds;        /* rad, rad/km                                 */
  double k[3];
} caminho;

static void meio_fn(void *user, double s, dae_numag_meio *m)
{
  const caminho *c = (const caminho *)user;
  double e1[3], e2[3], kh[3], nk;
  const double phi = c->phi0 + c->dphi_ds * s;
  int i;
  dae_numag_base_transversa(c->k, e1, e2);
  nk = sqrt(c->k[0] * c->k[0] + c->k[1] * c->k[1] + c->k[2] * c->k[2]);
  for (i = 0; i < 3; ++i) kh[i] = c->k[i] / nk;
  double b[3];
  m->rho = RHO;
  m->ye = c->ye0 + c->dye_ds * (s - c->s_ref);
  m->razao_energia = 1.0;
  m->ynu = 0.0;
  for (i = 0; i < 3; ++i)
    b[i] = c->b_perp * (cos(phi) * e1[i] + sin(phi) * e2[i]) + c->b_par * kh[i];
  /* Raio RETO: k fixo, e a base fixa de dae_numag_base_transversa é a
     transportada. É o único caso em que ela pode ser usada. */
  dae_numag_projeta_campo(b, e1, e2, m->b_t);
}

/* Y_e que dá o desajuste Delta (eV) com rho = RHO. */
static double ye_para(double delta) { return (2.0 * delta / VR + 1.0) / 3.0; }

static double exata(double delta, double omega, double t)
{
  const double w = sqrt(delta * delta + 4.0 * omega * omega);
  const double s = sin(0.5 * t * w);
  return (w == 0.0) ? 0.0 : 4.0 * omega * omega / (w * w) * s * s;
}

static void caminho_padrao(caminho *c)
{
  c->ye0 = 1.0 / 3.0; c->dye_ds = 0.0; c->s_ref = 0.0;
  c->b_perp = 1.0e15; c->b_par = 0.0;
  c->phi0 = 0.0; c->dphi_ds = 0.0;
  c->k[0] = 1.0; c->k[1] = 2.0; c->k[2] = 2.0;      /* oblíquo de propósito */
}

/* nu_eL em s0, propagado até s1 em `nf` fatias. Devolve P_{L->R} e a norma. */
static dae_status roda(const dae_numag_params *P, caminho *c, double s0, double s1,
                       int32_t nf, double tol, double *p_r, double *norma)
{
  double re[4], im[4], pre[2], pim[2];
  int32_t i;
  dae_status st = DAE_OK;
  dae_nu_inicial(0, 2, 1, re, im);
  for (i = 0; i < nf && st == DAE_OK; ++i)
    st = dae_numag_propaga(P, meio_fn, c, s0 + (s1 - s0) * i / nf,
                           s0 + (s1 - s0) * (i + 1) / nf, tol, re, im, NULL);
  if (st != DAE_OK) return st;
  dae_nu_projeta(re, im, 2, 1, pre, pim);
  *p_r = pre[1] * pre[1] + pim[1] * pim[1];
  *norma = *p_r + pre[0] * pre[0] + pim[0] * pim[0];
  return DAE_OK;
}

int main(void)
{
  dae_test T;
  dae_numag_params P;
  caminho c;
  double omega, p, nrm, pior, pior_norma = 0.0;
  const double TK = DAE_NU_T_POR_KM;
  int a, b;
  dae_status st;

  dae_test_begin(&T, "t107 dae_numag: dois estados Dirac contra as respostas exatas");
  dae_numag_params_default(&P);
  P.mu_re[0] = 1.0e-12;                                     /* mu_B */
  omega = 1.0e-12 * DAE_NUMAG_MU_B_EV_POR_T * 1.0e15 * DAE_NUMAG_T_POR_GAUSS;

  /* ---------------------------------------------- 1. camada constante ---- */
  {
    static const double Q[3] = { 0.0, 1.0, 5.0 };
    static const double X[4] = { 0.25, 0.5, 0.9, 1.7 };
    pior = 0.0;
    for (a = 0; a < 3; ++a) for (b = 0; b < 4; ++b) {
      const double delta = 2.0 * omega * Q[a];
      const double L = X[b] * PI / omega / TK;              /* km */
      caminho_padrao(&c); c.ye0 = ye_para(delta);
      st = roda(&P, &c, 0.0, L, 1, 1e-12, &p, &nrm);
      if (st != DAE_OK) { dae_test_ok(&T, 0, "camada: %s", dae_strerror(st)); return dae_test_end(&T); }
      if (fabs(p - exata(delta, omega, L * TK)) > pior) pior = fabs(p - exata(delta, omega, L * TK));
      if (fabs(nrm - 1.0) > pior_norma) pior_norma = fabs(nrm - 1.0);
    }
    dae_test_ok(&T, pior < 1e-10, "1. camada constante, q = 0, 1, 5 (pior %.2e; 12 pontos)", pior);
    caminho_padrao(&c);
    roda(&P, &c, 0.0, 0.5 * PI / omega / TK, 1, 1e-12, &p, &nrm);
    dae_test_near(&T, p, 1.0, 1e-10, "1. companheira: Delta = 0 atinge P = 1 em L = %.3f cm",
                  0.5 * PI / omega / TK * 1e5);
  }

  /* ------------------------------------------------ 2. rampa linear ---- */
  {
    static const double G[3] = { 0.1, 1.0, 3.0 };
    const double K = 4000.0;
    pior = 0.0;
    for (a = 0; a < 3; ++a) {
      const double dD_nat = 2.0 * PI * omega * omega / G[a];   /* eV^2   */
      const double dD = dD_nat * TK;                            /* eV/km  */
      const double meia = K * 2.0 * omega / dD;                 /* km     */
      double d;
      caminho_padrao(&c);
      c.ye0 = 1.0 / 3.0; c.dye_ds = dD / (1.5 * VR); c.s_ref = 0.0;
      st = roda(&P, &c, -meia, meia, 400, 1e-10, &p, &nrm);
      if (st != DAE_OK) { dae_test_ok(&T, 0, "rampa: %s", dae_strerror(st)); return dae_test_end(&T); }
      d = fabs(p - (1.0 - exp(-G[a])));
      if (d > pior) pior = d;
      if (fabs(nrm - 1.0) > pior_norma) pior_norma = fabs(nrm - 1.0);
      dae_test_note("   . gamma = %.1f: P = %.6f contra LZ %.6f", G[a], p, 1.0 - exp(-G[a]));
    }
    /* Resolução: com extremos a K Omega, a mistura residual nos extremos é
       ~ 1/(2K) em amplitude, e entra em P com no máximo ~ 2/K. */
    dae_test_ok(&T, pior < 2.0 / K, "2. rampa linear contra Landau-Zener (pior %.2e, resolucao 2/K = %.1e)",
                pior, 2.0 / K);
    {
      const double dD = 2.0 * PI * omega * omega * TK;          /* gamma = 1 */
      const double meia = 2.0 * omega / dD;                     /* K = 1      */
      caminho_padrao(&c); c.dye_ds = dD / (1.5 * VR);
      roda(&P, &c, -meia, meia, 50, 1e-10, &p, &nrm);
      dae_test_ok(&T, fabs(p - (1.0 - exp(-1.0))) > 0.05,
                  "2. companheira: extremos a 2 Omega divergem de LZ (P = %.4f contra %.4f)",
                  p, 1.0 - exp(-1.0));
    }
  }

  /* ------------------------------------------- 3. mu = 0 e mu != 0 ---- */
  {
    const double L = 0.5 * PI / omega / TK;
    P.mu_re[0] = 0.0;
    caminho_padrao(&c);
    roda(&P, &c, 0.0, L, 1, 1e-12, &p, &nrm);
    dae_test_ok(&T, p == 0.0, "3. mu = 0: P = 0 exatamente (%.2e)", p);
    P.mu_re[0] = 1.0e-12;
    roda(&P, &c, 0.0, L, 1, 1e-12, &p, &nrm);
    dae_test_ok(&T, p > 0.99, "3. companheira: mu != 0 converte (P = %.4f)", p);
  }

  /* ---------------------------------------- 4. B paralelo a k ---- */
  {
    const double L = 0.5 * PI / omega / TK;
    caminho_padrao(&c); c.b_perp = 0.0; c.b_par = 1.0e15;
    roda(&P, &c, 0.0, L, 1, 1e-12, &p, &nrm);
    dae_test_ok(&T, p < 1e-20, "4. B paralelo a k: P = %.2e", p);
    caminho_padrao(&c);
    roda(&P, &c, 0.0, L, 1, 1e-12, &p, &nrm);
    dae_test_ok(&T, p > 0.99, "4. companheira: B perpendicular converte (P = %.4f)", p);
  }

  /* --------------------------- 5. hermiticidade, traço, norma ---- */
  {
    double hre[4], him[4];
    dae_numag_meio m;
    caminho_padrao(&c); c.ye0 = 0.1; c.phi0 = 0.7;
    meio_fn(&c, 0.0, &m);
    st = dae_numag_hamiltoniano(&P, &m, hre, him);
    dae_test_ok(&T, st == DAE_OK && hre[1] == hre[2] && him[1] == -him[2] &&
                    him[0] == 0.0 && him[3] == 0.0,
                "5. H hermitiano");
    dae_test_ok(&T, hre[0] + hre[3] == 0.0, "5. traco nulo no espaco inteiro");
    dae_test_near(&T, hre[0] - hre[3], dae_numag_desajuste(&m), 1e-15,
                  "5. H_LL - H_RR = Delta = V rho (3 Ye - 1)/2");
    dae_test_ok(&T, pior_norma < 1e-10, "5. norma preservada em todas as corridas (%.2e)",
                pior_norma);
  }

  /* ------------------------------ 6. posição da ressonância ---- */
  {
    const double delta = DAE_NU_V_POR_RHO * RHO * (3.0 * 0.1 - 1.0) * 0.5;
    const double amp = 4.0 * omega * omega / (delta * delta + 4.0 * omega * omega);
    const double L = PI / sqrt(delta * delta + 4.0 * omega * omega) / TK;
    caminho_padrao(&c); c.ye0 = 0.1;
    roda(&P, &c, 0.0, L, 1, 1e-14, &p, &nrm);
    dae_test_near(&T, p / amp, 1.0, 1e-6,
                  "6. Ye = 0,1: P no primeiro maximo = %.4e contra a supressao %.4e", p, amp);
    dae_test_near(&T, amp, 1.91e-11, 0.01e-11, "6. a amplitude maxima analitica");
  }

  /* --------------------------- 7. campo girante e o sinal ---- */
  {
    const double dphi_nat = 2.0 * omega;                 /* phi' = 2 Omega */
    const double dphi = dphi_nat * TK;                   /* rad/km         */
    static const double X[4] = { 0.25, 0.5, 0.9, 1.7 };
    double pior_menos = 0.0, pior_mais = 0.0, p_ress, p_mais_max = 0.0;

    for (b = 0; b < 4; ++b) {
      const double L = X[b] * PI / omega / TK;
      /* Delta = -phi'/T: Delta_ef = 0 */
      caminho_padrao(&c); c.ye0 = ye_para(-dphi_nat); c.dphi_ds = dphi;
      st = roda(&P, &c, 0.0, L, 20, 1e-12, &p, &nrm);
      if (st != DAE_OK) { dae_test_ok(&T, 0, "girante: %s", dae_strerror(st)); return dae_test_end(&T); }
      if (fabs(p - exata(0.0, omega, L * TK)) > pior_menos) pior_menos = fabs(p - exata(0.0, omega, L * TK));
      /* Delta = +phi'/T: Delta_ef = 2 phi' = 4 Omega */
      caminho_padrao(&c); c.ye0 = ye_para(+dphi_nat); c.dphi_ds = dphi;
      roda(&P, &c, 0.0, L, 20, 1e-12, &p, &nrm);
      if (fabs(p - exata(2.0 * dphi_nat, omega, L * TK)) > pior_mais)
        pior_mais = fabs(p - exata(2.0 * dphi_nat, omega, L * TK));
      if (p > p_mais_max) p_mais_max = p;
    }
    dae_test_ok(&T, pior_menos < 1e-9,
                "7. Delta = -phi'/T segue a camada com Delta_ef = 0 (pior %.2e)", pior_menos);
    dae_test_ok(&T, pior_mais < 1e-9,
                "7. Delta = +phi'/T segue a camada com Delta_ef = 2 phi'/T (pior %.2e)", pior_mais);

    /* O sinal, dito sem fórmula: ressonância completa onde Delta != 0, e
       supressão no sinal oposto. */
    caminho_padrao(&c); c.ye0 = ye_para(-dphi_nat); c.dphi_ds = dphi;
    roda(&P, &c, 0.0, 0.5 * PI / omega / TK, 20, 1e-12, &p_ress, &nrm);
    dae_test_near(&T, p_ress, 1.0, 1e-9,
                  "7. SINAL: ressonancia completa em Delta = -phi'/T = %.3e eV (Delta != 0)",
                  -dphi_nat);
    dae_test_ok(&T, p_mais_max < 0.21,
                "7. SINAL: em Delta = +phi'/T a amplitude fica em 0,2 (max %.4f)", p_mais_max);

    /* Companheira: phi constante nao nulo = phi nulo. */
    {
      double p0, p1;
      const double L = 0.9 * PI / omega / TK;
      caminho_padrao(&c); c.ye0 = ye_para(omega);
      roda(&P, &c, 0.0, L, 1, 1e-12, &p0, &nrm);
      c.phi0 = 1.234;
      roda(&P, &c, 0.0, L, 1, 1e-12, &p1, &nrm);
      dae_test_ok(&T, fabs(p1 - p0) < 1e-12,
                  "7. companheira: phi_B constante nao muda P (%.2e)", fabs(p1 - p0));
    }
  }

  /* ------------------------------------ escopo do modulo ---- */
  {
    double hre[4], him[4];
    dae_numag_meio m;
    dae_numag_params Q = P;
    caminho_padrao(&c); meio_fn(&c, 0.0, &m);
    Q.sabor.antineutrino = 1;
    dae_test_ok(&T, dae_numag_hamiltoniano(&Q, &m, hre, him) == DAE_ERR_UNIMPLEMENTED,
                "escopo: antineutrino recusado");
    Q = P; Q.natureza = DAE_NUMAG_MAJORANA;
    dae_test_ok(&T, dae_numag_hamiltoniano(&Q, &m, hre, him) == DAE_ERR_UNIMPLEMENTED,
                "escopo: Majorana recusado");
    Q = P; Q.n_ativos = 3;
    dae_test_ok(&T, dae_numag_hamiltoniano(&Q, &m, hre, him) == DAE_ERR_UNIMPLEMENTED,
                "escopo: tres sabores recusados");
  }

  return dae_test_end(&T);
}
