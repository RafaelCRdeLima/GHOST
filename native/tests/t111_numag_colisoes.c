/* t111_numag_colisoes.c — ACEITAÇÃO 111: colisões como defasagem, por
 * trajetórias. O modelo está em dae_numag.h, seção "colisões e coerência".
 *
 *  1. O FATOR. gamma [1/km] = Gamma [eV] T / 2: a taxa de salto do
 *     desdobramento simétrico é METADE da taxa de colisão, porque o operador
 *     físico sqrt(Gamma)|L><L| amortece a coerência a Gamma/2 e o simétrico
 *     sqrt(gamma)P_j a gamma.
 *  2. H = 0, partindo de (L + R)/sqrt(2): a coerência média tem de seguir
 *     (1/2) exp(-gamma s), dentro de 4 sigma em 4 pontos, e P_R ficar em 1/2.
 *     É a verificação direta do fator: com gamma = Gamma T (sem o 1/2) a
 *     coerência cairia como exp(-2 gamma s) e o teste falharia por dezenas de
 *     sigma. Companheira: Gamma = 0 dá coerência 1/2 com variância ZERO.
 *  3. REPRODUTIBILIDADE: mesma semente, mesmo resultado bit a bit; semente
 *     diferente, resultado diferente.
 *  4. CRUZAMENTO SEM COLISÃO: o ensemble com Gamma = 0 é a propagação coerente,
 *     sem variância.
 *  5. DEFASAGEM FORTE (Gamma/Omega = 40): o limite das equações de taxa,
 *     P = (1 - exp(-2 gamma_LZ f))/2 com f = (2/pi) arctan(2K / (Gamma/2Omega))
 *     (o domínio finito incluído), dentro de 4 sigma mais a correção de modelo
 *     medida pela referência de Liouville (~1e-3 neste K).
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define TK DAE_NU_T_POR_KM
#define RHO 1.0e14
#define VR (DAE_NU_V_POR_RHO * RHO)

int main(void)
{
  dae_test T;
  dae_numag_params P;
  const double om = 1e-12 * DAE_NUMAG_MU_B_EV_POR_T * 1e15 * DAE_NUMAG_T_POR_GAUSS;

  dae_test_begin(&T, "t111 colisoes por trajetorias: o fator Gamma/2, Liouville e os limites");
  dae_numag_params_default(&P);
  P.mu_re[0] = 1e-12;

  /* ---------------------------------------------------------- 1. fator -- */
  dae_test_ok(&T, dae_numag_salto_por_km(2.0e-6) == 0.5 * 2.0e-6 * TK,
              "1. gamma [1/km] = Gamma [eV] T / 2 (%.6e)", dae_numag_salto_por_km(2.0e-6));

  /* ------------------------------------------------- 2. H = 0, coerência -- */
  {
    const double L = 1e-3;
    double s[2] = { 0.0, L }, r[2] = { RHO, RHO }, y[2] = { 1.0 / 3.0, 1.0 / 3.0 };
    double b1[2] = { 1e15, 1e15 }, b2[2] = { 0.0, 0.0 };
    dae_numag_trajetoria Tr;
    dae_numag_params P0 = P;
    dae_numag_colisao c;
    const double pr0[2] = { 1.0 / sqrt(2.0), 1.0 / sqrt(2.0) }, pi0[2] = { 0.0, 0.0 };
    const int32_t N = 4000, nt = 4;
    double pR[4], var[4], cre[4], cim[4], sal;
    double Gamma;
    int k, ok_coh = 1, ok_pop = 1;
    double pior = 0.0;
    dae_status st;
    Tr.s_km = s; Tr.rho = r; Tr.ye = y; Tr.b1 = b1; Tr.b2 = b2; Tr.razao_energia = NULL; Tr.ynu = NULL; Tr.n = 2;
    P0.mu_re[0] = 0.0;                          /* H = 0 */
    Gamma = 2.0 * 2.0 / (TK * L);               /* gamma L = 2 */
    c.gamma0_km = dae_numag_salto_por_km(Gamma); c.gamma_km = NULL;
    st = dae_numag_ensemble(&P0, &Tr, 0.0, L, &c, 1e-10, N, 7, pr0, pi0, nt,
                            pR, var, cre, cim, &sal);
    dae_test_ok(&T, st == DAE_OK, "2. ensemble com H = 0 (%s)", dae_strerror(st));
    for (k = 0; k < nt; ++k) {
      const double sk = L * (k + 1) / nt, ex = 0.5 * exp(-c.gamma0_km * sk);
      /* |psi_L psi_R^*| <= 1/2: o desvio padrão da média é <= 0,5/sqrt(N) */
      const double sig = 0.5 / sqrt((double)N);
      const double z = fabs(cre[k] - ex) / sig;
      if (z > pior) pior = z;
      if (z > 4.0 || fabs(cim[k]) > 4.0 * sig) ok_coh = 0;
      if (fabs(pR[k] - 0.5) > 4.0 * sqrt(var[k] / N) + 1e-12) ok_pop = 0;
    }
    dae_test_ok(&T, ok_coh,
                "2. coerencia media = (1/2) exp(-gamma s), gamma = Gamma T/2, em 4 pontos (pior %.1f sigma)", pior);
    dae_test_ok(&T, ok_pop, "2. populacoes ficam em 1/2 (a defasagem nao as move)");
    dae_test_ok(&T, fabs(sal - c.gamma0_km * L) < 4.0 * sqrt(c.gamma0_km * L / N),
                "2. saltos por trajetoria = gamma L = %.2f (medido %.3f)", c.gamma0_km * L, sal);
    /* Sem o fator 1/2, a coerência no último ponto seria e^{-4}/2 = 0,009: */
    dae_test_note("   . com gamma = Gamma T (sem 1/2) o ultimo ponto seria %.4f; medido %.4f, esperado %.4f",
                  0.5 * exp(-2.0 * c.gamma0_km * L), cre[nt - 1], 0.5 * exp(-c.gamma0_km * L));

    c.gamma0_km = 0.0;
    dae_numag_ensemble(&P0, &Tr, 0.0, L, &c, 1e-10, 50, 7, pr0, pi0, nt, pR, var, cre, cim, &sal);
    /* a variância é acumulada por Welford: valores iguais dão zero EXATO */
    dae_test_ok(&T, fabs(cre[nt - 1] - 0.5) < 1e-14 && var[nt - 1] == 0.0 && sal == 0.0,
                "2. companheira: Gamma = 0 deixa a coerencia em 1/2 (%.17g), sem saltos, "
                "variancia zero (%.1e)", cre[nt - 1], var[nt - 1]);
  }

  /* ----------------------------------- 3, 4, 5: um cruzamento linear ---- */
  {
    const double g = 1.0, K = 30.0;
    const double dD = 2.0 * PI * om * om * TK / g;       /* eV/km */
    const double meia = K * 2.0 * om / dD;
    double s[3], r[3] = { RHO, RHO, RHO }, y[3], b1[3] = { 1e15, 1e15, 1e15 }, b2[3] = { 0, 0, 0 };
    dae_numag_trajetoria Tr;
    dae_numag_colisao c;
    const double pr0[2] = { 1.0, 0.0 }, pi0[2] = { 0.0, 0.0 };
    double a, b, va, vb, sal;
    dae_status st;
    s[0] = 0.0; s[1] = meia; s[2] = 2.0 * meia;
    y[0] = 1.0 / 3.0 - dD * meia / (1.5 * VR); y[1] = 1.0 / 3.0; y[2] = 1.0 / 3.0 + dD * meia / (1.5 * VR);
    Tr.s_km = s; Tr.rho = r; Tr.ye = y; Tr.b1 = b1; Tr.b2 = b2; Tr.razao_energia = NULL; Tr.ynu = NULL; Tr.n = 3;
    c.gamma_km = NULL;

    /* 4. Gamma = 0 = coerente */
    {
      double re[4], im[4], pr[2], pi[2], pc;
      c.gamma0_km = 0.0;
      st = dae_numag_ensemble(&P, &Tr, 0.0, 2.0 * meia, &c, 1e-9, 3, 1, pr0, pi0, 1, &a, &va, NULL, NULL, &sal);
      dae_nu_inicial(0, 2, 1, re, im);
      dae_numag_propaga_tabela(&P, &Tr, 0.0, 2.0 * meia, 1e-9, re, im, NULL);
      dae_nu_projeta(re, im, 2, 1, pr, pi);
      pc = pr[1] * pr[1] + pi[1] * pi[1];
      dae_test_ok(&T, st == DAE_OK && a == pc && va == 0.0,
                  "4. Gamma = 0: o ensemble E a propagacao coerente (%.10f, variancia %.1e)", a, va);
    }

    /* 3. reprodutibilidade */
    c.gamma0_km = dae_numag_salto_por_km(1.0 * om);
    dae_numag_ensemble(&P, &Tr, 0.0, 2.0 * meia, &c, 1e-8, 20, 42, pr0, pi0, 1, &a, &va, NULL, NULL, &sal);
    dae_numag_ensemble(&P, &Tr, 0.0, 2.0 * meia, &c, 1e-8, 20, 42, pr0, pi0, 1, &b, &vb, NULL, NULL, &sal);
    dae_test_ok(&T, a == b && va == vb, "3. mesma semente: identico bit a bit (%.12f)", a);
    dae_numag_ensemble(&P, &Tr, 0.0, 2.0 * meia, &c, 1e-8, 20, 43, pr0, pi0, 1, &b, &vb, NULL, NULL, &sal);
    dae_test_ok(&T, a != b, "3. companheira: outra semente, outro resultado (%.6f contra %.6f)", a, b);

    /* 5. defasagem forte */
    {
      const double GO = 40.0;
      const int32_t N = 200;
      const double f = 2.0 / PI * atan(2.0 * K / (0.5 * GO));
      const double forte = 0.5 * (1.0 - exp(-2.0 * g * f));
      c.gamma0_km = dae_numag_salto_por_km(GO * om);
      st = dae_numag_ensemble(&P, &Tr, 0.0, 2.0 * meia, &c, 1e-7, N, 5, pr0, pi0, 1, &a, &va, NULL, NULL, &sal);
      dae_test_ok(&T, st == DAE_OK && fabs(a - forte) < 4.0 * sqrt(va / N) + 2e-3,
                  "5. Gamma/Omega = 40: P_R = %.4f +- %.4f contra equacoes de taxa %.4f "
                  "(coerente seria %.4f; %.0f saltos por trajetoria)",
                  a, sqrt(va / N), forte, 1.0 - exp(-g), sal);
    }
  }

  return dae_test_end(&T);
}
