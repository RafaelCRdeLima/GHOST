/* t106_dois_estados_edgelist.c — ACEITAÇÃO 106: o modelo magnético mínimo de
 * dois estados (nu_L, nu_R de um sabor Dirac) roda sem módulo próprio, como
 * edgelist pelo caminho completo spec -> dae_run.
 *
 * O Hamiltoniano sem traço
 *
 *     H = [ Delta/2   Omega  ]
 *         [ Omega   -Delta/2 ]
 *
 * entra como laços [0,0,-Delta/2], [1,1,+Delta/2] e aresta [0,1,-Omega],
 * porque o pipeline monta H = -gamma*A com gamma = 1. A resposta exata numa
 * camada constante é
 *
 *     P_{L->R}(t) = 4 Omega^2/(Delta^2 + 4 Omega^2)
 *                   * sin^2( t/2 * sqrt(Delta^2 + 4 Omega^2) ).
 *
 * O QUE ISTO VALIDA E O QUE NÃO VALIDA. Valida a forma da precessão e a
 * escala de tempo, inclusive em unidades físicas (Omega em eV, L em km via
 * DAE_NU_T_POR_KM). NÃO valida sinais: P é par em Delta e em Omega, então
 * um sinal trocado passaria aqui. O sinal do acoplamento só fica testável
 * com o campo girante (Delta_ef = Delta -+ dphi_B/ds).
 *
 * ANTI-VACUIDADE. Omega = 0 tem de dar P = 0 exatamente, e Delta = 0 tem de
 * atingir P = 1 no primeiro máximo: um teste que só olhasse P pequeno passaria
 * com um propagador que não faz nada.
 */
#include "harness.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

/* mu_B em eV/T; B = 1e15 G = 1e11 T; mu_nu = 1e-12 mu_B. */
#define MU_B_EV_POR_T 5.7883818060e-5

static double exata(double delta, double omega, double t)
{
  const double w = sqrt(delta * delta + 4.0 * omega * omega);
  if (w == 0.0) return 0.0;
  { const double s = sin(0.5 * t * w);
    return 4.0 * omega * omega / (w * w) * s * s; }
}

/* Roda o spec edgelist e devolve P_R = |psi_1|^2 no instante final. */
static int roda(dae_test *T, double delta, double omega, double t1, int32_t nt,
                double *p_r, double *norma)
{
  char texto[1024];
  dae_spec S;
  dae_graph G;
  dae_series R;
  dae_error err;
  dae_status st;

  snprintf(texto, sizeof texto,
    "{\"format_version\":3,\"seed\":1,"
    "\"graph\":{\"generator\":\"edgelist\",\"n\":2,"
    "\"edges\":[[0,0,%.17g],[1,1,%.17g],[0,1,%.17g]]},"
    "\"hamiltonian\":{\"kind\":\"adjacency\",\"gamma\":1.0,\"normalization\":\"none\"},"
    "\"initial\":{\"site\":0},"
    "\"time\":{\"t1\":%.17g,\"nt\":%d},"
    "\"observables\":{\"target\":1,\"population\":false},\"realizations\":1}",
    -0.5 * delta, 0.5 * delta, -omega, t1, (int)nt);

  st = dae_spec_parse(&S, texto, &err);
  if (st != DAE_OK) { dae_test_ok(T, 0, "spec: %s", err.msg); return 1; }
  st = dae_run(&S, S.seed, &G, NULL, &R, NULL, NULL);
  if (st != DAE_OK) {
    dae_test_ok(T, 0, "run: %s", dae_strerror(st));
    dae_spec_free(&S);
    return 1;
  }
  *p_r = R.psi_re[1] * R.psi_re[1] + R.psi_im[1] * R.psi_im[1];
  *norma = *p_r + R.psi_re[0] * R.psi_re[0] + R.psi_im[0] * R.psi_im[0];
  dae_series_free(&R); dae_graph_free(&G); dae_spec_free(&S);
  return 0;
}

int main(void)
{
  dae_test T;
  /* q = Delta/(2 Omega). */
  static const double Q[3] = { 0.0, 1.0, 5.0 };
  /* comprimentos em unidades de pi/Omega */
  static const double X[5] = { 0.25, 0.5, 0.9, 1.7, 2.6 };
  double pior = 0.0, p, nrm, pior_norma = 0.0;
  int32_t a, b;

  dae_test_begin(&T, "t106 dois estados nu_L/nu_R por edgelist contra a formula exata");

  /* 1. Unidades adimensionais, Omega = 1. */
  for (a = 0; a < 3; ++a) {
    const double omega = 1.0, delta = 2.0 * omega * Q[a];
    for (b = 0; b < 5; ++b) {
      const double t = X[b] * PI / omega;
      if (roda(&T, delta, omega, t, 50, &p, &nrm)) return dae_test_end(&T);
      { const double d = fabs(p - exata(delta, omega, t));
        if (d > pior) pior = d; }
      if (fabs(nrm - 1.0) > pior_norma) pior_norma = fabs(nrm - 1.0);
    }
  }
  dae_test_ok(&T, pior < 1e-13,
              "P_{L->R} contra a formula exata, q = 0, 1, 5 (pior %.3e)", pior);
  dae_test_ok(&T, pior_norma < 1e-13, "norma preservada (pior %.3e)", pior_norma);
  dae_test_note("   . resolucao: 15 pontos (3 valores de q x 5 comprimentos), 50 passos cada");

  /* 2. Anti-vacuidade: Omega = 0 nao converte; Delta = 0 converte tudo. */
  if (roda(&T, 2.0, 0.0, 3.7, 50, &p, &nrm)) return dae_test_end(&T);
  dae_test_ok(&T, p == 0.0, "Omega = 0 da P = 0 exatamente (%.3e)", p);
  if (roda(&T, 0.0, 1.0, 0.5 * PI, 50, &p, &nrm)) return dae_test_end(&T);
  dae_test_near(&T, p, 1.0, 1e-13, "Delta = 0 atinge P = 1 em t = pi/(2 Omega)");

  /* 3. Unidades fisicas: Omega = mu_nu B em eV, comprimento em km.
        Na ressonancia o primeiro maximo esta em L = pi hbar c/(2 Omega),
        ~5,4 cm para mu_nu = 1e-12 mu_B e B_perp = 1e15 G. */
  {
    const double omega = 1e-12 * MU_B_EV_POR_T * 1e11;          /* eV */
    const double l_max_km = PI / (2.0 * omega) / DAE_NU_T_POR_KM;
    if (roda(&T, 0.0, omega, l_max_km * DAE_NU_T_POR_KM, 50, &p, &nrm))
      return dae_test_end(&T);
    dae_test_near(&T, p, 1.0, 1e-12,
                  "unidades fisicas: P = 1 em L = %.4f cm", l_max_km * 1e5);
    dae_test_near(&T, l_max_km * 1e5, 5.356, 5e-3,
                  "o primeiro maximo fica a ~5,36 cm (metade de L_per)");

    /* Fora da ressonancia: rho = 1e14, Y_e = 0,1. */
    {
      const double delta = DAE_NU_V_POR_RHO * 0.5 * 1.0e14 * (3.0 * 0.1 - 1.0);
      const double amp = 4.0 * omega * omega / (delta * delta + 4.0 * omega * omega);
      const double w = sqrt(delta * delta + 4.0 * omega * omega);
      const double t = PI / w;                    /* primeiro maximo */
      if (roda(&T, delta, omega, t, 50, &p, &nrm)) return dae_test_end(&T);
      dae_test_ok(&T, fabs(p - amp) < 1e-13 * 1.0 + 1e-6 * amp,
                  "supressao pela materia: P_max = %.4e contra %.4e", p, amp);
      dae_test_near(&T, amp, 1.91e-11, 0.01e-11,
                    "amplitude maxima analitica (1,91e-11)");
    }
  }

  return dae_test_end(&T);
}
