/* t02_line_bessel.c — ACEITAÇÃO 2: linha infinita contra Bessel.
 *
 * Com H = -gamma A e o pacote longe das bordas, a solução exata da cadeia
 * infinita é psi_j(t) = i^d J_d(2 gamma t), d = |j - j0|. O teste compara a
 * amplitude COMPLEXA, não só a população: |psi|^2 sozinho perderia um erro de
 * fase, que é justo o tipo de erro que o fator (-i)^k dos coeficientes de
 * Chebyshev pode introduzir.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

static void ipow(int d, double *re, double *im)
{
  switch (d & 3) {
    case 0:  *re = 1.0;  *im = 0.0;  break;
    case 1:  *re = 0.0;  *im = 1.0;  break;
    case 2:  *re = -1.0; *im = 0.0;  break;
    default: *re = 0.0;  *im = -1.0; break;
  }
}

static int run(dae_test *T, int32_t nsteps, const char *label)
{
  dae_gen_params p;
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  const int32_t n = 2001, j0 = 1000, dmax = 200;
  const double gamma = 1.0, tfinal = 20.0;
  double *re, *im, *jref, worst = 0.0;
  int32_t i, s;
  dae_status st;

  dae_gen_params_default(&p);
  p.kind = DAE_G_PATH;
  p.n = n;
  st = dae_graph_build(&G, &p);
  dae_test_ok(T, st == DAE_OK, "linha: %s", dae_strerror(st));
  if (st != DAE_OK) return 1;

  st = dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, gamma, DAE_NORM_NONE, 0, NULL);
  dae_test_ok(T, st == DAE_OK, "hamiltoniano: %s", dae_strerror(st));
  st = dae_cheb_init(&W, &H, 0);          /* Gershgorin puro: [-2, 2] */
  dae_test_ok(T, st == DAE_OK, "cheb_init: %s", dae_strerror(st));
  dae_test_near(T, W.a, 2.0 * gamma, 1e-15, "semi-largura de Gershgorin: ");

  re = (double *)calloc((size_t)n, sizeof(double));
  im = (double *)calloc((size_t)n, sizeof(double));
  jref = (double *)malloc((size_t)(dmax + 1) * sizeof(double));
  if (!re || !im || !jref) { free(re); free(im); free(jref); return 1; }
  re[j0] = 1.0;

  for (s = 0; s < nsteps; ++s) {
    st = dae_cheb_step(&W, &H, tfinal / (double)nsteps, re, im, NULL);
    if (st != DAE_OK) { dae_test_ok(T, 0, "passo: %s", dae_strerror(st)); break; }
  }

  dae_bessel_j_array(2.0 * gamma * tfinal, dmax, jref);
  for (i = 0; i <= dmax; ++i) {
    double pr, pi;
    ipow((int)i, &pr, &pi);
    { const double d1 = fabs(re[j0 + i] - pr * jref[i]);
      const double d2 = fabs(im[j0 + i] - pi * jref[i]);
      if (d1 > worst) worst = d1;
      if (d2 > worst) worst = d2; }
    /* a linha é simétrica em torno de j0 */
    { const double d3 = fabs(re[j0 - i] - pr * jref[i]);
      if (d3 > worst) worst = d3; }
  }
  dae_test_ok(T, worst < 1e-10, "%s: pior desvio da amplitude = %.3g", label, worst);
  dae_test_note("%s: pior desvio da amplitude = %.3g (limite 1e-10)", label, worst);

  free(re); free(im); free(jref);
  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
  return 0;
}

int main(void)
{
  dae_test T;
  dae_test_begin(&T, "aceitacao 2: linha infinita = J_j(2 gamma t)");
  run(&T, 1,   "um passo dt=20");
  run(&T, 200, "200 passos dt=0.1");
  return dae_test_end(&T);
}
