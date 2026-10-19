/* t03_cycle.c — ACEITAÇÃO 3: ciclo C_N contra a solução fechada.
 *
 * Espectro E_k = -2 gamma cos(2 pi k / N) e
 * psi_j(t) = (1/N) sum_k exp(-i E_k t) exp(2 pi i k (j - j0) / N).
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

int main(void)
{
  dae_test T;
  dae_gen_params p;
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  const int32_t n = 64, j0 = 0;
  const double gamma = 1.0;
  const double times[3] = { 0.7, 3.7, 13.0 };
  double *re, *im;
  int32_t i, k, ti;
  dae_status st;

  dae_test_begin(&T, "aceitacao 3: ciclo C_N");

  dae_gen_params_default(&p);
  p.kind = DAE_G_CYCLE;
  p.n = n;
  st = dae_graph_build(&G, &p);
  dae_test_ok(&T, st == DAE_OK, "ciclo: %s", dae_strerror(st));
  dae_test_ok(&T, G.A.nnz == 2 * n, "|E| = N: nnz esperado %d, obtido %d",
              (int)(2 * n), (int)G.A.nnz);

  st = dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, gamma, DAE_NORM_NONE, 0, NULL);
  dae_test_ok(&T, st == DAE_OK, "hamiltoniano: %s", dae_strerror(st));

  re = (double *)calloc((size_t)n, sizeof(double));
  im = (double *)calloc((size_t)n, sizeof(double));
  if (!re || !im) { free(re); free(im); return 1; }

  for (ti = 0; ti < 3; ++ti) {
    const double t = times[ti];
    double worst = 0.0;
    st = dae_cheb_init(&W, &H, 0);
    dae_test_ok(&T, st == DAE_OK, "cheb_init: %s", dae_strerror(st));
    for (i = 0; i < n; ++i) { re[i] = 0.0; im[i] = 0.0; }
    re[j0] = 1.0;
    st = dae_cheb_step(&W, &H, t, re, im, NULL);
    dae_test_ok(&T, st == DAE_OK, "passo t=%g: %s", t, dae_strerror(st));

    for (i = 0; i < n; ++i) {
      double sr = 0.0, si = 0.0;
      for (k = 0; k < n; ++k) {
        const double q = 6.283185307179586 * (double)k / (double)n;
        const double E = -2.0 * gamma * cos(q);
        const double ph = -E * t + q * (double)(i - j0);
        sr += cos(ph);
        si += sin(ph);
      }
      sr /= (double)n; si /= (double)n;
      if (fabs(re[i] - sr) > worst) worst = fabs(re[i] - sr);
      if (fabs(im[i] - si) > worst) worst = fabs(im[i] - si);
    }
    dae_test_ok(&T, worst < 1e-10, "t=%g: pior desvio = %.3g", t, worst);
    dae_test_note("t=%g: pior desvio = %.3g (limite 1e-10)", t, worst);
    dae_cheb_free(&W);
  }

  /* Autovetor de momento: a evolução tem de ser fase pura exp(-i E_k t). */
  {
    const int32_t kmode = 7;
    const double q = 6.283185307179586 * (double)kmode / (double)n;
    const double E = -2.0 * gamma * cos(q), t = 2.3;
    double worst = 0.0;
    st = dae_cheb_init(&W, &H, 0);
    for (i = 0; i < n; ++i) {
      re[i] = cos(q * (double)i) / sqrt((double)n);
      im[i] = sin(q * (double)i) / sqrt((double)n);
    }
    dae_cheb_step(&W, &H, t, re, im, NULL);
    for (i = 0; i < n; ++i) {
      const double ph = q * (double)i - E * t;
      const double wr = cos(ph) / sqrt((double)n), wi = sin(ph) / sqrt((double)n);
      if (fabs(re[i] - wr) > worst) worst = fabs(re[i] - wr);
      if (fabs(im[i] - wi) > worst) worst = fabs(im[i] - wi);
    }
    dae_test_ok(&T, worst < 1e-12, "modo k=7: fase pura, desvio = %.3g", worst);
    dae_test_note("modo k=7: fase pura, desvio = %.3g (limite 1e-12)", worst);
    dae_cheb_free(&W);
  }

  free(re); free(im);
  dae_csr_free(&H); dae_graph_free(&G);
  return dae_test_end(&T);
}
