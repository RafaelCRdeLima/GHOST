/* t04_complete.c — ACEITAÇÃO 4: K_N reduz a dois níveis.
 *
 * A = J - I, com J o operador de todos-uns. |s> = (1/sqrt N) sum_j |j> tem
 * H|s> = -gamma(N-1)|s>, e todo vetor ortogonal a |s> tem H|v> = +gamma|v>.
 * Partindo de |0>:
 *
 *   <j|psi(t)> = (1/N) exp(+i gamma (N-1) t) + exp(-i gamma t) (delta_j0 - 1/N)
 *
 * Também exercita o aperto de Lanczos: Gershgorin dá [-49, 49] em N=50, e o
 * espectro verdadeiro é {-49, +1}. É o caso em que a cota rigorosa é frouxa
 * por um fator 2 e o refinamento vale a pena.
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
  dae_rng g;
  const int32_t n = 50;
  const double gamma = 1.0, t = 1.7;
  double *re, *im, worst = 0.0;
  int32_t i;
  dae_status st;

  dae_test_begin(&T, "aceitacao 4: grafo completo K_N");

  dae_gen_params_default(&p);
  p.kind = DAE_G_COMPLETE;
  p.n = n;
  st = dae_graph_build(&G, &p);
  dae_test_ok(&T, st == DAE_OK, "K_N: %s", dae_strerror(st));
  dae_test_ok(&T, G.A.nnz == n * (n - 1), "nnz = N(N-1): %d", (int)G.A.nnz);

  st = dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, gamma, DAE_NORM_NONE, 0, NULL);
  dae_test_ok(&T, st == DAE_OK, "hamiltoniano: %s", dae_strerror(st));

  dae_rng_seed(&g, 4242ULL);
  st = dae_cheb_init(&W, &H, 20);
  dae_test_ok(&T, st == DAE_OK, "cheb_init: %s", dae_strerror(st));
  dae_test_ok(&T, W.lanczos_used, "Lanczos apertou Gershgorin");
  dae_test_ok(&T, W.lo <= -gamma * (double)(n - 1) + 1e-9 && W.hi >= gamma - 1e-9,
              "intervalo [%.6g, %.6g] contem o espectro {%g, %g}",
              W.lo, W.hi, -gamma * (double)(n - 1), gamma);
  dae_test_ok(&T, W.a < 0.75 * (gamma * (double)(n - 1)),
              "aperto efetivo: a = %.6g contra %.6g de Gershgorin",
              W.a, gamma * (double)(n - 1));

  re = (double *)calloc((size_t)n, sizeof(double));
  im = (double *)calloc((size_t)n, sizeof(double));
  if (!re || !im) { free(re); free(im); return 1; }
  re[0] = 1.0;

  st = dae_cheb_step(&W, &H, t, re, im, NULL);
  dae_test_ok(&T, st == DAE_OK, "passo: %s", dae_strerror(st));

  for (i = 0; i < n; ++i) {
    const double a1 = gamma * (double)(n - 1) * t;   /* fase do setor simetrico */
    const double a2 = -gamma * t;                    /* fase do complemento     */
    const double dj = (i == 0) ? 1.0 : 0.0;
    const double wr = cos(a1) / (double)n + cos(a2) * (dj - 1.0 / (double)n);
    const double wi = sin(a1) / (double)n + sin(a2) * (dj - 1.0 / (double)n);
    if (fabs(re[i] - wr) > worst) worst = fabs(re[i] - wr);
    if (fabs(im[i] - wi) > worst) worst = fabs(im[i] - wi);
  }
  dae_test_ok(&T, worst < 1e-10, "pior desvio da amplitude = %.3g", worst);
  dae_test_note("pior desvio = %.3g (limite 1e-10); a=%.4g contra %.4g de Gershgorin", worst, W.a, 49.0);

  free(re); free(im);
  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
  return dae_test_end(&T);
}
