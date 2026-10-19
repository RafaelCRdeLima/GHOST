/* t05_cheb_vs_exact.c — ACEITAÇÃO 5: Chebyshev contra diagonalização exata.
 *
 * Grafo aleatório com N = 200, H denso montado a partir da CSR, autossistema
 * por Jacobi (oráculo O(N^3), fora do núcleo) e soma espectral
 * psi(t) = V diag(exp(-i E_k t)) V^T psi(0). Erro exigido: < 1e-10 no vetor.
 */
#include "harness.h"
#include "jacobi.h"

#include <math.h>
#include <stdlib.h>

static int compare(dae_test *T, dae_ham_kind kind, dae_norm_kind norm,
                   const char *label)
{
  const int32_t n = 200;
  const double gamma = 1.0, t = 5.0;
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  dae_rng g;
  double *dense, *vec, *val, *re, *im, *xr, *xi, worst = 0.0, scale = 0.0;
  int32_t i, j, k, p;
  dae_status st;

  st = dae_test_random_graph(&G, n, 6, 31337ULL);
  dae_test_ok(T, st == DAE_OK, "grafo: %s", dae_strerror(st));
  if (st != DAE_OK) return 1;

  dae_rng_seed(&g, 99ULL);
  st = dae_hamiltonian(&H, &G.A, kind, gamma, norm, 40, &scale);
  dae_test_ok(T, st == DAE_OK, "hamiltoniano: %s", dae_strerror(st));

  dense = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
  vec   = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
  val   = (double *)malloc((size_t)n * sizeof(double));
  re = (double *)calloc((size_t)n, sizeof(double));
  im = (double *)calloc((size_t)n, sizeof(double));
  xr = (double *)calloc((size_t)n, sizeof(double));
  xi = (double *)calloc((size_t)n, sizeof(double));
  if (!dense || !vec || !val || !re || !im || !xr || !xi) return 1;

  for (i = 0; i < n; ++i)
    for (p = H.rowptr[i]; p < H.rowptr[i + 1]; ++p)
      dense[(size_t)i * (size_t)n + (size_t)H.colind[p]] = H.val[p];

  dae_test_ok(T, dae_jacobi(dense, vec, val, n) == 0, "Jacobi convergiu");

  /* estado inicial deslocalizado: exercita todos os autovetores */
  {
    double nrm = 0.0;
    dae_rng gg;
    dae_rng_seed(&gg, 5ULL);
    for (i = 0; i < n; ++i) {
      re[i] = 2.0 * dae_rng_uniform(&gg) - 1.0;
      im[i] = 2.0 * dae_rng_uniform(&gg) - 1.0;
      nrm += re[i] * re[i] + im[i] * im[i];
    }
    nrm = sqrt(nrm);
    for (i = 0; i < n; ++i) { re[i] /= nrm; im[i] /= nrm; }
  }

  /* referência: V exp(-i E t) V^T psi0 */
  for (k = 0; k < n; ++k) {
    double cr = 0.0, ci = 0.0;
    const double ph = -val[k] * t;
    for (j = 0; j < n; ++j) {
      cr += vec[(size_t)j * (size_t)n + (size_t)k] * re[j];
      ci += vec[(size_t)j * (size_t)n + (size_t)k] * im[j];
    }
    { const double ar = cr * cos(ph) - ci * sin(ph);
      const double ai = cr * sin(ph) + ci * cos(ph);
      for (i = 0; i < n; ++i) {
        const double v = vec[(size_t)i * (size_t)n + (size_t)k];
        xr[i] += v * ar;
        xi[i] += v * ai;
      } }
  }

  st = dae_cheb_init(&W, &H, 40);
  dae_test_ok(T, st == DAE_OK, "cheb_init: %s", dae_strerror(st));
  st = dae_cheb_step(&W, &H, t, re, im, NULL);
  dae_test_ok(T, st == DAE_OK, "passo: %s", dae_strerror(st));

  for (i = 0; i < n; ++i) {
    if (fabs(re[i] - xr[i]) > worst) worst = fabs(re[i] - xr[i]);
    if (fabs(im[i] - xi[i]) > worst) worst = fabs(im[i] - xi[i]);
  }
  dae_test_ok(T, worst < 1e-10, "%s: pior desvio = %.3g (escala aplicada %.6g)",
              label, worst, scale);
  dae_test_note("%s: pior desvio = %.3g (limite 1e-10), escala %.6g", label, worst, scale);

  free(dense); free(vec); free(val); free(re); free(im); free(xr); free(xi);
  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
  return 0;
}

int main(void)
{
  dae_test T;
  dae_test_begin(&T, "aceitacao 5: Chebyshev contra diagonalizacao exata");
  compare(&T, DAE_H_ADJACENCY, DAE_NORM_NONE,        "adjacencia, sem normalizar");
  compare(&T, DAE_H_LAPLACIAN, DAE_NORM_NONE,        "laplaciana, sem normalizar");
  compare(&T, DAE_H_ADJACENCY, DAE_NORM_SPECTRAL,    "adjacencia, raio espectral");
  compare(&T, DAE_H_LAPLACIAN, DAE_NORM_MEAN_DEGREE, "laplaciana, grau medio");
  return dae_test_end(&T);
}
