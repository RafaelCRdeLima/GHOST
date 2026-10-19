/* t01_norm.c — ACEITAÇÃO 1: desvio da norma < 1e-12 ao longo de 1e4 passos.
 *
 * A conservação tem de valer para QUALQUER limite espectral válido, não só
 * para o que der sorte: a deriva pode ser muito sensível a `a`, a ponto de uma
 * mudança de 0,4% (trocar o vetor inicial do Lanczos) alterá-la por um fator
 * 16. Um teste com um único limite mede a sorte daquele limite. Este varre
 * cinco.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

static double run_one(dae_test *T, const dae_csr *H, int32_t n, int32_t lanczos,
                      int32_t nsteps, double *a_out, int32_t *k_out)
{
  dae_cheb W;
  dae_cheb_info info;
  double *re, *im, dt, worst = 0.0;
  int32_t i, s;
  dae_status st;

  st = dae_cheb_init(&W, H, lanczos);
  dae_test_ok(T, st == DAE_OK, "cheb_init: %s", dae_strerror(st));
  if (st != DAE_OK) return -1.0;

  re = (double *)calloc((size_t)n, sizeof(double));
  im = (double *)calloc((size_t)n, sizeof(double));
  if (!re || !im) { free(re); free(im); dae_cheb_free(&W); return -1.0; }
  re[n / 3] = 1.0;

  info.k_used = 0;
  dt = 5.0 / W.a;                       /* a*dt = 5: regime de trabalho */
  for (s = 0; s < nsteps; ++s) {
    st = dae_cheb_step(&W, H, dt, re, im, &info);
    if (st != DAE_OK) { dae_test_ok(T, 0, "passo %d: %s", (int)s, dae_strerror(st)); break; }
    if ((s % 100) == 0 || s == nsteps - 1) {
      double nrm = 0.0;
      for (i = 0; i < n; ++i) nrm += re[i] * re[i] + im[i] * im[i];
      if (fabs(nrm - 1.0) > worst) worst = fabs(nrm - 1.0);
    }
  }

  *a_out = W.a;
  *k_out = info.k_used;
  free(re); free(im);
  dae_cheb_free(&W);
  return worst;
}

/* ENVOLTÓRIA DA DERIVA. No piso de arredondamento, |norma-1| cresce como
   n_passos * K * c, com c ~ 5e-18 medido. Um truncamento cedo demais quebra
   essa lei por ordens de grandeza sem quebrar nenhum outro teste: um teto
   K ~ 1.2 alpha + 20, por exemplo, corta a serie com |2 J_K| ainda em 1.9e-11
   na faixa alpha ~ 20 a 100 (4.2e-13 de deriva POR PASSO em alpha = 50,
   contra 5.5e-17 em alpha = 5). Um teste de um alpha so nao ve isso, porque
   o defeito e nao-monotonico em alpha. */
static void varre_alpha(dae_test *T)
{
  const double alphas[5] = { 1.0, 5.0, 20.0, 50.0, 120.0 };
  const double C = 5e-17;              /* envoltoria, ~10x o medido */
  const int32_t n = 256, nsteps = 1000;
  dae_graph G;
  dae_csr H;
  int ai;

  if (dae_test_random_graph(&G, n, 4, 4242ULL) != DAE_OK) return;
  if (dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL) != DAE_OK) return;

  for (ai = 0; ai < 5; ++ai) {
    dae_cheb W;
    dae_cheb_info info;
    double *re, *im, worst = 0.0, lim;
    int32_t i, s;

    if (dae_cheb_init(&W, &H, 0) != DAE_OK) return;
    re = (double *)calloc((size_t)n, sizeof(double));
    im = (double *)calloc((size_t)n, sizeof(double));
    if (!re || !im) { free(re); free(im); dae_cheb_free(&W); return; }
    re[n / 3] = 1.0;
    info.k_used = 0;
    for (s = 0; s < nsteps; ++s) {
      if (dae_cheb_step(&W, &H, alphas[ai] / W.a, re, im, &info) != DAE_OK) break;
      if ((s % 50) == 0 || s == nsteps - 1) {
        double nrm = 0.0;
        for (i = 0; i < n; ++i) nrm += re[i] * re[i] + im[i] * im[i];
        if (fabs(nrm - 1.0) > worst) worst = fabs(nrm - 1.0);
      }
    }
    lim = C * (double)nsteps * (double)info.k_used;
    dae_test_ok(T, worst < lim,
                "alpha=%g: deriva %.3g acima da envoltoria %.3g (K=%d) — "
                "sinal de truncamento cedo demais",
                alphas[ai], worst, lim, (int)info.k_used);
    dae_test_note("alpha=%5.1f  K=%3d  deriva=%9.3g  envoltoria=%9.3g",
                  alphas[ai], (int)info.k_used, worst, lim);
    free(re); free(im);
    dae_cheb_free(&W);
  }
  dae_csr_free(&H); dae_graph_free(&G);
}

int main(void)
{
  dae_test T;
  dae_graph G;
  dae_csr H;
  const int32_t n = 512, nsteps = 10000;
  const int32_t lanczos[5] = { 0, 10, 20, 30, 50 };
  double t0, pior_geral = 0.0;
  int li;
  dae_status st;

  dae_test_begin(&T, "aceitacao 1: conservacao da norma");

  st = dae_test_random_graph(&G, n, 4, 777ULL);
  dae_test_ok(&T, st == DAE_OK, "grafo aleatorio: %s", dae_strerror(st));
  if (st != DAE_OK) return dae_test_end(&T);

  st = dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  dae_test_ok(&T, st == DAE_OK, "hamiltoniano: %s", dae_strerror(st));

  t0 = dae_test_seconds();
  for (li = 0; li < 5; ++li) {
    double a = 0.0;
    int32_t k = 0;
    const double worst = run_one(&T, &H, n, lanczos[li], nsteps, &a, &k);
    if (worst > pior_geral) pior_geral = worst;
    dae_test_ok(&T, worst >= 0.0 && worst < 1e-12,
                "lanczos=%d: pior |norma-1| = %.3g em %d passos",
                (int)lanczos[li], worst, (int)nsteps);
    dae_test_note("lanczos=%2d  a=%.6g  K=%d  pior |norma-1| = %.3g",
                  (int)lanczos[li], a, (int)k, worst);
  }
  dae_test_note("pior caso entre os 5 limites: %.3g (limite 1e-12), %.1f s",
                pior_geral, dae_test_seconds() - t0);

  dae_csr_free(&H); dae_graph_free(&G);
  varre_alpha(&T);
  return dae_test_end(&T);
}
