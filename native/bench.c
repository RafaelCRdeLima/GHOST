/* bench.c — metas de desempenho da especificação.
 *
 *   N ~  1 300, 500 pontos no tempo  ->  < 0,2 s  (slider fluido)
 *   N ~ 10 000, 500 pontos no tempo  ->  < 3 s    ("clica e espera")
 */
#include "dae.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double seconds(void) { return (double)clock() / (double)CLOCKS_PER_SEC; }

static void run(int32_t n, double dt, double budget)
{
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  dae_cheb_info info;
  dae_obs_cfg cfg;
  dae_obs O;
  const int32_t nt = 500;
  double *re, *im, t0, el;
  int32_t i, s, ne;
  dae_rng g;

  {  /* rede de grau 4, como o lattice microtubular */
    const int32_t rows = 13, cols = (n + 12) / 13;
    dae_gen_params p;
    dae_gen_params_default(&p);
    p.kind = DAE_G_GRID2D; p.rows = rows; p.cols = cols;
    if (dae_graph_build(&G, &p) != DAE_OK) { printf("falhou\n"); return; }
  }
  ne = G.A.nnz / 2;
  dae_rng_seed(&g, 1ULL);
  dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  dae_cheb_init(&W, &H, 0);

  cfg.module_of = G.module_of; cfg.nmod = G.nmod; cfg.target = G.n - 1;
  cfg.want_pop = 1; cfg.want_conc_mod = 1; cfg.want_conc_full = 0;
  dae_obs_alloc(&O, &cfg, G.n);

  re = (double *)calloc((size_t)G.n, sizeof(double));
  im = (double *)calloc((size_t)G.n, sizeof(double));
  re[0] = 1.0;

  t0 = seconds();
  for (s = 0; s < nt; ++s) {
    if (dae_cheb_step(&W, &H, dt, re, im, &info) != DAE_OK) { printf("  ABORTOU\n"); break; }
    dae_obs_eval(re, im, G.n, &cfg, &O);
  }
  el = seconds() - t0;
  i = 0; (void)i;
  if (budget > 1e8)
    printf("  N=%6d  |E|=%7d  dt=%.2f  K=%3d  %d passos  %6.3f s   (sem meta: exportar)\n",
           (int)G.n, (int)ne, dt, (int)info.k_used, (int)nt, el);
  else
    printf("  N=%6d  |E|=%7d  dt=%.2f  K=%3d  %d passos  %6.3f s   meta %.1f s  %s\n",
           (int)G.n, (int)ne, dt, (int)info.k_used, (int)nt, el, budget,
           el < budget ? "ok" : "ESTOUROU");

  free(re); free(im);
  dae_obs_free(&O); dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
}

int main(void)
{
  printf("ghost %s  core %s  — desempenho\n", DAE_VERSION, DAE_CORE_HASH);
  /* Duas grades: a fina, com dt pequeno, é limitada pelo piso K >= 20; a
     larga é o regime em que o custo real do Chebyshev aparece. */
  printf(" grade fina (dt = 0.1, janela t = 50):\n");
  run(1300,  0.1, 0.2);
  run(10000, 0.1, 3.0);
  printf(" grade larga (dt = 1.0, janela t = 500 — travessia do sistema):\n");
  run(1300,  1.0, 0.2);
  run(10000, 1.0, 3.0);
  printf(" acima do teto interativo:\n");
  run(50000, 1.0, 1e9);
  return 0;
}
