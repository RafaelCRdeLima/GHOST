/* dump.c — roda um cenario e imprime tudo que o teste 6 compara.
 *
 * Compilado DUAS vezes a partir do mesmo fonte: nativo pelo gcc e WebAssembly
 * pelo emcc. O teste 6 confronta as duas saidas. Nao ha "versao do navegador"
 * e "versao do cluster": ha um fonte e dois compiladores.
 */
#include "scenarios.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void imprime(const dae_scenario *S, int idx)
{
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  dae_cheb_info info;
  dae_obs_cfg cfg;
  dae_obs O;
  dae_metrics M;
  dae_metrics_cfg mcfg;
  double *re, *im, escala = 1.0, dt;
  int32_t i, s;
  dae_status st;

  st = dae_graph_build(&G, &S->gen);
  if (st != DAE_OK) { printf("ERRO grafo %s\n", dae_strerror(st)); return; }

  printf("# cenario %d %s\n", idx, S->nome);
  printf("grafo n %d nnz %d nmod %d descartadas %d religacoes_falhas %d\n",
         (int)G.n, (int)G.A.nnz, (int)G.nmod, (int)G.n_dropped, (int)G.n_rewire_failed);

  printf("grafo digital %llu\n", (unsigned long long)dae_graph_fingerprint(&G));

  dae_metrics_cfg_default(&mcfg);
  dae_metrics_compute(&G, &mcfg, &M);
  printf("metrica lambda2 %.17g residuo %.17g convergiu %d passos %d vazamento %.17g\n",
         M.lambda2, M.lambda2_residual, M.lambda2_converged,
         (int)M.lambda2_steps, M.lambda2_const_leak);
  printf("metrica Q %.17g grau %.17g caminho %.17g arestas %d componentes %d\n",
         M.modularity_Q, M.mean_degree, M.mean_path_len,
         (int)M.n_edges, (int)M.n_components);

  st = dae_hamiltonian(&H, &G.A, S->ham, S->gamma, S->norm, S->lanczos, &escala);
  if (st != DAE_OK) { printf("ERRO H %s\n", dae_strerror(st)); dae_graph_free(&G); return; }
  printf("hamiltoniano escala %.17g\n", escala);

  st = dae_cheb_init(&W, &H, S->lanczos);
  if (st != DAE_OK) { printf("ERRO cheb %s\n", dae_strerror(st)); return; }
  printf("cheb lo %.17g hi %.17g a %.17g b %.17g lanczos %d\n",
         W.lo, W.hi, W.a, W.b, W.lanczos_used);

  cfg.module_of = G.module_of; cfg.nmod = G.nmod; cfg.target = S->target;
  cfg.want_pop = 0; cfg.want_conc_mod = 1; cfg.want_conc_full = 0;
  dae_obs_alloc(&O, &cfg, G.n);

  re = (double *)calloc((size_t)G.n, sizeof(double));
  im = (double *)calloc((size_t)G.n, sizeof(double));
  if (!re || !im) return;
  re[S->init_site] = 1.0;

  dt = S->t1 / (double)S->nt;
  info.k_used = 0;
  for (s = 0; s < S->nt; ++s) {
    st = dae_cheb_step(&W, &H, dt, re, im, &info);
    if (st != DAE_OK) { printf("ERRO passo %d %s\n", (int)s, dae_strerror(st)); break; }
    dae_obs_eval(re, im, G.n, &cfg, &O);
    printf("passo %d K %d norma %.17g ipr %.17g coh %.17g palvo %.17g\n",
           (int)s, (int)info.k_used, O.norm, O.ipr, O.coh_l1, O.p_target);
    for (i = 0; i < O.nmod; ++i)
      printf("pmod %d %d %.17g %.17g\n", (int)s, (int)i, O.pmod[i], O.l1mod[i]);
  }
  for (i = 0; i < G.n; ++i) printf("psi %d %.17g %.17g\n", (int)i, re[i], im[i]);

  free(re); free(im);
  dae_obs_free(&O); dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
}

int main(int argc, char **argv)
{
  int i;
  printf("# ghost %s nucleo %s\n", DAE_VERSION, DAE_CORE_HASH);
  if (argc > 1) {
    const int idx = atoi(argv[1]);
    const dae_scenario *S = dae_scenario_get(idx);
    if (!S) { printf("ERRO cenario inexistente\n"); return 1; }
    imprime(S, idx);
    return 0;
  }
  for (i = 0; i < dae_scenario_count(); ++i) imprime(dae_scenario_get(i), i);
  return 0;
}
