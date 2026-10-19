/* t93_connectivity.c — religar contra acrescentar, e SBM.
 *
 * A disciplina de CONVENTIONS.md 3.1 diz que o modo padrão mantém |E| FIXO.
 * Se uma religação caísse sobre aresta existente e fosse aceita, a
 * deduplicação da CSR a descartaria e |E| cairia em silêncio — e transporte
 * pior por menos arestas seria lido como efeito de topologia.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

int main(void)
{
  dae_test T;
  const double ps[4] = { 0.0, 0.1, 0.5, 1.0 };
  int ip;

  dae_test_begin(&T, "conectividade: religar mantem |E|, acrescentar cresce");

  for (ip = 0; ip < 4; ++ip) {
    dae_gen_params p;
    dae_graph base, rew, add;
    dae_gen_params_default(&p);
    p.kind = DAE_G_MICROTUBULE; p.n_par = 60; p.n_perp = 13; p.seed = 7ULL;

    p.ws_p = 0.0;              dae_graph_build(&base, &p);
    p.ws_p = ps[ip]; p.conn_mode = DAE_REWIRE; dae_graph_build(&rew, &p);
    p.ws_p = ps[ip]; p.conn_mode = DAE_ADD;    dae_graph_build(&add, &p);

    dae_test_ok(&T, rew.A.nnz == base.A.nnz,
                "p=%.1f RELIGAR: nnz %d contra %d da base",
                ps[ip], (int)rew.A.nnz, (int)base.A.nnz);
    dae_test_ok(&T, rew.n_dropped == 0,
                "p=%.1f RELIGAR: nenhuma aresta perdida por duplicata (%d)",
                ps[ip], (int)rew.n_dropped);
    if (ps[ip] > 0.0)
      dae_test_ok(&T, add.A.nnz > base.A.nnz,
                  "p=%.1f ACRESCENTAR: nnz %d > %d",
                  ps[ip], (int)add.A.nnz, (int)base.A.nnz);
    dae_test_note("p=%.1f  base |E|=%d  religado |E|=%d  acrescentado |E|=%d  falhas=%d",
                  ps[ip], (int)(base.A.nnz / 2), (int)(rew.A.nnz / 2),
                  (int)(add.A.nnz / 2), (int)rew.n_rewire_failed);

    if (ps[ip] > 0.0) {
      int diff = 0;
      int32_t i;
      for (i = 0; i < base.A.nnz && i < rew.A.nnz; ++i)
        if (base.A.colind[i] != rew.A.colind[i]) diff = 1;
      dae_test_ok(&T, diff, "p=%.1f: a religacao mudou de fato a topologia", ps[ip]);
    }
    dae_graph_free(&base); dae_graph_free(&rew); dae_graph_free(&add);
  }

  /* SBM: casos limite exatos */
  {
    dae_gen_params p;
    dae_graph G;
    dae_gen_params_default(&p);
    p.kind = DAE_G_SBM; p.n = 120; p.n_modules = 3; p.seed = 11ULL;

    p.p_in = 1.0; p.p_out = 0.0;
    dae_graph_build(&G, &p);
    dae_test_ok(&T, G.A.nnz == 3 * 40 * 39,
                "p_in=1, p_out=0: tres K_40 disjuntos, nnz=%d", (int)G.A.nnz);
    { dae_metrics m; dae_metrics_cfg c; dae_metrics_cfg_default(&c);
      dae_metrics_compute(&G, &c, &m);
      dae_test_ok(&T, m.n_components == 3, "tres componentes (%d)", (int)m.n_components);
      dae_test_near(&T, m.lambda2, 0.0, 0.0, "desconexo => lambda2 = 0 exato: ");
      dae_test_ok(&T, m.lambda2_converged, "desconexo: convergencia trivial, sem iterar");
      dae_test_near(&T, m.modularity_Q, 2.0 / 3.0, 1e-12,
                    "Q de tres blocos completos disjuntos = 1 - 3*(1/3)^2: "); }
    dae_graph_free(&G);

    p.p_in = 0.0; p.p_out = 0.0;
    dae_graph_build(&G, &p);
    dae_test_ok(&T, G.A.nnz == 0, "p_in=p_out=0: grafo vazio (%d)", (int)G.A.nnz);
    dae_graph_free(&G);
  }

  /* determinismo do SBM */
  {
    dae_gen_params p;
    dae_graph A, B, C;
    int32_t i;
    int ab = 1, ac = 1;
    dae_gen_params_default(&p);
    p.kind = DAE_G_SBM; p.n = 200; p.n_modules = 4; p.p_in = 0.3; p.p_out = 0.02;
    p.seed = 5ULL; dae_graph_build(&A, &p);
    p.seed = 5ULL; dae_graph_build(&B, &p);
    p.seed = 6ULL; dae_graph_build(&C, &p);
    if (A.A.nnz != B.A.nnz) ab = 0;
    else for (i = 0; i < A.A.nnz; ++i) if (A.A.colind[i] != B.A.colind[i]) ab = 0;
    if (A.A.nnz != C.A.nnz) ac = 0;
    dae_test_ok(&T, ab, "SBM: mesma semente, mesmo grafo");
    dae_test_ok(&T, !ac, "SBM: semente diferente, grafo diferente");
    dae_graph_free(&A); dae_graph_free(&B); dae_graph_free(&C);
  }

  return dae_test_end(&T);
}
