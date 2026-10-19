/* t92_microtubule.c — lattice microtubular: costura e pontas.
 *
 * Conta arestas contra a fórmula fechada e verifica, sítio a sítio, ONDE a
 * costura helicoidal deixa buracos. Ver CONVENTIONS.md, parte 5.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

static int32_t grau(const dae_csr *A, int32_t j)
{
  return A->rowptr[j + 1] - A->rowptr[j];
}

static void caso(dae_test *T, int32_t np, int32_t nq, int32_t shift, int fechado)
{
  dae_gen_params p;
  dae_graph G;
  int32_t esperado, m, q, def_baixo = 0, def_alto = 0, outros = 0;
  dae_status st;

  dae_gen_params_default(&p);
  p.kind = DAE_G_MICROTUBULE;
  p.n_par = np; p.n_perp = nq;
  p.seam_shift = shift;
  p.longitudinal_closed = fechado;
  p.n_modules = 4;
  st = dae_graph_build(&G, &p);
  dae_test_ok(T, st == DAE_OK, "montagem: %s", dae_strerror(st));
  if (st != DAE_OK) return;

  dae_test_ok(T, G.n == np * nq, "N = N_par * N_perp = %d", (int)G.n);
  dae_test_ok(T, G.n_dropped == 0, "sem duplicatas descartadas (%d)", (int)G.n_dropped);

  /* longitudinais + transversais internas + costura */
  esperado = (fechado ? np : np - 1) * nq          /* longitudinais */
           + np * (nq - 1)                          /* transversais internas */
           + (fechado ? np : np - shift);           /* costura */
  dae_test_ok(T, G.A.nnz == 2 * esperado,
              "shift=%d fechado=%d: |E| esperado %d, nnz %d",
              (int)shift, fechado, (int)esperado, (int)G.A.nnz);

  /* ONDE faltam as ligações de costura. Uma aresta (m, nq-1)-(m+s, 0) some
     quando m+s >= N_par: some o extremo em COLUNA nq-1 do lado de m alto E o
     parceiro em COLUNA 0 do lado de m baixo. As duas pontas ficam deficientes,
     mas em protofilamentos DIFERENTES — é essa a inequivalência. */
  for (m = 0; m < np; ++m)
    for (q = 0; q < nq; ++q) {
      const int32_t j = m * nq + q;
      const int32_t d = grau(&G.A, j);
      const int32_t esperado_int = (fechado ? 4 : ((m == 0 || m == np - 1) ? 3 : 4));
      if (d == esperado_int) continue;
      if (q == 0)            ++def_baixo;
      else if (q == nq - 1)  ++def_alto;
      else                   ++outros;
    }
  dae_test_ok(T, outros == 0, "deficiencia so nas colunas 0 e nq-1 (%d fora)", (int)outros);
  if (!fechado) {
    dae_test_ok(T, def_baixo == shift,
                "shift=%d: %d sitios deficientes na coluna 0 (esperado %d)",
                (int)shift, (int)def_baixo, (int)shift);
    dae_test_ok(T, def_alto == shift,
                "shift=%d: %d sitios deficientes na coluna nq-1 (esperado %d)",
                (int)shift, (int)def_alto, (int)shift);
  }
  dae_test_note("shift=%d fechado=%d: |E|=%d  deficientes col0=%d col%d=%d",
                (int)shift, fechado, (int)esperado, (int)def_baixo,
                (int)(nq - 1), (int)def_alto);

  /* módulos: blocos contíguos ao longo do eixo longitudinal */
  { int ok = 1;
    for (m = 0; m < np; ++m)
      for (q = 0; q < nq; ++q)
        if (G.module_of[m * nq + q] != G.module_of[m * nq])
          ok = 0;                                  /* mesmo m => mesmo modulo */
    dae_test_ok(T, ok, "modulo constante ao longo de cada anel");
    dae_test_ok(T, G.module_of[0] == 0 && G.module_of[G.n - 1] == p.n_modules - 1,
                "modulos cobrem [0, M-1]"); }

  dae_graph_free(&G);
}

int main(void)
{
  dae_test T;
  dae_test_begin(&T, "microtubulo: costura, pontas abertas, modulos");

  caso(&T, 20, 13, 0, 0);       /* periodico ideal, pontas abertas */
  caso(&T, 20, 13, 3, 0);       /* helicoidal realista             */
  caso(&T, 20, 13, 3, 1);       /* helicoidal, cilindro fechado    */
  caso(&T, 50, 13, 5, 0);

  /* A costura muda o espectro: é o ponto de poder comparar os dois. */
  {
    dae_gen_params p;
    dae_graph G0, G3;
    dae_metrics m0, m3;
    dae_metrics_cfg cfg;
    dae_metrics_cfg_default(&cfg);
    dae_gen_params_default(&p);
    p.kind = DAE_G_MICROTUBULE; p.n_par = 40; p.n_perp = 13; p.n_modules = 4;
    p.seam_shift = 0; dae_graph_build(&G0, &p);
    p.seam_shift = 3; dae_graph_build(&G3, &p);
    dae_metrics_compute(&G0, &cfg, &m0);
    dae_metrics_compute(&G3, &cfg, &m3);
    dae_test_ok(&T, m0.lambda2_converged && m3.lambda2_converged, "lambda2 convergiu");
    dae_test_ok(&T, fabs(m0.lambda2 - m3.lambda2) > 1e-9,
                "seam_shift muda o espectro: lambda2 = %.9g contra %.9g",
                m0.lambda2, m3.lambda2);
    dae_test_note("seam=0: lambda2=%.6g  |E|=%d   seam=3: lambda2=%.6g  |E|=%d",
                  m0.lambda2, (int)m0.n_edges, m3.lambda2, (int)m3.n_edges);
    dae_graph_free(&G0); dae_graph_free(&G3);
  }

  /* determinismo: mesma semente, mesma CSR */
  {
    dae_gen_params p;
    dae_graph A, B;
    int igual = 1;
    int32_t i;
    dae_gen_params_default(&p);
    p.kind = DAE_G_MICROTUBULE; p.n_par = 30; p.n_perp = 13; p.ws_p = 0.2; p.seed = 99ULL;
    dae_graph_build(&A, &p);
    dae_graph_build(&B, &p);
    if (A.A.nnz != B.A.nnz) igual = 0;
    else for (i = 0; i < A.A.nnz; ++i)
           if (A.A.colind[i] != B.A.colind[i] || A.A.val[i] != B.A.val[i]) igual = 0;
    dae_test_ok(&T, igual, "mesma semente reproduz o microtubulo religado");
    dae_graph_free(&A); dae_graph_free(&B);
  }

  return dae_test_end(&T);
}
