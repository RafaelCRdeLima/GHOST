/* t94_metrics.c — métricas de rede, com atenção ao lambda_2.
 *
 * lambda_2 vira resultado de manchete (reportado contra p), e o regime que
 * interessa — rede fortemente modular — é o mais lento para convergir. Um
 * numero fixo de passos devolveria lambda_2 SUPERESTIMADO ali, achatando a
 * curva justamente na regiao interessante, e um plato numerico e
 * indistinguivel de um plato fisico. Estes testes cobrem: valor contra forma
 * fechada, valor contra diagonalizacao exata em grafo modular, o residuo como
 * cota rigorosa, e o comportamento da bandeira quando NAO converge.
 */
#include "harness.h"
#include "jacobi.h"

#include <math.h>
#include <stdlib.h>

/* lambda_2 exato por Jacobi denso: segundo menor autovalor de L. */
static double fiedler_exato(const dae_graph *G)
{
  const int32_t n = G->A.n;
  double *dense, *vec, *val, l1, l2;
  int32_t i, p;
  dae_csr L;

  dae_hamiltonian(&L, &G->A, DAE_H_LAPLACIAN, 1.0, DAE_NORM_NONE, 0, NULL);
  dense = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
  vec   = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
  val   = (double *)malloc((size_t)n * sizeof(double));
  for (i = 0; i < n; ++i)
    for (p = L.rowptr[i]; p < L.rowptr[i + 1]; ++p)
      dense[(size_t)i * (size_t)n + (size_t)L.colind[p]] = L.val[p];
  dae_jacobi(dense, vec, val, n);
  l1 = 1e300; l2 = 1e300;
  for (i = 0; i < n; ++i) {
    if (val[i] < l1) { l2 = l1; l1 = val[i]; }
    else if (val[i] < l2) l2 = val[i];
  }
  free(dense); free(vec); free(val);
  dae_csr_free(&L);
  return l2;
}

int main(void)
{
  dae_test T;
  dae_metrics_cfg cfg;
  dae_metrics M;
  dae_gen_params p;
  dae_graph G;
  const double TAU = 6.283185307179586;

  dae_test_begin(&T, "metricas: lambda_2 com convergencia, Q, caminho medio");
  dae_metrics_cfg_default(&cfg);

  /* --- formas fechadas --- */
  dae_gen_params_default(&p);
  p.kind = DAE_G_CYCLE; p.n = 64;
  dae_graph_build(&G, &p);
  dae_metrics_compute(&G, &cfg, &M);
  dae_test_ok(&T, M.lambda2_converged, "ciclo: convergiu em %d passos", (int)M.lambda2_steps);
  dae_test_near(&T, M.lambda2, 2.0 * (1.0 - cos(TAU / 64.0)), 1e-10, "ciclo C_64 lambda2: ");
  dae_test_near(&T, M.mean_path_len, (64.0 / 2.0) * (64.0 / 2.0) / 63.0, 1e-12,
                "ciclo C_64 caminho medio: ");
  dae_test_near(&T, M.mean_degree, 2.0, 1e-15, "ciclo: grau medio: ");
  dae_test_ok(&T, M.n_edges == 64, "ciclo: |E| = 64 (%d)", (int)M.n_edges);
  dae_test_note("ciclo C_64: lambda2=%.10g  residuo=%.2g  passos=%d",
                M.lambda2, M.lambda2_residual, (int)M.lambda2_steps);
  dae_graph_free(&G);

  p.kind = DAE_G_PATH; p.n = 64;
  dae_graph_build(&G, &p);
  dae_metrics_compute(&G, &cfg, &M);
  dae_test_near(&T, M.lambda2, 2.0 * (1.0 - cos(3.141592653589793 / 64.0)), 1e-10,
                "linha P_64 lambda2: ");
  dae_test_near(&T, M.mean_path_len, 65.0 / 3.0, 1e-12, "linha P_64 caminho medio: ");
  dae_graph_free(&G);

  p.kind = DAE_G_COMPLETE; p.n = 50;
  dae_graph_build(&G, &p);
  dae_metrics_compute(&G, &cfg, &M);
  dae_test_near(&T, M.lambda2, 50.0, 1e-9, "K_50 lambda2 = N: ");
  dae_test_near(&T, M.mean_path_len, 1.0, 1e-15, "K_50 caminho medio = 1: ");
  dae_graph_free(&G);

  p.kind = DAE_G_HYPERCUBE; p.dim = 6;
  dae_graph_build(&G, &p);
  dae_metrics_compute(&G, &cfg, &M);
  dae_test_near(&T, M.lambda2, 2.0, 1e-10, "hipercubo Q_6 lambda2 = 2: ");
  dae_graph_free(&G);

  /* --- regime modular: contra diagonalizacao exata --- */
  {
    const double pouts[3] = { 0.05, 0.01, 0.002 };
    int i;
    for (i = 0; i < 3; ++i) {
      double exato;
      dae_gen_params_default(&p);
      p.kind = DAE_G_SBM; p.n = 200; p.n_modules = 4;
      p.p_in = 0.35; p.p_out = pouts[i]; p.seed = 2026ULL;
      dae_graph_build(&G, &p);
      dae_metrics_compute(&G, &cfg, &M);
      exato = fiedler_exato(&G);
      if (M.n_components == 1) {
        dae_test_ok(&T, M.lambda2_converged,
                    "SBM p_out=%.3f: convergiu (%d passos, residuo %.2g)",
                    pouts[i], (int)M.lambda2_steps, M.lambda2_residual);
        dae_test_near(&T, M.lambda2, exato, 1e-8,
                      "SBM p_out=%.3f lambda2 contra Jacobi: ", pouts[i]);
        /* o residuo tem de ser cota RIGOROSA do erro */
        dae_test_ok(&T, fabs(M.lambda2 - exato) <= M.lambda2_residual + 1e-12,
                    "SBM p_out=%.3f: erro %.3g dentro do residuo %.3g",
                    pouts[i], fabs(M.lambda2 - exato), M.lambda2_residual);
      }
      dae_test_note("SBM p_out=%.3f: lambda2=%.8g exato=%.8g  Q=%.4f  passos=%d comp=%d",
                    pouts[i], M.lambda2, exato, M.modularity_Q,
                    (int)M.lambda2_steps, (int)M.n_components);
      dae_graph_free(&G);
    }
  }

  /* --- o caso que motiva a bandeira: poucos passos NAO convergem, e o valor
         devolvido e um LIMITE SUPERIOR, nunca uma medida --- */
  {
    dae_metrics_cfg curta;
    dae_metrics Mc;
    double exato;
    dae_gen_params_default(&p);
    p.kind = DAE_G_SBM; p.n = 200; p.n_modules = 4;
    p.p_in = 0.35; p.p_out = 0.002; p.seed = 2026ULL;
    dae_graph_build(&G, &p);
    exato = fiedler_exato(&G);

    dae_metrics_cfg_default(&curta);
    curta.lambda2_max_steps = 2;
    dae_metrics_compute(&G, &curta, &Mc);
    if (Mc.n_components == 1) {
      dae_test_ok(&T, !Mc.lambda2_converged,
                  "2 passos: bandeira diz que NAO convergiu");
      dae_test_ok(&T, Mc.lambda2 >= exato - 1e-12,
                  "2 passos: %.6g e limite SUPERIOR de %.6g (Rayleigh)",
                  Mc.lambda2, exato);
      dae_test_ok(&T, Mc.lambda2 > exato * 1.5,
                  "2 passos superestima mesmo (%.4g contra %.4g) — e por isso "
                  "que passos fixos achatariam a curva lambda2(p)",
                  Mc.lambda2, exato);
      dae_test_note("2 passos: lambda2=%.6g (residuo %.2g) contra exato %.6g",
                    Mc.lambda2, Mc.lambda2_residual, exato);
    }
    dae_graph_free(&G);
  }

  /* --- Q do microtubulo particionado --- */
  {
    dae_gen_params_default(&p);
    p.kind = DAE_G_MICROTUBULE; p.n_par = 80; p.n_perp = 13; p.n_modules = 4;
    dae_graph_build(&G, &p);
    dae_metrics_compute(&G, &cfg, &M);
    dae_test_ok(&T, M.modularity_Q > 0.6 && M.modularity_Q < 0.76,
                "microtubulo em 4 modulos: Q = %.4f", M.modularity_Q);
    dae_test_ok(&T, M.n_components == 1, "microtubulo conexo");
    dae_test_note("microtubulo 80x13, M=4: Q=%.4f  lambda2=%.6g  grau medio=%.3f  L=%.2f",
                  M.modularity_Q, M.lambda2, M.mean_degree, M.mean_path_len);
    dae_graph_free(&G);
  }

  /* --- ATAQUE AO MECANISMO, nao ao sintoma ---
   *
   * lambda_2 errado e sintoma tardio e mascarado: os quatro pontos de deflacao
   * de dae_fiedler sao redundantes entre si quanto a ELE, entao remover um so
   * nao quebra nada. A quantidade que todos existem para segurar e a projecao
   * dos vetores de Lanczos sobre o vetor constante. Olhando ela, a erosao
   * aparece cedo, e um refactor que apague qualquer um dos pontos morde aqui.
   */
  {
    struct { const char *nome; dae_gen_kind k; int32_t n, np, nq, dim, M;
             double pin, pout, ws; } casos[] = {
      { "linha P_64",        DAE_G_PATH,        64, 0,  0,  0, 1, 0, 0, 0.0 },
      { "ciclo C_64",        DAE_G_CYCLE,       64, 0,  0,  0, 1, 0, 0, 0.0 },
      { "K_50",              DAE_G_COMPLETE,    50, 0,  0,  0, 1, 0, 0, 0.0 },
      { "hipercubo Q_7",     DAE_G_HYPERCUBE,    0, 0,  0,  7, 1, 0, 0, 0.0 },
      { "SBM modular",       DAE_G_SBM,        200, 0,  0,  0, 4, 0.35, 0.002, 0.0 },
      { "SBM quase disjunto",DAE_G_SBM,        200, 0,  0,  0, 4, 0.40, 0.006, 0.0 },
      { "microtubulo 40x13", DAE_G_MICROTUBULE,  0, 40, 13, 0, 4, 0, 0, 0.0 },
      { "microtubulo religado", DAE_G_MICROTUBULE, 0, 40, 13, 0, 4, 0, 0, 0.1 }
    };
    const int ncasos = (int)(sizeof(casos) / sizeof(casos[0]));
    int c;
    for (c = 0; c < ncasos; ++c) {
      dae_metrics Mm;
      dae_gen_params_default(&p);
      p.kind = casos[c].k; p.n = casos[c].n; p.n_par = casos[c].np;
      p.n_perp = casos[c].nq; p.dim = casos[c].dim; p.n_modules = casos[c].M;
      p.p_in = casos[c].pin; p.p_out = casos[c].pout; p.ws_p = casos[c].ws;
      p.seed = 314159ULL;
      if (dae_graph_build(&G, &p) != DAE_OK) continue;
      dae_metrics_compute(&G, &cfg, &Mm);
      /* Envoltoria, nao constante: o vazamento acumula com os passos, entao o
         limite tem de ser uma lei. Medido no codigo sadio, o pior caso e
         ~1.2e-17 por passo; a envoltoria da ~4x de folga e ainda separa por
         uma ordem de grandeza do codigo com qualquer deflacao removida. */
      { const double env = 5e-17 * (double)Mm.lambda2_steps + 2e-16;
        dae_test_ok(&T, Mm.lambda2_const_leak < env,
                    "%s: vazamento da direcao constante = %.3g acima da "
                    "envoltoria %.3g em %d passos — alguma das deflacoes de "
                    "dae_fiedler sumiu",
                    casos[c].nome, Mm.lambda2_const_leak, env,
                    (int)Mm.lambda2_steps);
        dae_test_note("%-24s vazamento=%8.2g  envoltoria=%8.2g  passos=%3d  lambda2=%.6g",
                      casos[c].nome, Mm.lambda2_const_leak, env,
                      (int)Mm.lambda2_steps, Mm.lambda2); }
      dae_graph_free(&G);
    }
  }

  /* --- REGIME DE EXAUSTAO ---
   *
   * Os casos acima convergem antes de esgotar o Krylov, e por isso nao exercem
   * as defesas: o vazamento so tem chance de aparecer quando a iteracao chega
   * ao fim do espaco. Com tolerancia inatingivel, P_64 roda ate o teto — que
   * TEM de ser n-1, porque o espaco deflacionado tem dimensao n-1 e a n-esima
   * direcao so pode ser ruido contendo o vetor constante.
   */
  {
    dae_metrics_cfg exausta;
    dae_metrics Me;
    double exato;
    const int32_t n = 64;
    dae_metrics_cfg_default(&exausta);
    exausta.lambda2_tol = 1e-16;           /* inatingivel de proposito */
    dae_gen_params_default(&p);
    p.kind = DAE_G_PATH; p.n = n;
    dae_graph_build(&G, &p);
    exato = 2.0 * (1.0 - cos(3.141592653589793 / (double)n));
    dae_metrics_compute(&G, &exausta, &Me);

    dae_test_ok(&T, Me.lambda2_steps <= n - 1,
                "exaustao: %d passos, teto tem de ser n-1 = %d",
                (int)Me.lambda2_steps, (int)(n - 1));
    dae_test_ok(&T, Me.lambda2_const_leak < 5e-17 * (double)Me.lambda2_steps + 2e-16,
                "exaustao: vazamento = %.3g apos %d passos",
                Me.lambda2_const_leak, (int)Me.lambda2_steps);
    dae_test_ok(&T, Me.lambda2 >= exato - 1e-12,
                "exaustao: %.10g continua sendo limite superior de %.10g",
                Me.lambda2, exato);
    dae_test_near(&T, Me.lambda2, exato, 1e-9,
                  "exaustao: mesmo sem bater a tolerancia, o valor esta certo: ");
    dae_test_note("exaustao P_64: passos=%d  vazamento=%.2g  lambda2=%.10g  residuo=%.2g",
                  (int)Me.lambda2_steps, Me.lambda2_const_leak, Me.lambda2,
                  Me.lambda2_residual);
    dae_graph_free(&G);
  }

  return dae_test_end(&T);
}
