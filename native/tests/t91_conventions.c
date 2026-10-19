/* t91_conventions.c — invariantes que o CONVENTIONS.md promete.
 *
 * Cada bloco aqui corresponde a uma decisão de projeto que, se for quebrada,
 * quebra em silêncio: nenhum teste de física falha, e os números saem
 * plausíveis e errados.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

int main(void)
{
  dae_test T;
  dae_test_begin(&T, "convencoes: dedup, concurrence, NaN, cota espectral, cache de dt");

  /* --------------------------------------------------------------------
   * 1. Aresta repetida é DESCARTADA, nunca somada. Somar dobraria j_perp em
   *    alguns sítios sem nenhum aviso.
   * ------------------------------------------------------------------ */
  {
    const int32_t ei[5] = { 0, 1, 1, 0, 2 };
    const int32_t ej[5] = { 1, 2, 0, 1, 3 };   /* (0,1) tres vezes, contando (1,0) */
    const double  w[5]  = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    dae_csr A;
    int32_t dropped = -1, p;
    double v01 = 0.0;
    dae_test_ok(&T, dae_csr_from_edges(&A, 4, ei, ej, w, 5, &dropped) == DAE_OK, "montagem");
    dae_test_ok(&T, dropped == 4, "4 entradas dirigidas descartadas, obtido %d", (int)dropped);
    dae_test_ok(&T, A.nnz == 6, "nnz = 6 (3 arestas simples), obtido %d", (int)A.nnz);
    for (p = A.rowptr[0]; p < A.rowptr[1]; ++p) if (A.colind[p] == 1) v01 = A.val[p];
    dae_test_near(&T, v01, 1.0, 0.0, "peso NAO somado em (0,1): ");
    /* linhas ordenadas e sem repetição */
    { int32_t i, ok = 1;
      for (i = 0; i < A.n; ++i)
        for (p = A.rowptr[i] + 1; p < A.rowptr[i + 1]; ++p)
          if (A.colind[p] <= A.colind[p - 1]) ok = 0;
      dae_test_ok(&T, ok, "colind estritamente crescente em cada linha"); }
    dae_csr_free(&A);
  }

  /* --------------------------------------------------------------------
   * 2. Concurrence por módulo: o bloco diagonal não pode conter o termo i == j
   *    do produto externo. A identidade que fecha é
   *      soma do triangulo superior COM a diagonal  ==  C_l1.
   * ------------------------------------------------------------------ */
  {
    const int32_t n = 12, nmod = 3;
    int32_t mo[12], i, j;
    double re[12], im[12], acc = 0.0, nrm = 0.0;
    dae_obs_cfg cfg;
    dae_obs O;
    dae_rng g;
    dae_rng_seed(&g, 314ULL);
    for (i = 0; i < n; ++i) {
      mo[i] = i % nmod;
      re[i] = dae_rng_uniform(&g) - 0.5;
      im[i] = dae_rng_uniform(&g) - 0.5;
      nrm += re[i] * re[i] + im[i] * im[i];
    }
    nrm = sqrt(nrm);
    for (i = 0; i < n; ++i) { re[i] /= nrm; im[i] /= nrm; }

    cfg.module_of = mo; cfg.nmod = nmod; cfg.target = 5;
    cfg.want_pop = 1; cfg.want_conc_mod = 1; cfg.want_conc_full = 1;
    dae_test_ok(&T, dae_obs_alloc(&O, &cfg, n) == DAE_OK, "obs_alloc");
    dae_test_ok(&T, dae_obs_eval(re, im, n, &cfg, &O) == DAE_OK, "obs_eval");
    dae_test_near(&T, O.norm, 1.0, 1e-14, "norma: ");

    for (i = 0; i < nmod; ++i)
      for (j = i; j < nmod; ++j)
        acc += O.conc_mod[(size_t)i * (size_t)nmod + (size_t)j];
    dae_test_near(&T, acc, O.coh_l1, 1e-12,
                  "soma do triangulo superior da conc_mod = C_l1: ");

    /* e a matriz completa, somada sobre pares nao-ordenados, dá o mesmo */
    { double full = 0.0;
      for (i = 0; i < n; ++i)
        for (j = i + 1; j < n; ++j) full += O.conc_full[(size_t)i * (size_t)n + (size_t)j];
      dae_test_near(&T, full, O.coh_l1, 1e-12, "soma de C_ij sobre i<j = C_l1: "); }

    /* o bloco diagonal sem o termo espurio */
    { const double s0 = O.l1mod[0], q0 = O.pmod[0];
      dae_test_near(&T, O.conc_mod[0], s0 * s0 - q0, 1e-15,
                    "C_MM = s_M^2 - q_M (sem o termo i==j): "); }
    dae_obs_free(&O);
  }

  /* --------------------------------------------------------------------
   * 3. Sem alvo, p_alvo é NaN — nunca zero, que é valor fisicamente válido e
   *    acabaria plotado como dado.
   * ------------------------------------------------------------------ */
  {
    const int32_t n = 4;
    double re[4] = { 1.0, 0.0, 0.0, 0.0 }, im[4] = { 0.0, 0.0, 0.0, 0.0 };
    dae_obs_cfg cfg;
    dae_obs O;
    cfg.module_of = NULL; cfg.nmod = 1; cfg.target = -1;
    cfg.want_pop = 0; cfg.want_conc_mod = 0; cfg.want_conc_full = 0;
    dae_obs_alloc(&O, &cfg, n);
    dae_obs_eval(re, im, n, &cfg, &O);
    dae_test_ok(&T, O.p_target != O.p_target, "p_alvo = NaN sem alvo");
    cfg.target = 2;
    dae_obs_eval(re, im, n, &cfg, &O);
    dae_test_near(&T, O.p_target, 0.0, 0.0, "p_alvo = 0 no sitio 2, que esta vazio: ");
    dae_obs_free(&O);
  }

  /* --------------------------------------------------------------------
   * 4. A cota de Gershgorin é rigorosa, e a rede de segurança da norma
   *    realmente dispara quando `a` fica pequeno demais. Sem esta segunda
   *    parte, um limite espectral apertado devolveria lixo silencioso.
   * ------------------------------------------------------------------ */
  {
    dae_gen_params p;
    dae_graph G;
    dae_csr H;
    dae_cheb W;
    double re[64], im[64];
    int32_t i;
    dae_status st;

    dae_gen_params_default(&p);
    p.kind = DAE_G_CYCLE; p.n = 64;
    dae_graph_build(&G, &p);
    dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
    dae_cheb_init(&W, &H, 0);
    dae_test_near(&T, W.lo, -2.0, 1e-15, "Gershgorin inferior do ciclo: ");
    dae_test_near(&T, W.hi,  2.0, 1e-15, "Gershgorin superior do ciclo: ");

    for (i = 0; i < 64; ++i) { re[i] = 0.0; im[i] = 0.0; }
    re[0] = 1.0;
    W.a *= 0.5;                    /* sabota o limite: espectro sai de [-1,1] */
    W.coef_valid = 0;
    st = dae_cheb_step(&W, &H, 8.0, re, im, NULL);
    dae_test_ok(&T, st == DAE_ERR_NORM,
                "limite espectral sabotado aborta com DAE_ERR_NORM, veio '%s'",
                dae_strerror(st));
    dae_cheb_free(&W);

    /* laplaciana: linhas somam zero, por construção */
    { dae_csr L;
      double worst = 0.0;
      int32_t q;
      dae_hamiltonian(&L, &G.A, DAE_H_LAPLACIAN, 1.0, DAE_NORM_NONE, 0, NULL);
      for (i = 0; i < L.n; ++i) {
        double s = 0.0;
        for (q = L.rowptr[i]; q < L.rowptr[i + 1]; ++q) s += L.val[q];
        if (fabs(s) > worst) worst = fabs(s);
      }
      dae_test_ok(&T, worst < 1e-14, "linhas da laplaciana somam zero (%.3g)", worst);
      dae_csr_free(&L); }

    dae_csr_free(&H); dae_graph_free(&G);
  }

  /* --------------------------------------------------------------------
   * 5. O cache dos c_k é premissa de GRADE UNIFORME. Trocar dt entre passos
   *    tem de recalcular os coeficientes, não reutilizar os antigos.
   * ------------------------------------------------------------------ */
  {
    dae_gen_params p;
    dae_graph G;
    dae_csr H;
    dae_cheb W1, W2;
    double a_re[64], a_im[64], b_re[64], b_im[64], worst = 0.0;
    int32_t i;

    dae_gen_params_default(&p);
    p.kind = DAE_G_CYCLE; p.n = 64;
    dae_graph_build(&G, &p);
    dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);

    for (i = 0; i < 64; ++i) { a_re[i] = b_re[i] = 0.0; a_im[i] = b_im[i] = 0.0; }
    a_re[3] = 1.0; b_re[3] = 1.0;

    dae_cheb_init(&W1, &H, 0);
    dae_cheb_step(&W1, &H, 1.5, a_re, a_im, NULL);   /* dt alternado */
    dae_cheb_step(&W1, &H, 0.4, a_re, a_im, NULL);
    dae_cheb_step(&W1, &H, 1.5, a_re, a_im, NULL);

    dae_cheb_init(&W2, &H, 0);
    dae_cheb_step(&W2, &H, 3.4, b_re, b_im, NULL);   /* 1.5 + 0.4 + 1.5 */

    for (i = 0; i < 64; ++i) {
      if (fabs(a_re[i] - b_re[i]) > worst) worst = fabs(a_re[i] - b_re[i]);
      if (fabs(a_im[i] - b_im[i]) > worst) worst = fabs(a_im[i] - b_im[i]);
    }
    dae_test_ok(&T, worst < 1e-12, "dt alternado recalcula os c_k (desvio %.3g)", worst);

    dae_cheb_free(&W1); dae_cheb_free(&W2);
    dae_csr_free(&H); dae_graph_free(&G);
  }

  return dae_test_end(&T);
}
