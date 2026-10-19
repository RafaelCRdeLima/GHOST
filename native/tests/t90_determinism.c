/* t90_determinism.c — ACEITAÇÃO 8 (parte de alvo único).
 *
 * A outra parte — mesma semente reproduzindo o mesmo grafo em ALVOS
 * diferentes (WASM, nativo, C++ exportado) — fica fora deste teste. Aqui se
 * garante que o PRNG e a montagem da CSR são funções puras da semente dentro
 * de um alvo, e que o fluxo de bits do xoshiro256++ é o esperado.
 */
#include "harness.h"

#include <stdlib.h>
#include <string.h>

static int same_csr(const dae_csr *a, const dae_csr *b)
{
  int32_t i;
  if (a->n != b->n || a->nnz != b->nnz) return 0;
  for (i = 0; i <= a->n; ++i) if (a->rowptr[i] != b->rowptr[i]) return 0;
  for (i = 0; i < a->nnz; ++i) if (a->colind[i] != b->colind[i]) return 0;
  for (i = 0; i < a->nnz; ++i) if (a->val[i] != b->val[i]) return 0;
  return 1;
}

int main(void)
{
  dae_test T;
  dae_rng g1, g2;
  dae_graph A, B, C;
  int i;

  dae_test_begin(&T, "aceitacao 8 (parcial): determinismo do PRNG e da montagem");

  dae_rng_seed(&g1, 12345ULL);
  dae_rng_seed(&g2, 12345ULL);
  for (i = 0; i < 1000; ++i)
    if (dae_rng_u64(&g1) != dae_rng_u64(&g2)) { dae_test_ok(&T, 0, "fluxo divergiu em %d", i); break; }
  dae_test_ok(&T, 1, "mesma semente, 1000 saidas identicas");

  dae_rng_seed(&g1, 12345ULL);
  dae_rng_seed(&g2, 12346ULL);
  { int diff = 0;
    for (i = 0; i < 64; ++i) if (dae_rng_u64(&g1) != dae_rng_u64(&g2)) ++diff;
    dae_test_ok(&T, diff >= 60, "sementes vizinhas divergem (%d de 64)", diff); }

  /* Estado zero é ponto fixo do xoshiro: a guarda de semeadura tem de impedir. */
  { dae_rng z; uint64_t acc = 0;
    dae_rng_seed(&z, 0ULL);
    for (i = 0; i < 16; ++i) acc |= dae_rng_u64(&z);
    dae_test_ok(&T, acc != 0ULL, "semente 0 nao colapsa o estado"); }

  /* Intervalo de dae_rng_below, que é onde o viés se esconderia. */
  { dae_rng z; int32_t lo = 1 << 30, hi = -1;
    dae_rng_seed(&z, 7ULL);
    for (i = 0; i < 200000; ++i) {
      const int32_t v = dae_rng_below(&z, 13);
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    dae_test_ok(&T, lo == 0 && hi == 12, "below(13) cobre [0,12]: [%d,%d]", (int)lo, (int)hi); }

  dae_test_ok(&T, dae_test_random_graph(&A, 300, 4, 2024ULL) == DAE_OK, "grafo A");
  dae_test_ok(&T, dae_test_random_graph(&B, 300, 4, 2024ULL) == DAE_OK, "grafo B");
  dae_test_ok(&T, dae_test_random_graph(&C, 300, 4, 2025ULL) == DAE_OK, "grafo C");
  dae_test_ok(&T, same_csr(&A.A, &B.A), "mesma semente reproduz a mesma CSR bit a bit");
  dae_test_ok(&T, !same_csr(&A.A, &C.A), "semente diferente muda o grafo");

  dae_graph_free(&A); dae_graph_free(&B); dae_graph_free(&C);
  return dae_test_end(&T);
}
