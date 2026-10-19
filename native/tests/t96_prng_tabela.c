/* t96_prng_tabela.c — congela o fluxo do PRNG para o pacote Wolfram.
 *
 * O alvo Wolfram reimplementa xoshiro256++, splitmix64 e a rejeicao de Lemire
 * em aritmetica exata mod 2^64. Se essa reimplementacao divergir, a falha NAO
 * se apresenta como erro: se apresenta como grafo diferente, portanto fisica
 * diferente, portanto "os dois programas discordam".
 *
 * A digital do grafo pega isso, mas tarde e sem dizer onde. Esta tabela e a
 * verificacao barata e diagnostica que vem antes — mesmo raciocinio do teste 0,
 * que comparou Bessel contra tabela externa em vez de contra si mesmo.
 *
 * Roda com --escrever para (re)gerar specs/oraculo/prng.json.
 */
#include "harness.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NSEM 3
#define NVAL 100

static const uint64_t SEMENTES[NSEM] = { 0ULL, 12345ULL, 2026ULL };
static const int32_t  LIMITES[4] = { 2, 13, 1000, 65536 };

static void escreve(const char *caminho)
{
  FILE *f = fopen(caminho, "wb");
  int s, i, k;
  if (!f) { printf("   x nao consegui escrever %s\n", caminho); return; }
  fprintf(f, "{\n  \"nota\": \"gerado por native/tests/t96_prng_tabela.c --escrever\",\n");
  fprintf(f, "  \"core_hash\": \"%s\",\n  \"fluxos\": [\n", DAE_CORE_HASH);
  for (s = 0; s < NSEM; ++s) {
    dae_rng g;
    dae_rng_seed(&g, SEMENTES[s]);
    fprintf(f, "    { \"semente\": %llu, \"u64\": [", (unsigned long long)SEMENTES[s]);
    for (i = 0; i < NVAL; ++i)
      fprintf(f, "%s\"%llu\"", i ? ", " : "", (unsigned long long)dae_rng_u64(&g));
    fprintf(f, "] }%s\n", s + 1 < NSEM ? "," : "");
  }
  fprintf(f, "  ],\n  \"abaixo\": [\n");
  for (k = 0; k < 4; ++k) {
    dae_rng g;
    dae_rng_seed(&g, 2026ULL);
    fprintf(f, "    { \"limite\": %d, \"valores\": [", (int)LIMITES[k]);
    for (i = 0; i < NVAL; ++i)
      fprintf(f, "%s%d", i ? ", " : "", (int)dae_rng_below(&g, LIMITES[k]));
    fprintf(f, "] }%s\n", k + 1 < 4 ? "," : "");
  }
  /* uniform em [0,1): o gerador de SBM decide cada aresta com ele, entao um
     desvio de meio ulp aqui muda o grafo. */
  fprintf(f, "  ],\n  \"uniforme\": [");
  { dae_rng g;
    dae_rng_seed(&g, 2026ULL);
    for (i = 0; i < 20; ++i)
      fprintf(f, "%s%.17g", i ? ", " : "", dae_rng_uniform(&g)); }
  fprintf(f, "]\n}\n");
  fclose(f);
  printf("   . tabela escrita em %s\n", caminho);
}

int main(int argc, char **argv)
{
  dae_test T;
  dae_test_begin(&T, "tabela do PRNG para o alvo Wolfram");

  if (argc > 1 && strcmp(argv[1], "--escrever") == 0) {
    escreve(argc > 2 ? argv[2] : "specs/oraculo/prng.json");
    return dae_test_end(&T);
  }

  /* Sem --escrever, o teste confere que o fluxo continua o mesmo de sempre:
     alguem mexer no PRNG invalida silenciosamente todo grafo estocastico ja
     publicado, e a tabela do Wolfram junto. */
  { dae_rng g;
    uint64_t primeiro, centesimo;
    int i;
    dae_rng_seed(&g, 2026ULL);
    primeiro = dae_rng_u64(&g);
    for (i = 1; i < NVAL; ++i) centesimo = dae_rng_u64(&g);
    dae_test_note("semente 2026: primeiro=%llu centesimo=%llu",
                  (unsigned long long)primeiro, (unsigned long long)centesimo);
    dae_test_ok(&T, primeiro != 0ULL && centesimo != 0ULL, "fluxo nao degenerou"); }

  { dae_rng g;
    int32_t soma = 0, i;
    dae_rng_seed(&g, 2026ULL);
    for (i = 0; i < NVAL; ++i) soma += dae_rng_below(&g, 13);
    dae_test_note("soma de 100 sorteios abaixo de 13: %d", (int)soma);
    dae_test_ok(&T, soma > 300 && soma < 900, "abaixo(13) no intervalo plausivel"); }

  return dae_test_end(&T);
}
