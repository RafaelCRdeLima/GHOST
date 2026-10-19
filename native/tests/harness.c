#include "harness.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void dae_test_begin(dae_test *T, const char *name)
{
  T->name = name; T->checks = 0; T->fails = 0;
  printf("== %s  (ghost %s, core %s)\n", name, DAE_VERSION, DAE_CORE_HASH);
}

int dae_test_end(dae_test *T)
{
  if (T->fails == 0) { printf("   PASSOU  %d verificacoes\n", T->checks); return 0; }
  printf("   FALHOU  %d de %d verificacoes\n", T->fails, T->checks);
  return 1;
}

/* GHOST_VERBOSO=1 imprime também as verificações que PASSAM, com os valores
   que elas carregam na mensagem. */
static int verboso(void)
{
  const char *v = getenv("GHOST_VERBOSO");
  return v && v[0] && v[0] != '0';
}

void dae_test_ok(dae_test *T, int cond, const char *fmt, ...)
{
  ++T->checks;
  if (cond) {
    if (verboso()) { va_list ap; va_start(ap, fmt); printf("   v "); vprintf(fmt, ap); printf("\n"); va_end(ap); }
    return;
  }
  ++T->fails;
  { va_list ap; va_start(ap, fmt); printf("   x "); vprintf(fmt, ap); printf("\n"); va_end(ap); }
}

void dae_test_near(dae_test *T, double got, double want, double tol,
                   const char *fmt, ...)
{
  const double d = fabs(got - want);
  ++T->checks;
  if (d <= tol && (d == d)) return;
  ++T->fails;
  { va_list ap; va_start(ap, fmt); printf("   x "); vprintf(fmt, ap); va_end(ap); }
  printf("  obtido=%.17g esperado=%.17g desvio=%.3g tol=%.3g\n", got, want, d, tol);
}

void dae_test_note(const char *fmt, ...)
{
  va_list ap; va_start(ap, fmt); printf("   . "); vprintf(fmt, ap); printf("\n"); va_end(ap);
}

dae_status dae_test_random_graph(dae_graph *G, int32_t n, int32_t mean_deg,
                                 uint64_t seed)
{
  const int32_t extra = n * mean_deg / 2;
  const int32_t cap = (n - 1) + (extra > 0 ? extra : 0);
  int32_t *ei = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  int32_t *ej = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  dae_rng g;
  int32_t e = 0, i;
  dae_status st;

  if (!ei || !ej) { free(ei); free(ej); return DAE_ERR_ALLOC; }
  dae_rng_seed(&g, seed);
  for (i = 0; i + 1 < n; ++i) { ei[e] = i; ej[e] = i + 1; ++e; }   /* conexo */
  while (e < cap) {
    int32_t a = dae_rng_below(&g, n);
    int32_t b = dae_rng_below(&g, n);
    if (a == b) continue;
    ei[e] = a; ej[e] = b; ++e;
  }
  st = dae_graph_from_edges(G, n, ei, ej, NULL, cap);
  free(ei); free(ej);
  return st;
}

double dae_test_seconds(void)
{
  return (double)clock() / (double)CLOCKS_PER_SEC;
}
