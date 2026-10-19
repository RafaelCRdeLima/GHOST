/* harness.h — mínimo para os testes de aceitação. Fora do núcleo: pode usar
 * stdio, e usa. */
#ifndef DAE_TEST_HARNESS_H
#define DAE_TEST_HARNESS_H

#include "dae.h"

typedef struct { const char *name; int checks; int fails; } dae_test;

void dae_test_begin(dae_test *T, const char *name);
int  dae_test_end(dae_test *T);                 /* código de saída do processo */
void dae_test_ok(dae_test *T, int cond, const char *fmt, ...);
/* |got - want| <= tol, com o valor impresso quando falha */
void dae_test_near(dae_test *T, double got, double want, double tol,
                   const char *fmt, ...);
void dae_test_note(const char *fmt, ...);

/* Grafo aleatório G(n, grau médio) com caminho geradora embutida, para os
 * testes 1 e 5. Usa o PRNG do núcleo, então é reprodutível. */
dae_status dae_test_random_graph(dae_graph *G, int32_t n, int32_t mean_deg,
                                 uint64_t seed);

double dae_test_seconds(void);

#endif /* DAE_TEST_HARNESS_H */
