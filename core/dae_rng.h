/* dae_rng.h — xoshiro256++, semeado explicitamente.
 *
 * CONTRATO DE REPRODUTIBILIDADE. Todo consumo de aleatoriedade do Daedalus
 * passa por aqui. Nunca Math.random() do JS: o grafo visto na tela tem de ser
 * o grafo que o cluster rodou.
 *
 * A ORDEM das chamadas faz parte do contrato, não só a semente. Todo sorteio
 * varre índices em ordem crescente; nenhum gerador pode sortear iterando sobre
 * uma tabela de dispersão ou um Set, porque a ordem mudaria entre alvos e a
 * semente deixaria de reproduzir o grafo.
 */
#ifndef DAE_RNG_H
#define DAE_RNG_H

#include "dae_types.h"

typedef struct { uint64_t s[4]; } dae_rng;

/* splitmix64 sobre a semente preenche o estado; nenhuma semente leva a zero. */
void     dae_rng_seed(dae_rng *g, uint64_t seed);
uint64_t dae_rng_u64(dae_rng *g);

/* [0,1) com 53 bits de mantissa. */
double   dae_rng_uniform(dae_rng *g);

/* Inteiro em [0,n) sem viés, por rejeição. Determinístico em qualquer alvo:
 * não usa __int128 nem ponto flutuante. */
int32_t  dae_rng_below(dae_rng *g, int32_t n);

#endif /* DAE_RNG_H */
