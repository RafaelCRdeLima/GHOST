#include "dae_rng.h"

static uint64_t dae_rotl(uint64_t x, int k)
{
  return (x << k) | (x >> (64 - k));
}

void dae_rng_seed(dae_rng *g, uint64_t seed)
{
  int i;
  uint64_t z = seed;
  for (i = 0; i < 4; ++i) {
    /* splitmix64 */
    z += 0x9E3779B97F4A7C15ULL;
    {
      uint64_t t = z;
      t = (t ^ (t >> 30)) * 0xBF58476D1CE4E5B9ULL;
      t = (t ^ (t >> 27)) * 0x94D049BB133111EBULL;
      g->s[i] = t ^ (t >> 31);
    }
  }
  /* estado todo-zero é ponto fixo do xoshiro; splitmix não produz isso para
     nenhuma semente, mas a guarda é barata e o custo de errar é silencioso. */
  if ((g->s[0] | g->s[1] | g->s[2] | g->s[3]) == 0ULL) g->s[0] = 0x9E3779B97F4A7C15ULL;
}

uint64_t dae_rng_u64(dae_rng *g)
{
  const uint64_t result = dae_rotl(g->s[0] + g->s[3], 23) + g->s[0];
  const uint64_t t = g->s[1] << 17;
  g->s[2] ^= g->s[0];
  g->s[3] ^= g->s[1];
  g->s[1] ^= g->s[2];
  g->s[0] ^= g->s[3];
  g->s[2] ^= t;
  g->s[3] = dae_rotl(g->s[3], 45);
  return result;
}

double dae_rng_uniform(dae_rng *g)
{
  /* 2^-53 escrito como divisão exata para não depender de literal hexadecimal
     de ponto flutuante. */
  return (double)(dae_rng_u64(g) >> 11) / 9007199254740992.0;
}

int32_t dae_rng_below(dae_rng *g, int32_t n)
{
  uint64_t range, blocks, bound, r;
  if (n <= 1) return 0;
  range  = (uint64_t)n;
  blocks = UINT64_MAX / range;   /* quantos blocos inteiros de `range` cabem */
  bound  = blocks * range;       /* aceita só [0, bound): uniforme, sem viés  */
  do { r = dae_rng_u64(g); } while (r >= bound);
  return (int32_t)(r % range);
}
