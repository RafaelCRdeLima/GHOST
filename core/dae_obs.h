/* dae_obs.h — observáveis O(N) por passo, estado puro no setor de uma excitação.
 *
 * CONVENÇÃO DA CONCURRENCE (precisa bater com C_l1, e bate — ver abaixo).
 * Para estado puro, C_ij = 2|psi_i||psi_j|. Escrevendo s_M = soma de |psi_i|
 * sobre i no módulo M e q_M = soma de |psi_i|^2 no mesmo módulo, a soma de
 * C_ij sobre PARES NÃO-ORDENADOS {i,j}, i != j, dá
 *
 *     C_MN = 2 s_M s_N          (M != N)
 *     C_MM = s_M^2 - q_M        (o produto externo inclui i == j, que não é
 *                                par: o termo espúrio é exatamente q_M)
 *
 * A matriz é armazenada simétrica, e a soma do triângulo superior COM a
 * diagonal reproduz C_l1 = (soma |psi_j|)^2 - 1. Isso custa O(N + M^2), não
 * O(N^2): a versão agregada por módulo é um produto externo, e só a matriz
 * completa é cara.
 */
#ifndef DAE_OBS_H
#define DAE_OBS_H

#include "dae_types.h"

typedef struct {
  const int32_t *module_of;   /* n entradas em [0, nmod), ou NULL */
  int32_t nmod;
  int32_t target;             /* sítio-alvo; -1 = nenhum          */
  int     want_pop;
  int     want_conc_mod;
  int     want_conc_full;     /* exige n <= DAE_FULL_CONC_MAX     */
} dae_obs_cfg;

typedef struct {
  double  norm;       /* soma |psi_j|^2 — diagnóstico, deve ficar em 1        */
  double  ipr;        /* soma |psi_j|^4                                       */
  double  coh_l1;     /* (soma |psi_j|)^2 - 1                                 */
  double  p_target;   /* |psi_alvo|^2, ou NaN quando target < 0 — nunca zero, */
                      /* que é valor fisicamente válido e viraria dado falso  */
  double *pop;        /* n     */
  double *pmod;       /* nmod  */
  double *l1mod;      /* nmod: s_M = soma |psi_i| no módulo                   */
  double *conc_mod;   /* nmod^2, opcional */
  double *conc_full;  /* n^2, opcional    */
  int32_t n, nmod;
} dae_obs;

dae_status dae_obs_alloc(dae_obs *O, const dae_obs_cfg *c, int32_t n);
void       dae_obs_free(dae_obs *O);
dae_status dae_obs_eval(const double *psire, const double *psiim, int32_t n,
                        const dae_obs_cfg *c, dae_obs *O);

#endif /* DAE_OBS_H */
