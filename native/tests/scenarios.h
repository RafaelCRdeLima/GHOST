/* scenarios.h — cenarios compartilhados pelo teste 6 (WASM contra nativo).
 *
 * UMA definicao, dois alvos. O mesmo arquivo e compilado pelo gcc e pelo emcc,
 * e e ele que garante que os dois lados rodam a MESMA simulacao — pela mesma
 * razao que o nucleo e amalgamado em vez de reescrito.
 */
#ifndef DAE_TEST_SCENARIOS_H
#define DAE_TEST_SCENARIOS_H

#include "dae.h"

typedef struct {
  const char     *nome;
  dae_gen_params  gen;
  dae_ham_kind    ham;
  dae_norm_kind   norm;
  double          gamma;
  double          t1;
  int32_t         nt;
  int32_t         lanczos;
  int32_t         init_site;
  int32_t         target;
} dae_scenario;

int                 dae_scenario_count(void);
const dae_scenario *dae_scenario_get(int i);

#endif /* DAE_TEST_SCENARIOS_H */
