/* dae_ham.h — do grafo ao Hamiltoniano, com a escala de energia fixada.
 *
 * DISCIPLINA METODOLÓGICA (ver CONVENTIONS.md, parte 3): ao comparar
 * topologias diferentes é obrigatório normalizar ||H||, senão "mais coerência"
 * pode ser só "hopping maior". A normalização vem LIGADA por padrão; desligar
 * é escolha explícita e fica registrada no spec.json.
 *
 * Semântica da escala: o fator é medido em H_cru (gamma = 1) e depois
 * H = gamma * H_cru / escala. Assim gamma continua sendo o botão físico e a
 * normalização remove só a escala imposta pela topologia.
 */
#ifndef DAE_HAM_H
#define DAE_HAM_H

#include "dae_csr.h"
#include "dae_types.h"

typedef enum { DAE_H_ADJACENCY = 0, DAE_H_LAPLACIAN = 1 } dae_ham_kind;

typedef enum {
  DAE_NORM_NONE        = 0,
  DAE_NORM_SPECTRAL    = 1,   /* raio espectral de H_cru = 1 (via Lanczos) */
  DAE_NORM_MEAN_DEGREE = 2    /* grau médio ponderado = 1 (exato)          */
} dae_norm_kind;

/* H = -gamma A  (adjacência)  ou  H = gamma L = gamma (D - A)  (laplaciana),
 * dividido pelo fator de normalização, que sai em *scale_out para a interface
 * exibir.
 *
 * DAE_NORM_SPECTRAL estima o raio espectral por lanczos_steps passos de
 * Lanczos com vetor inicial DETERMINÍSTICO — sem dae_rng. Se puxasse do fluxo
 * compartilhado, a escala de energia dependeria da ordem das chamadas e de
 * qual semente o grafo gastou, contaminando exatamente a comparação a escala
 * fixa entre topologias que a normalização existe para garantir. Ver
 * CONVENTIONS.md, parte 4. */
dae_status dae_hamiltonian(dae_csr *H, const dae_csr *A, dae_ham_kind kind,
                           double gamma, dae_norm_kind norm,
                           int32_t lanczos_steps, double *scale_out);

#endif /* DAE_HAM_H */
