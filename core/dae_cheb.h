/* dae_cheb.h — propagador de Chebyshev com produto matriz-vetor esparso.
 *
 *   exp(-iHt)|psi> = exp(-i b t) * sum_k c_k T_k(H~)|psi>,   H~ = (H - b)/a
 *   c_0 = J_0(a t),  c_k = 2 (-i)^k J_k(a t)
 *
 * Recorrência de três termos, sem nunca materializar os T_k. Nada de expm
 * denso: em N = 2600 seriam ~1e11 flops por chamada.
 */
#ifndef DAE_CHEB_H
#define DAE_CHEB_H

#include "dae_csr.h"
#include "dae_types.h"

/* Limiar da cauda: o último coeficiente retido fica abaixo disto. É ele, e
 * não a fórmula de dae_cheb_order, que determina a ordem efetiva. */
#define DAE_CHEB_TAIL 1e-16

typedef struct {
  int32_t  n;
  double   a, b;              /* H~ = (H - b)/a                        */
  double   lo, hi;            /* intervalo espectral efetivamente usado */
  int      lanczos_used;      /* 1 se Lanczos apertou Gershgorin        */

  double  *t0re, *t0im;       /* T_{k-1}                                */
  double  *t1re, *t1im;       /* T_k                                    */
  double  *t2re, *t2im;       /* rascunho                               */
  double  *acre, *acim;       /* acumulador da soma                     */

  double  *ckre, *ckim;       /* coeficientes complexos                 */
  double  *jbuf;              /* J_k(|a dt|)                            */
  int32_t  kcap, k_used;
  double   cached_dt;
  int      coef_valid;
} dae_cheb;

typedef struct {
  double  lo, hi;             /* intervalo espectral usado              */
  double  c_tail;             /* |c_K| no corte — diagnóstico            */
  int32_t k_used;             /* ordem efetiva do último passo          */
} dae_cheb_info;

/* Limites espectrais. Gershgorin é RIGOROSA e vira o teto; com
 * lanczos_steps > 0 e um PRNG, o intervalo é apertado por dentro dela e nunca
 * além dela — subestimar `a` põe autovalores de H~ fora de [-1,1], os T_k
 * explodem e o resultado vira lixo silencioso, enquanto superestimar só custa
 * alguns termos a mais em K. O erro é assimétrico e a escolha reflete isso.
 * O Lanczos usado aqui é determinístico e não consome PRNG nenhum. */
dae_status dae_cheb_init(dae_cheb *W, const dae_csr *H, int32_t lanczos_steps);
void       dae_cheb_free(dae_cheb *W);

/* CHUTE INICIAL da ordem, K ~ 1.2 (a |dt|) + 20. A ordem que o propagador de
 * fato usa é decidida pela cauda (|2 J_K| < DAE_CHEB_TAIL) e sai em
 * dae_cheb_info.k_used; para alpha entre ~20 e ~100 ela é MAIOR que este
 * chute. Ver a medição em dae_cheb.c. */
int32_t    dae_cheb_order(const dae_cheb *W, double dt);

/* psi <- exp(-i H dt) psi, no lugar.
 *
 * PREMISSA DA GRADE UNIFORME: os c_k dependem só de dt e ficam em cache; um dt
 * repetido custa zero. Uma grade de tempo LOGARÍTMICA anula esse cache e
 * recalcula Bessel a cada passo — funciona, mas custa.
 *
 * Devolve DAE_ERR_NORM se a norma variar mais que DAE_NORM_TOL em relação à
 * entrada — rede de segurança contra limite espectral apertado demais. */
dae_status dae_cheb_step(dae_cheb *W, const dae_csr *H, double dt,
                         double *psire, double *psiim, dae_cheb_info *info);

#endif /* DAE_CHEB_H */
