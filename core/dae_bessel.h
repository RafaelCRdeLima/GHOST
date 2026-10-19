/* dae_bessel.h — J_0..J_kmax(x) num só passe.
 *
 * Recorrência REGRESSIVA de Miller com a normalização J_0 + 2*sum J_2m = 1.
 * A recorrência progressiva J_{k+1} = (2k/x) J_k - J_{k-1} é numericamente
 * instável para k > x: a solução Y_k, que cresce, contamina e destrói a cauda.
 * A cauda é exatamente onde o truncamento de Chebyshev decide o erro, então
 * usar a progressiva aqui seria trocar o critério de parada por ruído.
 *
 * Validado contra tabela de alta precisão em native/tests/t00_bessel_table.c
 * (gerada no Wolfram, conferida contra mpmath). Sem esse teste, o teste 2
 * compararia Bessel com Bessel — a mesma rotina dos dois lados da igualdade.
 */
#ifndef DAE_BESSEL_H
#define DAE_BESSEL_H

#include "dae_types.h"

/* out[k] = J_k(x), k = 0..kmax. Exige x >= 0; o chamador aplica
 * J_k(-x) = (-1)^k J_k(x). Valores cujo verdadeiro módulo esteja abaixo da
 * faixa do double saem como 0. */
dae_status dae_bessel_j_array(double x, int32_t kmax, double *out);

#endif /* DAE_BESSEL_H */
