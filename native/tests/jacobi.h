/* jacobi.h — autossolver simétrico denso.
 *
 * Existe SÓ para o teste 5, como oráculo independente do Chebyshev. Não é do
 * núcleo e nunca deve entrar nele: é O(N^3) e o ponto do núcleo é não fazer
 * isso. */
#ifndef DAE_TEST_JACOBI_H
#define DAE_TEST_JACOBI_H

#include <stdint.h>

/* a: n x n simétrica, linha-maior, DESTRUÍDA. Ao voltar, w[] tem os
 * autovalores e v[] os autovetores em COLUNAS (v[i*n+k] = componente i do
 * autovetor k). Devolve 0 se convergiu. */
int dae_jacobi(double *a, double *v, double *w, int32_t n);

#endif /* DAE_TEST_JACOBI_H */
