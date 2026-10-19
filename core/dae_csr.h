/* dae_csr.h — matriz esparsa real simétrica e limites espectrais. */
#ifndef DAE_CSR_H
#define DAE_CSR_H

#include "dae_types.h"

dae_status dae_csr_alloc(dae_csr *A, int32_t n, int32_t nnz);
void       dae_csr_free(dae_csr *A);

/* Lista de arestas NÃO-DIRECIONADAS → CSR simétrico.
 *
 * DEDUPLICAÇÃO: aresta repetida é DESCARTADA, nunca somada. Com seam_shift e
 * N_perp pequeno dá para gerar duplicata legítima; somar os pesos dobraria
 * j_perp em alguns sítios sem nenhum aviso. O número de entradas descartadas
 * sai em *n_dropped (pode ser NULL) para a interface avisar.
 *
 * w == NULL significa peso 1 em todas as arestas. Laço (i==i) é aceito e vira
 * entrada diagonal, que é como se representa a energia de sítio.
 */
dae_status dae_csr_from_edges(dae_csr *A, int32_t n,
                              const int32_t *ei, const int32_t *ej,
                              const double *w, int32_t ne,
                              int32_t *n_dropped);

/* y = A x, com a mesma matriz real aplicada às duas componentes. Laço quente. */
void dae_csr_spmv(const dae_csr *A,
                  const double *xre, const double *xim,
                  double *yre, double *yim);

/* y = A x, versão puramente real (Lanczos, métricas). */
void dae_csr_spmv_real(const dae_csr *A, const double *x, double *y);

/* Cota de Gershgorin: RIGOROSA, o espectro está contido em [*lo, *hi] para
 * qualquer matriz simétrica. É o teto que dae_cheb_init nunca ultrapassa. */
void dae_csr_gershgorin(const dae_csr *A, double *lo, double *hi);

/* Estimativa por m passos de Lanczos, para APERTAR a cota de Gershgorin.
 *
 * Devolve [theta_min - beta_m, theta_max + beta_m]: os valores de Ritz são
 * interiores ao espectro (theta_max <= lambda_max sempre), e beta_m = ||r|| é
 * a margem rigorosa em aritmética exata. Em ponto flutuante a ortogonalidade
 * se perde, então quem chama AINDA precisa inflar a semi-largura e recortar
 * contra Gershgorin — é o que dae_cheb_init faz.
 *
 * SEM dae_rng: o vetor inicial é DETERMINÍSTICO, gerado internamente a partir
 * de uma constante fixa. Consumir o fluxo compartilhado do PRNG faria a ordem
 * das chamadas virar parte do contrato de forma invisível — bastaria o .cpp
 * exportado invocar isto antes ou depois da religação para o grafo mudar. E
 * com DAE_NORM_SPECTRAL a escala de energia passaria a depender de uma
 * estimativa estocástica, contaminando justamente a comparação a escala fixa
 * entre topologias. Ver CONVENTIONS.md, parte 4.
 */
dae_status dae_csr_lanczos_bounds(const dae_csr *A, int32_t m,
                                  double *lo, double *hi);

/* Vetor inicial determinístico do Lanczos, exposto porque dae_metrics usa o
 * mesmo. splitmix64 sobre uma constante fixa: nada de dae_rng. */
void dae_csr_lanczos_start(double *v, int32_t n);

/* --- tridiagonal simétrica (al[0..m-1], be[0..m-2]) ---
 * Públicas porque dae_cheb e dae_metrics usam as mesmas; duplicá-las como
 * `static` em dois arquivos quebraria a amalgamação num tradutor só. */

/* Número de autovalores estritamente menores que mu (sequência de Sturm). */
int32_t dae_tridiag_count(const double *al, const double *be, int32_t m, double mu);

/* Bissecção pelo limiar `target` da contagem de Sturm dentro de [lo, hi]:
 * target = 0 dá o menor autovalor, target = m-1 dá o maior. */
double dae_tridiag_bisect(const double *al, const double *be, int32_t m,
                          int32_t target, double lo, double hi);

/* Resolve (T - sigma I) x = b pelo algoritmo de Thomas, x entra como b e sai
 * como solução. `d` é rascunho de m doubles. */
void dae_tridiag_solve(const double *al, const double *be, int32_t m,
                       double sigma, double *x, double *d);

#endif /* DAE_CSR_H */
