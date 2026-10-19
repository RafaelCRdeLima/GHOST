/* dae_metrics.h — métricas de rede.
 *
 * O parâmetro de controle é `p`, mas os resultados são reportados contra
 * lambda_2 e Q. Isso põe lambda_2 na posição de resultado de manchete, e é
 * justamente onde a numérica é traiçoeira: rede fortemente modular tem
 * lambda_2 -> 0, e é ali que a convergência é mais lenta. Um número FIXO de
 * passos devolveria lambda_2 SUPERESTIMADO no regime modular, achatando a
 * curva na região mais interessante — e um platô numérico é indistinguível de
 * um platô físico quando se olha só a curva.
 *
 * Por isso aqui não existe "rodou m passos e pronto": existe critério de
 * convergência, teto de passos, resíduo devolvido e bandeira dizendo se
 * convergiu. Se `lambda2_converged` for 0, o valor é um LIMITE SUPERIOR, não
 * uma medida.
 */
#ifndef DAE_METRICS_H
#define DAE_METRICS_H

#include "dae_graph.h"
#include "dae_types.h"

typedef struct {
  int32_t lambda2_max_steps;  /* 0 = padrão (min(300, n), limitado por memória) */
  double  lambda2_tol;        /* 0 = padrão 1e-8, relativo a ||L||              */
  int32_t bfs_sources;        /* 0 = todas as fontes (exato)                    */
} dae_metrics_cfg;

typedef struct {
  /* valor de Fiedler da laplaciana L = D - A */
  double  lambda2;
  double  lambda2_residual;   /* ||L x - lambda2 x||, ||x||=1: COTA RIGOROSA —
                                 existe autovalor de L a menos disto          */
  int     lambda2_converged;  /* resíduo abaixo de tol*||L||                  */
  int32_t lambda2_steps;      /* passos de Lanczos gastos                     */
  double  lambda2_const_leak; /* MAIOR |<v, 1/sqrt(n)>| observado ao longo de
                                 toda a iteracao, sobre os vetores de Lanczos e
                                 os de Ritz. E a quantidade que os quatro pontos
                                 de deflacao existem para segurar, e por isso o
                                 teste olha ELA, e nao so o lambda_2 final: o
                                 sintoma aparece tarde e mascarado, o mecanismo
                                 aparece cedo. Sadio: < 1e-15.               */

  double  modularity_Q;       /* Newman, ponderada, sobre a partição do grafo  */
  double  mean_degree;
  double  mean_path_len;
  int     path_len_exact;     /* 0 = estimado por amostragem de fontes         */
  int32_t n_edges;
  int32_t n_components;
} dae_metrics;

void       dae_metrics_cfg_default(dae_metrics_cfg *c);

/* Layout espectral: autovetores 2 e 3 da laplaciana como coordenadas (x, y),
 * escritos em `xy` (2*n floats). Reaproveita o Lanczos deflacionado do
 * lambda_2 — e literalmente o vetor de Fiedler mais o proximo.
 *
 * Para grafo modular o embedding SEPARA os modulos no plano. Para grafo sem
 * estrutura ele devolve uma nuvem sem separacao, e isso e resposta, nao falha:
 * a fixture em native/tests/t95_layout.c verifica os DOIS casos, porque uma
 * asserção de separacao que passasse tambem no grafo homogeneo nao estaria
 * testando nada.
 *
 * Coordenadas CRUAS: normalizar para a caixa do canvas e de quem desenha. */
dae_status dae_layout_espectral(const dae_graph *G, float *xy);
dae_status dae_metrics_compute(const dae_graph *G, const dae_metrics_cfg *cfg,
                               dae_metrics *out);

#endif /* DAE_METRICS_H */
