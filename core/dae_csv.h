/* dae_csv.h — o formato de saída, e o único.
 *
 * O navegador, o binário nativo e o .cpp exportado emitem EXATAMENTE o mesmo
 * texto. Se cada um serializasse do seu jeito, o teste 7 estaria comparando
 * duas serializações e não duas simulações — e uma diferença de formato
 * apareceria como discordância física.
 *
 * DUAS CLASSES DE COMENTÁRIO, e a distinção importa:
 *
 *   `#!`  reprodutibilidade. Tem de bater byte a byte entre alvos: versão,
 *         hash do núcleo, impressão digital do grafo, spec.json canônico,
 *         alpha, a, b, escala. É o que o comparador confere primeiro.
 *   `#`   informativo. Data, linha de comando, quem gerou. Muda entre execuções
 *         por construção, e o comparador ignora.
 *
 * Sem essa separação, ou o CSV não carrega a data (e perde rastreabilidade) ou
 * o comparador precisa de uma lista de exceções por nome de campo, que apodrece.
 */
#ifndef DAE_CSV_H
#define DAE_CSV_H

#include "dae_metrics.h"
#include "dae_run.h"
#include "dae_spec.h"

/* Escreve em buf e devolve quantos bytes escreveria (como snprintf): chame com
 * cap = 0 para dimensionar, aloque, chame de novo. */
int32_t dae_csv(const dae_spec *S, const dae_graph *G, const dae_metrics *M,
                const dae_series *R, int incluir_estado, char *buf, int32_t cap);

/* Cabeçalho de procedência do ACERVO DE PSI do modo `archive_psi`. Mesmas
 * linhas `#!` do CSV — spec canônico, core_hash, implementacao — mais a
 * semente-base e o formato do bloco binário que vem depois.
 *
 * É o mesmo emissor de propósito: um acervo com procedência própria, escrita
 * noutro lugar, divergiria do CSV sem que nada quebrasse, e um arquivo que diz
 * de qual simulação veio é a única coisa que separa cache de lixo.
 *
 * O acervo é CACHE REPRODUZÍVEL, não acervo: apagá-lo é seguro, e regenerá-lo
 * custa tempo de CPU, não informação. */
int32_t dae_csv_cabecalho_psi(const dae_spec *S, uint64_t semente_base,
                              int32_t n, int32_t n_amostras, int32_t n_traj,
                              char *buf, int32_t cap);

/* CSV do setor magnético (dae_numag_run). Mesmas linhas de procedência `#!`
 * — versão, core_hash, implementação, spec canônico com as tabelas da
 * trajetória — mais o máximo de subpassos. Colunas: s_km (fim da fatia, no
 * comprimento local do observador estático), norm, p_L, p_R, subdiv. */
int32_t dae_numag_csv(const dae_spec *S, const dae_numag_serie *R,
                      char *buf, int32_t cap);

#endif /* DAE_CSV_H */
