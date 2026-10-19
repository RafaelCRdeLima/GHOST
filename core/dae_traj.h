/* dae_traj.h — ensemble de trajetórias com defasagem pura, e os DOIS modos de
 * saída.
 *
 * ---------------------------------------------------------------------------
 * POR QUE HÁ DOIS MODOS, E POR QUE NENHUM É PADRÃO
 *
 * O observável central é rho_ij = <psi_i psi_j*> promediado sobre
 * trajetórias, e só DEPOIS tomado em módulo. Calcular a coerência por trajetória
 * e promediar mediria espalhamento de amplitude, não coerência: cada trajetória
 * permanece pura.
 *
 *   ACUMULAR_RHO  soma rho durante a execução e descarta cada psi no fim da
 *                 trajetória. Memória: n_amostras * N^2 * 16 B por thread.
 *                 Escolha para varredura ampla, com os observáveis já definidos.
 *
 *   ARQUIVAR_PSI  entrega cada psi ao chamador, que grava; rho é montado na
 *                 análise. Disco: n_amostras * N * 16 B * n_traj por célula.
 *                 Escolha quando a definição do observável ainda puder mudar.
 *
 * Não há padrão implícito: `DAE_SAIDA_INDEFINIDA` é ERRO quando há trajetórias.
 * Um padrão silencioso escolheria por quem não sabe que estava escolhendo, e a
 * escolha errada só aparece meses depois, quando o observável mudar e os psi
 * não existirem mais.
 *
 * ---------------------------------------------------------------------------
 * A ORDEM DA SOMA É CONTRATO, NÃO CONSEQUÊNCIA DO LAÇO
 *
 * Soma de ponto flutuante não é associativa. Para os dois modos darem o MESMO
 * resultado bit a bit, as contribuições têm de entrar na mesma ordem — e essa
 * ordem é o índice crescente da trajetória, exigido por `dae_rho_acc_somar`,
 * que devolve DAE_ERR_PARAM se o índice não for o esperado. Não é comentário
 * pedindo boa-fé: quem paralelizar dentro de uma célula e somar fora de ordem
 * recebe erro, em vez de receber outro número.
 *
 * POR ISSO A PARALELIZAÇÃO É POR CÉLULA DA GRADE, não por trajetória dentro de
 * uma célula. Acumular por thread e reduzir em ordem de thread pareceria
 * resolver, e não resolve: o número de threads entraria no resultado, e o mesmo
 * spec com OMP_NUM_THREADS diferente daria números diferentes. Com uma célula
 * por thread, cada acumulador é privado, a soma interna é serial e ordenada, e
 * a reprodutibilidade não depende do escalonador.
 *
 * ---------------------------------------------------------------------------
 * O DESDOBRAMENTO
 *
 * Defasagem pura, L_j = sqrt(gamma)|j><j|. Como sum_j L_j^dag L_j = gamma * I, a
 * taxa total de salto NÃO depende do estado: os tempos de espera são
 * exponenciais de taxa gamma, a evolução entre saltos é unitária, e o salto
 * projeta em |j> com probabilidade |psi_j|^2. Não há aproximação de grade aqui —
 * os saltos caem em tempos arbitrários e o propagador aceita dt arbitrário.
 *
 * Isto NÃO dispensa a companheira de verificação: a média de trajetórias tem de
 * reproduzir o Liouvilliano vetorizado em N pequeno antes que qualquer física
 * seja afirmada.
 */
#ifndef DAE_TRAJ_H
#define DAE_TRAJ_H

#include "dae_cheb.h"
#include "dae_csr.h"
#include "dae_types.h"

typedef enum {
  DAE_SAIDA_INDEFINIDA = 0,   /* nao ha padrao implicito aceitavel */
  DAE_SAIDA_ACUMULAR_RHO,
  DAE_SAIDA_ARQUIVAR_PSI
} dae_saida_traj;

/* Acumulador de rho. Uma matriz N x N complexa por amostra de tempo. */
typedef struct {
  int32_t  n, n_amostras;
  double  *re, *im;      /* n_amostras * n * n, ordem linha-maior */
  /* Segundo momento, REAL: m2_ij = media_k |psi_i^k|^2 |psi_j^k|^2. Existe so
     quando o acumulador foi criado por `dae_rho_acc_init_wk`, e serve para uma
     coisa: tornar |rho_ij|^2 NAO ENVIESADO. Ver o cabecalho de dae_traj.h. */
  double  *m2;
  int32_t  proxima;      /* indice da proxima trajetoria aceita   */
  int32_t  somadas;
  int32_t  n_traj_alvo;  /* informativo; nao usado pelo acumulador */
} dae_rho_acc;

dae_status dae_rho_acc_init(dae_rho_acc *A, int32_t n, int32_t n_amostras);

/* Como o anterior, mas acumula tambem o segundo momento, o que permite a
 * correcao de vies de Wardle-Kronberg. Custa mais uma matriz N x N REAL por
 * amostra — metade da memoria de uma complexa. */
dae_status dae_rho_acc_init_wk(dae_rho_acc *A, int32_t n, int32_t n_amostras);

/* |rho_ij|^2 sem vies, pela U-estatistica de grau 2:
 *
 *     |rho_ij|^2_U = ( n |rho^_ij|^2 - m2_ij ) / (n - 1)
 *
 * Exatamente nao enviesado: E[n|rho^|^2] = (n-1)|rho|^2 + E|X|^2, e m2 estima
 * E|X|^2. E a mesma conta que a radioastronomia chama de correcao de
 * Wardle-Kronberg (subtraia a variancia do quadrado antes de tirar a raiz) e
 * que a literatura de medidas aleatorizadas chama de U-estatistica.
 *
 * Devolve 0 quando a estimativa sai negativa — o que acontece justamente onde
 * o valor verdadeiro e ~0. Esse corte reintroduz um vies positivo pequeno, e e
 * o vies residual que Simmons & Stewart (1985) mostraram que NENHUM metodo
 * conhecido elimina em SNR baixo. Nao se finge que ele nao existe: a
 * extrapolacao em n existe para medi-lo. */
double     dae_rho_acc_mod2_sem_vies(const dae_rho_acc *A, int32_t amostra,
                                     int32_t i, int32_t j, int32_t n_traj);
void       dae_rho_acc_free(dae_rho_acc *A);

/* Soma psi psi^dag da trajetoria `idx`. EXIGE idx == A->proxima. */
dae_status dae_rho_acc_somar(dae_rho_acc *A, int32_t idx,
                             const double *re, const double *im);

/* Divide por `somadas`, transformando a soma em media. Idempotente nao e:
 * chamar duas vezes divide duas vezes. */
void       dae_rho_acc_finalizar(dae_rho_acc *A);

/* Bytes de UMA trajetoria: n_amostras * n amplitudes, re e im intercalados por
 * vetor (primeiro todos os re da amostra 0, depois todos os im, e assim por
 * diante). O chamador grava; o nucleo nao abre arquivo (dae_types.h, regra 3). */
typedef dae_status (*dae_traj_sink)(void *user, int32_t idx,
                                    const double *re, const double *im,
                                    int32_t n_amostras, int32_t n);

typedef struct {
  double   gamma_deph;    /* taxa de defasagem                        */
  int32_t  n_traj;        /* numero de trajetorias                    */
  int32_t  rho_stride;    /* amostra o tempo a cada quantos passos    */
  int32_t  nt;            /* passos da grade                          */
  double   dt;
  int32_t  sitio_inicial;
  dae_saida_traj saida;

  /* --- GHOST: tres extensoes, e as tres tem a mesma origem --- */

  /* ESTADO INICIAL VETORIAL. NULL usa `sitio_inicial`. Existe porque com o
     mergulho 2n x 2n o estado inicial e (psi, -i psi)/sqrt(2), que nao e um
     sitio. Sem isto, mergulho e trajetorias sao incompativeis. */
  const double *init_re, *init_im;

  /* COLAPSO NO PAR. 0 = colapso em sitio, o caso generico.
     > 0 = o estado esta MERGULHADO, n == 2*nf_par, e o colapso e no PAR
     {j, nf_par + j}.
     A razao e geometrica: projetar num sitio SO do espaco de 2n destroi a
     forma (psi, -i psi) e joga o estado para fora do subespaco invariante —
     o resultado deixa de significar coisa alguma. Projetar no par mantem
       phi' = (psi_j e_j, -i psi_j e_j)/sqrt(2),
     que E da forma certa. E a probabilidade sai correta de graca:
       |phi_j|^2 + |phi_{n+j}|^2 = |psi_j|^2.
     Com isso, defasagem na base de SABOR funciona com fase CP ligada. */
  int32_t  nf_par;

  /* DEFASAGEM DEPENDENTE DO TEMPO. NULL usa `gamma_deph` constante.
     `gamma_tab[k]` vale no passo k, e tem `nt` entradas.
     Isto e barato, e a razao e estrutural: com gamma uniforme nos sitios a
     taxa TOTAL de salto vale gamma(t) * soma|psi_j|^2 = gamma(t), independente
     do estado, e o tempo de espera continua saindo de forma fechada por
     integracao. Uma taxa por SITIO (gamma_j) tornaria a taxa total dependente
     do estado e exigiria amostragem por rejeicao — outra coisa. */
  const double *gamma_tab;
} dae_traj_cfg;

/* Zera a config e poe os padroes.
 *
 * USE SEMPRE. A struct tem ponteiros (init_re, init_im, gamma_tab), e uma
 * `dae_traj_cfg` declarada na pilha e preenchida campo a campo deixa os que o
 * chamador nao conhece com LIXO — um ponteiro de lixo em gamma_tab e
 * segmentation fault na primeira trajetoria. */
void       dae_traj_cfg_default(dae_traj_cfg *cfg);

/* Quantas amostras de tempo a configuracao produz. */
int32_t    dae_traj_amostras(const dae_traj_cfg *cfg);

/* Roda o ensemble. Em ACUMULAR_RHO, `A` recebe as contribuicoes e `sink` e
 * ignorado; em ARQUIVAR_PSI, `sink` recebe cada trajetoria e `A` e ignorado.
 * A DINAMICA E A MESMA nos dois: o modo escolhe o destino de psi, nunca o
 * caminho que psi percorreu. E o que torna a igualdade bit a bit possivel. */
dae_status dae_traj_ensemble(const dae_traj_cfg *cfg, uint64_t semente_base,
                             dae_cheb *W, const dae_csr *H,
                             dae_rho_acc *A, dae_traj_sink sink, void *user);

/* Como `dae_traj_ensemble`, mas alimenta DOIS acumuladores: `A` recebe todas as
 * trajetorias e `Aq` recebe as primeiras `corte`. As duas estimativas saem do
 * MESMO fluxo de trajetorias, o que e o ponto: elas ficam correlacionadas, e a
 * combinacao 2*C(n) - C(corte) que cancela o vies tem variancia menor do que
 * teria com ensembles separados. */
dae_status dae_traj_ensemble_dupla(const dae_traj_cfg *cfg, uint64_t semente_base,
                                   dae_cheb *W, const dae_csr *H,
                                   dae_rho_acc *A, dae_rho_acc *Aq, int32_t corte);

/* Semente da trajetoria `idx`, derivada da base. Exposta porque a
 * reprodutibilidade sob OpenMP depende de cada trajetoria ter fluxo PROPRIO e
 * derivavel do indice — nunca de um gerador compartilhado, cujo consumo
 * dependeria da ordem em que as threads chegam. */
uint64_t   dae_traj_semente(uint64_t base, int32_t idx);

#endif /* DAE_TRAJ_H */
