/* dae_neutrino.h — oscilação de sabor: PMNS, matéria, e o mergulho real.
 *
 * Um neutrino de n sabores é um sistema de n níveis com H hermitiano. O núcleo
 * propaga exp(-iHt) para H esparsa e hermitiana; oscilação de sabor é um caso
 * particular. Este módulo existe para que o `spec.json` possa DIZER que é um
 * neutrino, em vez de esconder a física numa lista de arestas com pesos
 * numéricos, o que quebraria a procedência: de uma edgelist de 6 nós ninguém
 * recupera theta13, Dm2, E ou rho.
 *
 * UNIDADES. H sai em eV. O tempo do propagador é, portanto, 1/eV, e uma linha
 * de base em km converte por DAE_NU_T_POR_KM. Ver a nota lá embaixo.
 *
 * SEM _Complex. CONVENTIONS.md parte 1 exige que o núcleo seja C99 que também
 * é C++17 válido, e `_Complex` não é. A aritmética complexa da PMNS anda em
 * pares de double (hre/him), explícita e feia, de propósito.
 */
#ifndef DAE_NEUTRINO_H
#define DAE_NEUTRINO_H

#include "dae_csr.h"
#include "dae_types.h"

/* t = L[km] * DAE_NU_T_POR_KM, com H em eV.  É 1e5 cm/km dividido por
 * hbar*c = 1.9732698e-5 eV cm. */
#define DAE_NU_T_POR_KM 5067730727.952153

/* V = sqrt(2) G_F n_e, com n_e = Y_e rho / m_N, dá V[eV] por (Y_e * rho[g/cm3]).
 * Derivado das constantes, não tabelado: sqrt(2) G_F (hbar c)^3 / m_N. */
#define DAE_NU_V_POR_RHO 7.566899624856641e-14

#define DAE_NU_MAX_SABORES 3

typedef struct {
  int32_t n_flavors;          /* 2 ou 3                                      */
  double  theta12, theta13, theta23;   /* radianos; com n_flavors=2 só t12    */
  double  delta_cp;           /* radianos. Fora de {0, pi} torna H complexa. */
  double  dm21_sq, dm31_sq;   /* eV^2. dm31_sq < 0 é hierarquia invertida.   */
  double  energy_eV;
  int     antineutrino;       /* conjuga delta_cp e inverte o sinal de V     */
} dae_nu_params;

void dae_nu_params_default(dae_nu_params *P);

/* H no espaço de SABOR, em eV, com o traço removido (uma fase global não é
 * observável, e tirá-la mantém o espectro centrado — o que ajuda Chebyshev).
 * `hre` e `him` recebem n*n cada, row-major. `rho_g_cm3` e `ye` entram no
 * potencial de matéria de Wolfenstein, que afeta só nu_e. */
dae_status dae_nu_hamiltoniano(const dae_nu_params *P, double rho_g_cm3, double ye,
                               double *hre, double *him);

/* 1 se H tem parte imaginária que importa (logo precisa do mergulho). */
int dae_nu_precisa_mergulho(const double *him, int32_t nf);

/* A adjacência REAL que o resto do pipeline consome, com H = -gamma*A,
 * gamma = 1 e normalization = none.
 *
 * Se Im(H) = 0, A é n x n e *mergulhado = 0.
 *
 * Se Im(H) != 0, monta o MERGULHO: com H = X + iY, X simétrica e Y
 * ANTIssimétrica (é o que hermiticidade obriga),
 *
 *     H~ = [ X  -Y ]     é real e SIMÉTRICA, porque (-Y)^T = Y,
 *          [ Y   X ]
 *
 * e o subespaço {(psi, -i psi)} é invariante sob H~, onde H~ age como H. Então
 * A é 2n x 2n, *mergulhado = 1, e o chamador usa dae_nu_inicial/dae_nu_projeta
 * para entrar e sair do subespaço sem saber do truque.
 *
 * LIMITE: o mergulho vale para dinâmica UNITÁRIA. Um colapso em |j><j| no
 * espaço de 2n destrói a forma (psi, -i psi); com defasagem, o colapso tem de
 * ser no PAR {j, n+j}. */
dae_status dae_nu_adjacencia(const double *hre, const double *him, int32_t nf,
                             dae_csr *A, int *mergulhado);

/* |nu_sabor> no espaço em que A vive. `re`/`im` recebem nf (ou 2*nf). */
void dae_nu_inicial(int32_t sabor, int32_t nf, int mergulhado,
                    double *re, double *im);

/* Volta do espaço de A para nf amplitudes complexas. `pre`/`pim` recebem nf. */
void dae_nu_projeta(const double *re, const double *im, int32_t nf, int mergulhado,
                    double *pre, double *pim);

/* --------------------------------------------------------- perfil de materia --
 *
 * A densidade ao longo do caminho. `constante` e o caso de laboratorio; a lei
 * de potencia e o envelope de supernova, rho ~ r^-3.
 */
typedef enum {
  DAE_NU_PERFIL_CONSTANTE = 0,
  DAE_NU_PERFIL_POTENCIA,
  DAE_NU_PERFIL_TABELA
} dae_nu_perfil_kind;

typedef struct {
  dae_nu_perfil_kind kind;
  double ye;
  double rho;             /* CONSTANTE: g/cm3                                */
  double rho_ref;         /* POTENCIA: rho(r) = rho_ref * (r_ref/r)^expoente  */
  double r_ref;           /* km                                              */
  double expoente;
  double r_start;         /* km: onde o caminho comeca                       */

  /* TABELA: um perfil estelar de verdade, ponto a ponto.
   *
   * Existe porque progenitor real NAO e lei de potencia. O que a lei de
   * potencia nao tem e o que importa aqui: o GRADIENTE varia, e e o gradiente
   * que controla a adiabaticidade — a largura da ressonancia vale
   * 2H tan(2theta) com H = |dln n_e/dr|^-1. Num progenitor de 20 Msun a
   * interface O/He e ~3x mais abrupta que o entorno, e a ressonancia H de
   * 20 MeV cai dentro da camada de oxigenio.
   *
   * `ye_tab` pode ser NULL, e ai vale `ye` constante. A interpolacao e
   * LOG-LOG, que e a natural para estrutura estelar: em log-log o perfil e
   * quase reto por trechos, e interpolar linearmente em rho daria densidade
   * errada por ordens de grandeza entre dois pontos afastados. */
  const double *r_tab;    /* km, CRESCENTE                                   */
  const double *rho_tab;  /* g/cm3                                           */
  const double *ye_tab;   /* opcional; NULL usa `ye`                         */
  int32_t       n_tab;

  /* FLUTUACAO com espectro de potencia, sobre o perfil suave.
   *
   * rho(r) -> rho(r) * [1 + sum_n a_n cos(k_n r + phi_n)]
   *
   * com k_n = 2 pi n / L, amplitudes a_n ~ n^(-p/2) normalizadas para dar
   * `rms` de contraste, e fases de dae_rng semeado por `ruido_semente`. E o
   * "randomization method" que a literatura de turbulencia usa.
   *
   * DUAS COISAS QUE ESTE BLOCO EXISTE PARA PERMITIR, e a segunda e o ponto:
   *
   * 1. amplitude e expoente ANCORADOS em medida, nao escolhidos. Supor
   *    "alguns %" onde a simulacao da 3e-4 foi o que produziu, na literatura,
   *    um amortecimento duas ordens de grandeza maior que o real.
   * 2. `n_max` CORTA o espectro. Comparar com e sem corte responde se as
   *    escalas que simulacao nenhuma resolve mudam a resposta — que e a
   *    frouxidao que a area carrega desde 2013.
   *
   * `ruido_rms <= 0` desliga tudo. */
  /* Amplitudes e fases PRE-CALCULADAS, `ruido_n_max - ruido_n_min + 1` cada.
   * NULL faz o nucleo recalcular a cada avaliacao, o que funciona e e LENTO:
   * medido, o custo por corrida vai de 9 s (n_max=100) a 250 s (n_max=1000),
   * escalando como n_max^1.4, porque `pow` e o PRNG sao refeitos em cada uma
   * das centenas de milhares de avaliacoes que o guarda de adiabaticidade
   * pede. Com os arrays, sobra so o cosseno.
   * O chamador e dono da memoria — o nucleo nao guarda estado (parte 1). */
  const double *ruido_amp, *ruido_fase;

  double  ruido_rms;      /* contraste rms de densidade                      */
  double  ruido_p;        /* expoente do espectro (5/3 = Kolmogorov)         */
  int32_t ruido_n_min;    /* menor modo (1 = escala do dominio)              */
  int32_t ruido_n_max;    /* maior modo; e ELE o objeto do teste do corte    */
  uint64_t ruido_semente;
} dae_nu_perfil;

void   dae_nu_perfil_default(dae_nu_perfil *M);
double dae_nu_densidade(const dae_nu_perfil *M, double r_km);

/* Preenche `amp` e `fase` (n_max - n_min + 1 cada) e os liga ao perfil.
 * Chamar UMA vez, depois de fixar rms/expoente/faixa/semente. */
void dae_nu_ruido_prepara(dae_nu_perfil *M, double *amp, double *fase);

/* Y_e no raio. Com tabela e `ye_tab`, interpola; senao devolve `ye` constante.
 * Importa porque V = sqrt(2) G_F n_e depende de Y_e, e num progenitor real ele
 * varia entre as camadas. */
double dae_nu_ye(const dae_nu_perfil *M, double r_km);

/* Propaga psi de r0 a r1 (km) com CONTROLE DE PASSO ADAPTATIVO.
 *
 * O problema que isto resolve: uma escada de densidade grossa demais na
 * ressonancia gera transicoes nao-adiabaticas ESPURIAS, e o resultado pode
 * sair errado por um fator 8 com a curva na forma certa.
 * Um `expm` de referencia erra igual — nao e bug de propagador, e erro de
 * discretizacao, e por isso comparar com outro codigo nao pega.
 *
 * O criterio NAO e o parametro de Landau-Zener. E comparacao direta: resolve a
 * fatia com um passo e com dois, e se as duas respostas diferirem mais que
 * `tol`, subdivide e repete. Pega qualquer erro de discretizacao, nao so o
 * nao-adiabatico, e nao depende de a teoria de adiabaticidade valer no regime.
 *
 * Devolve DAE_ERR_TOOBIG se o teto de subdivisoes for atingido — RECUSAR e o
 * ponto: um numero errado em silencio e a classe de falha que CONVENTIONS.md
 * parte 10 existe para proibir. `subdiv` (pode ser NULL) recebe quantas
 * subdivisoes foram de fato necessarias, para o diagnostico sair no CSV.
 */
dae_status dae_nu_propaga(const dae_nu_params *P, const dae_nu_perfil *M,
                          double r0_km, double r1_km, double tol,
                          int merg, double *re, double *im, int32_t *subdiv);

/* Matriz PMNS, exposta porque o chamador precisa dela para a base de MASSA —
 * que é onde a convergência de um perfil se mede, e onde a decoerência
 * astrofísica age. `ure`/`uim` recebem nf*nf. */
dae_status dae_nu_pmns(const dae_nu_params *P, double *ure, double *uim);

#endif /* DAE_NEUTRINO_H */
