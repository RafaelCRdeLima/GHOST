/* dae_numag.h — conversão de helicidade pelo momento magnético do neutrino.
 *
 * Setor magnético. Um neutrino de Dirac com momento
 * magnético mu_nu, num campo B, acopla nu_L a nu_R pelo campo TRANSVERSAL ao
 * momento:
 *
 *     H = [ Delta/2            Omega e^{+i phi_B} ]     Omega = mu_nu B_perp
 *         [ Omega e^{-i phi_B}  -Delta/2          ]
 *
 * na base (nu_eL, nu_eR), com o traço removido do ESPAÇO INTEIRO
 * (CONVENTIONS.md 13.4). Delta é o potencial relativo ativo-estéril em matéria
 * neutra de e, p, n, com fundo de nu_e aprisionados:
 *
 *     Delta = sqrt(2) G_F (n_e - n_n/2 + 2 (n_nue - n_nuebar))
 *           = DAE_NU_V_POR_RHO * rho * [(3 Y_e - 1)/2 + 2 Y_nue],
 *
 * (Nötzold e Raffelt 1988), que se anula em Y_e + (4/3) Y_nue = 1/3; sem fundo,
 * em Y_e = 1/3. Y_nue = (n_nue - n_nuebar)/n_b. Assimetrias de nu_mu e nu_tau
 * entrariam com peso 1 em vez de 2; aqui são tomadas nulas (mu = 0) e não
 * estão no meio.
 *
 * POR QUE UM MÓDULO SEPARADO de dae_neutrino. As duas físicas têm traços
 * diferentes (por sabor ativo lá, pelo espaço inteiro aqui), mergulho
 * diferente (lá só quando delta_cp exige; aqui SEMPRE, porque phi_B varia ao
 * longo do caminho e a dimensão do estado não pode mudar no meio dele) e
 * limites diferentes. Generalizar DAE_NU_MAX_SABORES mudaria uma constante que
 * o código de sabor usa em arrays de pilha.
 *
 * ESCOPO. Um sabor de NEUTRINO de Dirac (n_ativos = 1). Antineutrino,
 * Majorana e mais de um sabor devolvem DAE_ERR_UNIMPLEMENTED: o sinal do
 * acoplamento para nu-barra e os momentos de transição não estão
 * implementados. O meio também é o mínimo: rho, Y_e e B já projetado na base
 * transversal transportada, mais o fundo de nu_e. Polarizações do meio não
 * entram.
 *
 * O SINAL DE Delta_ef. Com phi_B variando, a passagem ao referencial que gira
 * com B_perp, psi = diag(e^{i phi/2}, e^{-i phi/2}) psi~, soma
 * -i U^dag dU/ds = diag(phi'/2, -phi'/2) a H, e a ressonância passa a ser o
 * zero de
 *
 *     Delta_ef [eV] = Delta [eV] + (dphi_B/ds)[rad/km] / DAE_NU_T_POR_KM.
 *
 * O sinal é CONSEQUÊNCIA da convenção H_LR = Omega e^{+i phi_B}, e vale só
 * com ela: com H_LR = Omega e^{-i phi_B} seria Delta - phi'/T. t107 o
 * confirma com um campo girante em ressonância sem zero de Delta
 * (CONVENTIONS.md 13.3). phi_B é medido na base transportada (dae_numag_meio),
 * e a rotação dessa base, se houvesse, entraria aqui do mesmo jeito.
 */
#ifndef DAE_NUMAG_H
#define DAE_NUMAG_H

#include "dae_csr.h"
#include "dae_neutrino.h"
#include "dae_types.h"

#define DAE_NUMAG_MAX_ATIVOS  3
#define DAE_NUMAG_MAX_ESTADOS (2 * DAE_NUMAG_MAX_ATIVOS)

/* Magneton de Bohr em eV/T (CODATA 2018), e o gauss em tesla. */
#define DAE_NUMAG_MU_B_EV_POR_T 5.7883818060e-5
#define DAE_NUMAG_T_POR_GAUSS   1.0e-4

typedef enum { DAE_NUMAG_DIRAC = 0, DAE_NUMAG_MAJORANA } dae_numag_natureza;

typedef struct {
  dae_nu_params      sabor;     /* energia e antineutrino; com n_ativos = 1 a
                                   massa é a mesma nas duas helicidades e a
                                   mistura não entra                         */
  dae_numag_natureza natureza;
  int32_t            n_ativos;  /* só 1 está implementado                    */
  double mu_re[DAE_NUMAG_MAX_ATIVOS * DAE_NUMAG_MAX_ATIVOS];  /* em mu_B    */
  double mu_im[DAE_NUMAG_MAX_ATIVOS * DAE_NUMAG_MAX_ATIVOS];
} dae_numag_params;

/* O meio num ponto do caminho, no referencial ortonormal do observador
 * ESTÁTICO local.
 *
 * O CAMPO ENTRA JÁ PROJETADO. `b_t` são as componentes de B na base
 * transversal (e1, e2) TRANSPORTADA ao longo do raio — e é ela que define
 * phi_B = atan2(b_t[1], b_t[0]). O núcleo não constrói base nenhuma, e de
 * propósito: uma base reconstruída em cada ponto por uma regra local (por
 * exemplo z x k) gira em torno de k quando o raio se curva, e essa rotação
 * entra em Delta_ef como um dphi_B/ds que não existe — uma ressonância
 * artificial. Transportar a base é trabalho da geometria
 * (tools/neutrinos/trajetoria_ns.py), que fica auditável no spec.
 * CONVENTIONS.md 13.2. */
typedef struct {
  double rho;              /* g/cm3                                         */
  double ye;
  double b_t[2];           /* gauss, na base transversal transportada       */
  double razao_energia;    /* E_local / E_infinito = e^{-nu/2}. Guardado; um
                              sabor sem termo de vácuo não depende de E, e
                              ele só entraria em H com mistura de sabores  */
  double ynu;              /* Y_nue = (n_nue - n_nuebar)/n_b; 0 sem fundo.
                              POR ÚLTIMO: o núcleo zera o meio antes de chamar
                              o callback, e um callback que não o preenche
                              vale "sem fundo"                             */
} dae_numag_meio;

/* Dirac, um sabor, mu = 0, 10 MeV. */
void dae_numag_params_default(dae_numag_params *P);

/* Base transversa para uma direção k FIXA: e1, e2 ortonormais e ortogonais a
 * k, com e2 = k x e1; e1 = z x k normalizado, ou x x k se k for quase z.
 *
 * SÓ PARA RAIO RETO. Com k variando ao longo do caminho, esta regra produz
 * uma base que gira em torno de k e cria um dphi_B/ds espúrio. Para raio
 * curvo, a base vem transportada da geometria; ver dae_numag_meio. */
dae_status dae_numag_base_transversa(const double k[3], double e1[3], double e2[3]);

/* Projeta B (3 componentes) na base (e1, e2): b_t[0] = B.e1, b_t[1] = B.e2. */
void dae_numag_projeta_campo(const double b[3], const double e1[3],
                             const double e2[3], double b_t[2]);

/* B_perp (gauss) e phi_B (rad, em (-pi, pi]) do meio. */
dae_status dae_numag_acoplamento(const dae_numag_meio *m,
                                 double *b_perp_gauss, double *phi_rad);

/* Delta em eV: potencial relativo nu_L - nu_R. */
double dae_numag_desajuste(const dae_numag_meio *m);

/* H (n x n, n = 2 n_ativos) em eV, row-major, traço removido do espaço
 * inteiro. */
dae_status dae_numag_hamiltoniano(const dae_numag_params *P,
                                  const dae_numag_meio *m,
                                  double *hre, double *him);

/* A adjacência real que o pipeline consome (H = -A), SEMPRE mergulhada:
 * 2n x 2n, seja Im(H) nulo ou não. Ver o cabeçalho, "por que um módulo
 * separado". */
dae_status dae_numag_adjacencia(const double *hre, const double *him, int32_t n,
                                dae_csr *A);

/* O meio ao longo do caminho, em s (km). O chamador é dono do que `user`
 * aponta; o núcleo não guarda estado. */
typedef void (*dae_numag_meio_fn)(void *user, double s_km, dae_numag_meio *m);

/* ---------------------------------------------------- trajetória tabelada --
 *
 * O meio ao longo de um raio já integrado pela geometria, ponto a ponto.
 *
 * `s_km` é o COMPRIMENTO ESPACIAL LOCAL medido pelo observador estático ao
 * longo do raio nulo, dl^2 = e^lambda dr^2 + r^2 dOmega^2 — não um
 * parâmetro afim nem a coordenada r, e não o intervalo próprio da geodésica,
 * que é nulo. É nesse comprimento que H (em eV, no referencial do observador
 * estático) gera a evolução: i dpsi/dl = H psi, com t = l * DAE_NU_T_POR_KM.
 *
 * Interpolação: rho em log (linear em ln rho, que é como perfil estelar varia;
 * cai para linear se algum ponto não for positivo); Y_e, b_t e
 * razao_energia e Y_nue lineares. Interpolar as COMPONENTES de b_t, e não phi_B, evita
 * o salto de 2 pi. Fora da tabela, segura o valor da borda. */
typedef struct {
  const double *s_km;          /* estritamente crescente                     */
  const double *rho, *ye;
  const double *b1, *b2;       /* gauss, base transportada                   */
  const double *razao_energia; /* opcional: NULL vale 1                      */
  const double *ynu;           /* opcional: NULL vale 0 (sem fundo de nu_e)  */
  int32_t       n;
} dae_numag_trajetoria;

/* Callback de meio para uma trajetória tabelada: `user` aponta um
 * dae_numag_trajetoria. */
void dae_numag_meio_tabela(void *user, double s_km, dae_numag_meio *m);

/* dae_numag_propaga sobre uma trajetória tabelada, cortando [s0, s1] nos nós
 * internos da tabela. O meio interpolado tem QUINA nos nós (derivadas
 * descontínuas); dentro de um passo, ela derruba a regra do ponto médio para
 * 1a ordem e o controle de passo subdivide até o teto. `subdiv` recebe o
 * máximo sobre os pedaços. */
dae_status dae_numag_propaga_tabela(const dae_numag_params *P,
                                    const dae_numag_trajetoria *T,
                                    double s0_km, double s1_km, double tol,
                                    double *re, double *im, int32_t *subdiv);

/* Propaga o estado MERGULHADO (2n componentes; ver dae_nu_inicial e
 * dae_nu_projeta com mergulhado = 1) de s0 a s1, com o mesmo controle de passo
 * de dae_nu_propaga: piso de factibilidade e comparação de um passo contra dois.
 *
 * O piso lê o meio pela MESMA função `meio` que a propagação usa
 * (CONVENTIONS.md 13.1), nos dois extremos da fatia. Isso basta quando |H|
 * é monótono na fatia, que é o caso dos perfis simples (camada, rampa linear,
 * campo que gira com |B| fixo); um perfil com máximo interior de |H| precisa
 * de fatias menores, e o critério de erro as pede.
 *
 * Devolve DAE_ERR_TOOBIG se o teto de subdivisões for atingido. `subdiv` (pode
 * ser NULL) recebe quantos subpassos a fatia exigiu. */
dae_status dae_numag_propaga(const dae_numag_params *P, dae_numag_meio_fn meio,
                             void *user, double s0_km, double s1_km, double tol,
                             double *re, double *im, int32_t *subdiv);

/* ------------------------------------------------ integrador por janelas --
 *
 * Propaga ao longo de uma trajetória TABELADA resolvendo por força bruta só
 * onde é preciso, e adiabaticamente no resto.
 *
 * O REFERENCIAL. Todo trecho adiabático é tratado no referencial que gira com
 * B_perp, psi = diag(e^{i phi/2}, e^{-i phi/2}) psi~, onde
 *     H~ = [ Delta_ef/2  Omega ; Omega  -Delta_ef/2 ]  é REAL,
 *     Delta_ef = Delta + (dphi_B/ds)/T.
 * A ressonância é o zero de Delta_ef, não de Delta (CONVENTIONS.md 13.3). Os
 * autoestados de H~ são reais e não têm fase de Berry; a de H no laboratório
 * fica inteira na transformação R(s).
 *
 * O CRITÉRIO NÃO É |Delta| GRANDE. Fora da janela, o erro é a amplitude de
 * transição a = int theta~' e^{i int E~} ds. Com eps~ = theta~'/(T E~) (a
 * adiabaticidade local, que já inclui Omega' — |B| variando — e phi'' —
 * orientação variando — por meio de Delta_ef'), integrar por partes dá
 *     |a| <= |eps~(início)| + |eps~(fim)| + variação total de eps~
 *            + soma dos saltos de theta~ nos nós da tabela.
 * O código AMOSTRA essa variação (extremos e meio de cada célula): o que ele
 * calcula é um ESTIMADOR da cota, em geral um limite inferior dela, não uma
 * cota certificada (CONVENTIONS.md 13.5b). Uma oscilação de theta~' em
 * ressonância com o gap, resolvida pelas amostras, acumula variação e é pega
 * mesmo com eps~ local pequeno. Trecho cujo estimador passa de `eps_tol` vira
 * janela. É de primeira ordem em eps~; o teste de cada configuração nova é
 * contra a força bruta, não contra o estimador.
 *
 * JANELA (força bruta, dae_numag_propaga): |Delta_ef| < K 2 Omega, ou
 * eps~ > eps_max, ou Delta_ef troca de sinal na célula, ou B_perp = 0.
 *
 * TRECHO ADIABÁTICO: as amplitudes nos autoestados de H~ ganham as fases
 * -+ T int E~/2 ds (quadratura de Gauss-Legendre adaptativa, até `tol_fase`
 * rad). Nos nós da tabela o meio interpolado tem quina (phi' e Delta_ef
 * saltam): psi é contínuo e a base muda de uma vez, o que é EXATO para o meio
 * interpolado. A soma desses saltos sai em `saltos_tabela`: é o erro de
 * resolução da tabela, não do integrador.
 *
 * A AMOSTRAGEM. A classificação avalia o meio em `sub` pontos por intervalo da
 * tabela (mais o ponto médio de cada célula). Estrutura mais fina que isso não
 * é vista. Dentro de um intervalo o meio interpolado é suave, mas isso NÃO
 * garante que `sub` o resolva: phi' = N/|b|^2 tem um pico de largura
 * |b|_min/|b'| quando b passa perto de zero numa célula, e pode esconder dois
 * zeros de Delta_ef entre amostras. Tabelas de campo devem pôr nós nesses
 * mínimos de |b| (CONVENTIONS.md 13.5b). */
typedef struct {
  double  K;           /* janela por |Delta_ef| < K 2 Omega                 */
  double  eps_max;     /* janela onde eps~ local passa disto                */
  double  eps_tol;     /* limite do estimador de amplitude por trecho       */
  double  tol;         /* força bruta dentro da janela                       */
  double  tol_fase;    /* quadratura da fase, rad                            */
  int32_t sub;         /* pontos de classificação por intervalo da tabela    */
} dae_numag_janelas_cfg;

/* Diagnóstico, ACUMULADO sobre as chamadas (zere antes da primeira). */
typedef struct {
  int32_t n_janelas, n_trechos, n_cruzamentos, subdiv_max;
  double  comprimento_janela_km;
  double  eps_fora_max;      /* maior eps~ num trecho tratado como adiabático */
  double  cota_adiabatica;   /* soma dos estimadores aceitos (não certificada) */
  double  saltos_tabela;     /* soma |salto de theta~| nos nós aceitos        */
  double  cruzamento_s_km;   /* posição do último zero de Delta_ef achado     */
  /* Zeros que NÃO são cruzamentos isolados (classificação pelos sinais dos
     dois lados; um zero num ponto de amostragem é localizado nele, sem
     bisseção). Contato: Delta_ef toca zero sem trocar de sinal. Intervalo
     nulo: células com Delta_ef = 0 nos três pontos — ressonância extensa,
     não uma sequência de cruzamentos. A classificação usa os lados da
     TRAJETÓRIA, não os da chamada: várias chamadas contíguas (as fatias do
     CLI) contam o mesmo que uma. Ver CONVENTIONS.md 13.5b. */
  int32_t n_contatos;
  int32_t n_intervalos_nulos;
  double  comprimento_nulo_km;
  int     terminou_em_janela;/* a ÚLTIMA chamada terminou dentro de uma janela */
} dae_numag_janelas_info;

void dae_numag_janelas_default(dae_numag_janelas_cfg *c);

/* Modo COERENTE: re/im é o estado mergulhado (4 componentes), como em
 * dae_numag_propaga. */
dae_status dae_numag_propaga_janelas(const dae_numag_params *P,
                                     const dae_numag_trajetoria *T,
                                     double s0_km, double s1_km,
                                     const dae_numag_janelas_cfg *cfg,
                                     double *re, double *im,
                                     dae_numag_janelas_info *info);

/* Uma janela em andamento entre chamadas do modo incoerente: o operador
 * PARCIAL (laboratório, colunas = helicidades) desde o início dela, e a base
 * e as populações nesse início. Zere antes da primeira chamada. */
typedef struct {
  int    aberta;
  double ur[4], ui[4];
  double theta_a, phi_a;
  double pop_a[2];
} dae_numag_janela_aberta;

/* Modo INCOERENTE: `pop` são as populações dos autoestados de H~ (+, -) em
 * s0, e saem em s1. Nas janelas entram as probabilidades |U|^2 da força bruta;
 * entre janelas, a fase relativa é descartada. É uma escolha física, não
 * uma aproximação numérica.
 *
 * ENCADEAR chamadas por um ponto DENTRO de uma janela exige `ab`: sem ele,
 * cada chamada fecharia a janela em |M|^2, e dois pedaços de uma ressonância
 * combinados por probabilidade descartam a coerência DENTRO dela — resultado
 * errado, medido (P_R 0,839 contra 0,914). Com `ab`, uma janela que termina em
 * s1 fica ABERTA: o operador parcial continua na chamada seguinte, e `pop`
 * sai com as populações provisórias em s1 (para a linha de saída). Com
 * ab = NULL, a janela final é fechada em s1 — certo só para uma chamada única
 * sobre o trecho inteiro. */
dae_status dae_numag_propaga_janelas_inc(const dae_numag_params *P,
                                         const dae_numag_trajetoria *T,
                                         double s0_km, double s1_km,
                                         const dae_numag_janelas_cfg *cfg,
                                         double *pop,
                                         dae_numag_janela_aberta *ab,
                                         dae_numag_janelas_info *info);

/* theta~ (rad) e phi_B (rad) em s, no intervalo da tabela à DIREITA de s (o
 * de s0 num trecho que começa em s). Para converter helicidade <-> populações
 * adiabáticas: |+> = (cos th, sin th), |-> = (-sin th, cos th) no referencial
 * girante. */
dae_status dae_numag_angulos(const dae_numag_params *P, const dae_numag_trajetoria *T,
                             double s_km, int lado_esquerdo,
                             double *theta, double *phi);

/* ---------------------------------------------------- colisões e coerência --
 *
 * O OPERADOR FÍSICO. As colisões atuam só sobre nu_L (nu_R é estéril):
 *     L = sqrt(Gamma) |L><L|,
 * com Gamma a taxa TOTAL de interação do nu_L (espalhamento e absorção: as
 * duas tiram o neutrino do feixe coerente desta trajetória). Para dois níveis
 * ele deixa as populações e amortece a coerência:
 *     d rho_LR / dt = -i[H, rho]_LR - (Gamma/2) rho_LR.
 *
 * O DESDOBRAMENTO USADO. A defasagem simétrica do GHOST (dae_traj), com
 * L_j = sqrt(gamma) P_j sobre OS DOIS estados físicos, tem taxa total de salto
 * gamma INDEPENDENTE do estado e amortece a coerência a gamma. As duas dão o
 * MESMO Liouvilliano de dois níveis se
 *     gamma = Gamma / 2.
 * As trajetórias individuais diferem (no desdobramento simétrico, uma
 * "colisão" em R é um artefato da escolha); a média, não. Em unidades do
 * caminho, gamma [1/km] = Gamma [eV] * DAE_NU_T_POR_KM / 2
 * (dae_numag_salto_por_km).
 *
 * O PAR {j, n+j}. No mergulho (psi, -i psi), o projetor sobre o estado FÍSICO
 * j atua no par de componentes {j, n+j}; colapsar num componente só destruiria
 * a forma do mergulho. O par é a representação de P_j — não é o par nu_L/nu_R.
 * Aqui o salto é feito no espaço físico e o estado é remergulhado.
 *
 * A TAXA AO LONGO DO CAMINHO. Gamma(s) vem por nó da trajetória tabelada
 * (linear entre nós), ou constante. O próximo salto é o ponto onde a taxa
 * acumulada int gamma ds atinge -ln u: exato para gamma linear por intervalo.
 *
 * REPRODUTIBILIDADE. Cada trajetória tem semente própria,
 * dae_traj_semente(base, índice), e as somas do ensemble entram na ordem do
 * índice (CONVENTIONS.md 11.6). */
typedef struct {
  double        gamma0_km;    /* taxa de salto do desdobramento simétrico, 1/km */
  const double *gamma_km;     /* opcional, por nó da trajetória; NULL: gamma0  */
} dae_numag_colisao;

/* gamma [1/km] = Gamma [eV] * DAE_NU_T_POR_KM / 2 */
double dae_numag_salto_por_km(double gamma_colisao_eV);

/* Média do ensemble em `nt` pontos igualmente espaçados de (s0, s1]. Cada
 * array tem nt entradas; os que forem NULL não são preenchidos.
 *   p_r       média de |psi_R|^2
 *   p_r_var   variância AMOSTRAL de |psi_R|^2 (divida por n_traj para a da
 *             média)
 *   coh_re/im média de psi_L psi_R^* (a coerência rho_LR)
 * `pr0/pi0` é o estado físico inicial (2 complexos, normalizado).
 * `saltos_medio` (pode ser NULL) recebe o número médio de saltos. */
dae_status dae_numag_ensemble(const dae_numag_params *P,
                              const dae_numag_trajetoria *T,
                              double s0_km, double s1_km,
                              const dae_numag_colisao *col, double tol,
                              int32_t n_traj, uint64_t semente_base,
                              const double *pr0, const double *pi0,
                              int32_t nt, double *p_r, double *p_r_var,
                              double *coh_re, double *coh_im,
                              double *saltos_medio);

#endif /* DAE_NUMAG_H */
