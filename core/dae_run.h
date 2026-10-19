/* dae_run.h — do spec.json à série temporal.
 *
 * É o que o navegador, o binário nativo e o .cpp exportado chamam. Reentrante e
 * sem estado: o template exportado roda realizações sob `#pragma omp parallel
 * for`, e cada thread traz o seu próprio spec, grafo e séries.
 */
#ifndef DAE_RUN_H
#define DAE_RUN_H

#include "dae_metrics.h"
#include "dae_spec.h"
#include "dae_types.h"

/* Colunas de dae_series.scal, uma linha por passo. Definidas AQUI e usadas
 * também pela ponte do WASM: duas listas de colunas em arquivos diferentes
 * seriam uma fonte de divergência silenciosa, só que em índices. */
enum { DAE_S_NORMA = 0, DAE_S_IPR, DAE_S_COH_L1, DAE_S_P_ALVO, DAE_S_NCOL };

typedef struct {
  double   escala;       /* fator de normalização de ||H|| aplicado          */
  double   lo, hi, a, b; /* intervalo espectral e a transformação de Chebyshev */
  double   dt;
  /* alpha = a*dt. É ELE que a fórmula da ordem usa, e é nele que um teto de K
     falharia de forma não-monotônica (CONVENTIONS.md, parte 6.2). A grade
     exportada tem dt livre e `a` dependendo da normalização, então alpha varre
     uma faixa larga — por isso ele sai daqui e entra no cabeçalho do CSV, em
     vez de ficar implícito. */
  double   alpha;
  int32_t  k_used;
  uint64_t fingerprint;
  int      lanczos_used;
} dae_run_info;

typedef struct {
  int32_t  nt, n, nmod, npop;
  double  *t;         /* nt                    */
  double  *scal;      /* nt * DAE_S_NCOL       */
  double  *pmod;      /* nt * nmod             */
  double  *pop;       /* npop * n  (opcional)  */
  double  *conc_mod;  /* nt * nmod^2 (opcional)*/
  double  *psi_re;    /* n: estado final. Os observáveis são funções de |psi|,
                         então um erro de FASE não apareceria em nenhum deles —
                         e é justamente o tipo de erro que o fator (-i)^k pode
                         introduzir. O teste 7 compara o estado, não só o que
                         se deriva dele. */
  double  *psi_im;
  dae_run_info info;
  /* GHOST. Presentes so quando o spec tem bloco `neutrino`; NULL fora disso.
     Cada um tem nt * n_flavors entradas.

     `nu_sabor`  P(nu_alpha) no espaco fisico, ja projetado do mergulho.
     `nu_massa`  |<nu_i|psi>|^2 — a base de MASSA. E ela, e nao P de sabor,
                 que mede se um perfil convergiu: P oscila e engana.
     `nu_inc`    P INCOERENTE por sabor, sum_i |U_ai|^2 |<nu_i|psi>|^2.
                 Para fonte astrofisica e ESTE o observavel: os pacotes dos
                 autoestados de massa se separam em ~32 km, muito antes de
                 sair da galaxia, entao a fase que define o P coerente ja nao
                 existe quando o neutrino chega. */
  double  *nu_sabor, *nu_massa, *nu_inc;
  int32_t  nu_nf;
} dae_series;

/* Devolve 0 para cancelar. */
typedef int (*dae_progress_fn)(void *user, int32_t passo, int32_t total);

/* A semente vai SEPARADA do spec, e não dentro dele, porque a varredura de
 * realizações precisa variá-la sem clonar a estrutura (que tem ponteiros de
 * heap e não sobrevive a uma cópia rasa que alguém depois libere duas vezes).
 * Execução única passa S->seed. */
dae_status dae_run(const dae_spec *S, uint64_t seed,
                   dae_graph *G, dae_metrics *M, dae_series *R,
                   dae_progress_fn cb, void *user);

void dae_series_free(dae_series *R);

/* ------------------------------------------------------- setor magnético --
 *
 * Saída de dae_numag_run: uma linha por fatia de t1/nt km ao longo da
 * trajetória, a partir de s_km[0]. Populações das duas helicidades no espaço
 * FÍSICO (já projetadas do mergulho), a norma e quantos subpassos o guarda
 * exigiu na fatia. */
typedef struct {
  int32_t  nt;
  double  *s_km;            /* nt: fim de cada fatia                       */
  double  *p_l, *p_r;       /* nt                                          */
  double  *norma;           /* nt                                          */
  int32_t *subdiv;          /* nt                                          */
  int32_t  subdiv_max;
  double   psi_re[2], psi_im[2];   /* estado final, espaço físico; só no
                                      modo coerente (tem_estado = 1)       */
  int      tem_estado;
  /* Com colisões: erro padrão de p_R (desvio da média) e a coerência média
     rho_LR = <psi_L psi_R^*>; NULL sem colisões. */
  double  *p_r_erro, *coh_re, *coh_im;
  int32_t  n_traj;
  double   saltos_medio, salto_por_km;
  dae_numag_janelas_info janelas;  /* acumulado, com integrator windows   */
} dae_numag_serie;

/* Roda o spec com blocos magnetic + trajectory. Recusa (DAE_ERR_PARAM) spec
 * sem eles. Em erro, R sai liberado. */
dae_status dae_numag_run(const dae_spec *S, dae_numag_serie *R,
                         dae_progress_fn cb, void *user);
void       dae_numag_serie_free(dae_numag_serie *R);

#endif /* DAE_RUN_H */
