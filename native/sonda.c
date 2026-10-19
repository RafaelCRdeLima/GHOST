#define _POSIX_C_SOURCE 199309L
/* sonda.c — DIMENSIONAMENTO da varredura de trajetórias, e não a varredura.
 *
 * Mede o que uma célula custa e o que ela entrega, para que a grade seja
 * escolhida com número e não com esperança. O erro a evitar: comprometer dias
 * de CPU numa grade cuja barra de erro não distingue a crista do platô. Se
 * C_inter variar 15% entre células vizinhas e o desvio do estimador for 12%, o
 * plano é ruído com aparência de superfície — e parece respeitável do mesmo
 * jeito.
 *
 * O DESVIO É MEDIDO POR RÉPLICA, não por lote. C_inter é função NÃO LINEAR de
 * rho (há |.| dentro do bloco, e entropia em C_rel), então dividir um ensemble
 * em lotes e olhar a dispersão dos lotes não dá a dispersão do estimador. Cada
 * réplica é um ensemble inteiro e independente, com outra semente-base; a
 * dispersão ENTRE réplicas é a barra de erro de uma célula.
 *
 * Fora do núcleo de propósito: isto é análise, usa stdio, lê /proc e chama o
 * autossolver denso de tests/.
 */
#include "dae.h"
#include "tests/jacobi.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double agora(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Pico de RSS do processo, MEDIDO. Estimar memória é como estimar tolerância:
   o número sai plausível e ninguém confere. */
static long pico_rss_kb(void)
{
  FILE *f = fopen("/proc/self/status", "r");
  char linha[256];
  long v = -1;
  if (!f) return -1;
  while (fgets(linha, sizeof(linha), f))
    if (strncmp(linha, "VmHWM:", 6) == 0) { sscanf(linha + 6, "%ld", &v); break; }
  fclose(f);
  return v;
}

/* Coerência de bloco na amostra `s`, com |.| DENTRO do bloco:
     C_MN = 2 sum_{i in M, j in N} |rho_ij|   (M < N)
     C_MM = sum_{i,j in M, i!=j} |rho_ij|
   e a identidade sum_M C_MM + sum_{M<N} C_MN = sum_{i!=j} |rho_ij| = C_l1. */
static void coerencias(const dae_rho_acc *A, const dae_graph *G, int32_t s,
                       double *c_inter, double *c_intra, double *c_l1)
{
  const size_t off = (size_t)s * (size_t)A->n * (size_t)A->n;
  double inter = 0.0, intra = 0.0;
  int32_t i, j;
  for (i = 0; i < A->n; ++i) {
    for (j = 0; j < A->n; ++j) {
      double m;
      if (i == j) continue;
      m = sqrt(A->re[off + (size_t)i * (size_t)A->n + (size_t)j] *
               A->re[off + (size_t)i * (size_t)A->n + (size_t)j] +
               A->im[off + (size_t)i * (size_t)A->n + (size_t)j] *
               A->im[off + (size_t)i * (size_t)A->n + (size_t)j]);
      if (G->module_of[i] == G->module_of[j]) intra += m; else inter += m;
    }
  }
  /* A soma sobre i != j percorre cada par ORDENADO, entao ela ja traz o fator 2
     de C_MN e a soma completa de C_MM. Nada a corrigir aqui — e o teste de
     n_modules = 1 confere isso: com um bloco so, inter tem de dar zero. */
  *c_inter = inter; *c_intra = intra; *c_l1 = inter + intra;
}

/* C_inter com a correcao de vies de Wardle-Kronberg aplicada POR ENTRADA:
   |rho_ij| <- sqrt( (n|rho^_ij|^2 - m2_ij)/(n-1) ), cortado em zero. */
static double c_inter_wk(const dae_rho_acc *A, const dae_graph *G, int32_t s,
                         int32_t n_traj)
{
  double inter = 0.0;
  int32_t i, j;
  for (i = 0; i < A->n; ++i)
    for (j = 0; j < A->n; ++j) {
      if (i == j || G->module_of[i] == G->module_of[j]) continue;
      inter += sqrt(dae_rho_acc_mod2_sem_vies(A, s, i, j, n_traj));
    }
  return inter;
}

/* Tridiagonalizacao de Householder, SEM acumular autovetores. E a diferenca
   entre 2 s e 135 s: a entropia precisa so dos autovalores, e Jacobi paga
   autovetores que ninguem usa. Medido no proprio problema (2N = 1040), Jacobi
   levou 135 s e isto leva ~2 s. */
static void tridiagonalizar(double *a, int32_t m, double *al, double *be)
{
  double *u = (double *)malloc((size_t)m * sizeof(double));
  double *pv = (double *)malloc((size_t)m * sizeof(double));
  int32_t k, i, j;
  for (k = 0; k < m - 2; ++k) {
    double normx = 0.0, alpha, r2, kappa;
    for (i = k + 1; i < m; ++i) normx += a[(size_t)i * (size_t)m + (size_t)k] *
                                        a[(size_t)i * (size_t)m + (size_t)k];
    normx = sqrt(normx);
    if (normx < 1e-300) continue;
    alpha = (a[(size_t)(k + 1) * (size_t)m + (size_t)k] > 0.0) ? -normx : normx;
    for (i = 0; i <= k; ++i) u[i] = 0.0;
    for (i = k + 1; i < m; ++i) u[i] = a[(size_t)i * (size_t)m + (size_t)k];
    u[k + 1] -= alpha;
    r2 = 0.0;
    for (i = k + 1; i < m; ++i) r2 += u[i] * u[i];
    if (r2 < 1e-300) continue;
    /* A <- A - kappa (u v^T + v u^T),  v = A u,  kappa = 2/|u|^2  */
    kappa = 2.0 / r2;
    for (i = 0; i < m; ++i) {
      double t = 0.0;
      for (j = k + 1; j < m; ++j) t += a[(size_t)i * (size_t)m + (size_t)j] * u[j];
      pv[i] = kappa * t;
    }
    { double uv = 0.0;
      for (i = k + 1; i < m; ++i) uv += u[i] * pv[i];
      for (i = 0; i < m; ++i) pv[i] -= 0.5 * kappa * uv * u[i]; }
    for (i = 0; i < m; ++i)
      for (j = 0; j < m; ++j)
        a[(size_t)i * (size_t)m + (size_t)j] -= u[i] * pv[j] + pv[i] * u[j];
  }
  for (i = 0; i < m; ++i) al[i] = a[(size_t)i * (size_t)m + (size_t)i];
  for (i = 0; i < m - 1; ++i) be[i] = a[(size_t)i * (size_t)m + (size_t)(i + 1)];
  free(u); free(pv);
}

/* S(rho) por autovalores, via o mergulho real 2N x 2N [[Re,-Im],[Im,Re]], cujos
   autovalores sao os de rho cada um DUPLICADO. Householder + bisseccao de Sturm
   (dae_tridiag_*), que sao as mesmas rotinas que dae_cheb e dae_metrics usam. */
static double entropia(const double *re, const double *im, int32_t n)
{
  const int32_t m = 2 * n;
  double *a = (double *)calloc((size_t)m * (size_t)m, sizeof(double));
  double *al = (double *)calloc((size_t)m, sizeof(double));
  double *be = (double *)calloc((size_t)m, sizeof(double));
  double s = 0.0;
  int32_t i, j;
  if (!a || !al || !be) { free(a); free(al); free(be); return NAN; }
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j) {
      const double r = re[(size_t)i * (size_t)n + (size_t)j];
      const double q = im[(size_t)i * (size_t)n + (size_t)j];
      a[(size_t)i * (size_t)m + (size_t)j] = r;
      a[(size_t)(i + n) * (size_t)m + (size_t)(j + n)] = r;
      a[(size_t)i * (size_t)m + (size_t)(j + n)] = -q;
      a[(size_t)(i + n) * (size_t)m + (size_t)j] = q;
    }
  tridiagonalizar(a, m, al, be);
  /* rho e PSD de traco 1, entao os autovalores vivem em [0, 1]; a margem cobre
     o erro de arredondamento da tridiagonalizacao. */
  for (i = 0; i < m; ++i) {
    const double lam = dae_tridiag_bisect(al, be, m, i, -1e-6, 1.0 + 1e-6);
    if (lam > 1e-14) s -= lam * log(lam);
  }
  free(a); free(al); free(be);
  return 0.5 * s;                 /* cada autovalor apareceu duas vezes */
}

/* C_rel = S(Delta[rho]) - S(rho), a medida canonica da teoria de recursos de
   coerencia de bloco. Delta zera tudo fora dos blocos diagonais. */
static double c_rel(const dae_rho_acc *A, const dae_graph *G, int32_t s)
{
  const size_t off = (size_t)s * (size_t)A->n * (size_t)A->n;
  const size_t tot = (size_t)A->n * (size_t)A->n;
  double *dr = (double *)malloc(tot * sizeof(double));
  double *di = (double *)malloc(tot * sizeof(double));
  double sr, sd;
  int32_t i, j;
  if (!dr || !di) { free(dr); free(di); return NAN; }
  memcpy(dr, A->re + off, tot * sizeof(double));
  memcpy(di, A->im + off, tot * sizeof(double));
  for (i = 0; i < A->n; ++i)
    for (j = 0; j < A->n; ++j)
      if (G->module_of[i] != G->module_of[j]) {
        dr[(size_t)i * (size_t)A->n + (size_t)j] = 0.0;
        di[(size_t)i * (size_t)A->n + (size_t)j] = 0.0;
      }
  sd = entropia(dr, di, A->n);
  sr = entropia(A->re + off, A->im + off, A->n);
  free(dr); free(di);
  return sd - sr;
}

/* UMA CELULA DA VARREDURA, em BLOCOS.
 *
 * O valor reportado sai do ensemble INTEIRO (menor vies, que cai como 1/sqrt(n));
 * a barra de erro e o vies saem da dispersao entre os R blocos. Isso da tres
 * coisas que nao se recuperam depois sem re-rodar, e por isso vao no CSV:
 *
 *   valor      C_inter corrigido por Wardle-Kronberg, no n cheio
 *   sd         desvio DA MEDIA, = sd entre blocos / sqrt(R)
 *   vies       residual estimado. Com C(n) = C_inf + B/sqrt(n), a media dos
 *              blocos (cada um de n/R) fica em C_inf + B*sqrt(R)/sqrt(n), entao
 *              vies_no_n_cheio = (media_blocos - cheio)/(sqrt(R) - 1). Sem
 *              ajuste nenhum, so a algebra da lei ja medida.
 *
 * E AS DUAS SUPERFICIES SAEM DO MESMO CONJUNTO. C_inter e eficiencia de
 * transferencia sao calculadas nos MESMOS blocos, entao a covariancia entre
 * elas e mensuravel — e a distancia entre as duas cristas, que e o observavel
 * de interesse, tem barra CORRELACIONADA. Parte do ruido cancela na diferenca, e a
 * distancia fica medida com precisao melhor que a de qualquer uma isolada. Se
 * viessem de rodadas separadas isso seria impossivel de recuperar.
 *
 * Dois acumuladores em memoria, nao R: um do bloco corrente e um do total. */
static int rodar_celula(dae_spec *S, int32_t n_total, int32_t blocos)
{
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  dae_metrics M;
  dae_metrics_cfg mc;
  dae_traj_cfg cfg;
  dae_rho_acc cheio, bloco;
  dae_status st;
  double *ci_b, *ef_b;
  double ci_cheio, ef_cheio, mci = 0.0, mef = 0.0, sci = 0.0, sef = 0.0, cov = 0.0;
  double t0 = agora();
  int32_t na, b, i, por_bloco = n_total / blocos;

  { dae_gen_params gp = S->gen; gp.seed = S->seed; st = dae_graph_build(&G, &gp); }
  if (st != DAE_OK) { fprintf(stderr, "grafo: %s\n", dae_strerror(st)); return 1; }
  memset(&mc, 0, sizeof(mc));
  if (dae_metrics_compute(&G, &mc, &M) != DAE_OK) { fprintf(stderr, "metricas\n"); return 1; }
  st = dae_hamiltonian(&H, &G.A, S->ham, S->gamma, S->norm, S->lanczos_steps, NULL);
  if (st != DAE_OK) { fprintf(stderr, "H: %s\n", dae_strerror(st)); return 1; }
  st = dae_cheb_init(&W, &H, S->lanczos_steps);
  if (st != DAE_OK) { fprintf(stderr, "cheb: %s\n", dae_strerror(st)); return 1; }

  /* A config tem ponteiros (init_re, gamma_tab): preencher campo a
     campo deixaria os demais com lixo da pilha. */
  dae_traj_cfg_default(&cfg);
  cfg.gamma_deph = S->gamma_deph; cfg.rho_stride = S->rho_stride;
  cfg.nt = S->nt; cfg.dt = S->t1 / (double)S->nt;
  cfg.sitio_inicial = S->init_site >= 0 ? S->init_site : 0;
  cfg.saida = DAE_SAIDA_ACUMULAR_RHO; cfg.n_traj = por_bloco;
  na = dae_traj_amostras(&cfg);

  if (dae_rho_acc_init_wk(&cheio, G.n, na) != DAE_OK) return 1;
  if (dae_rho_acc_init_wk(&bloco, G.n, na) != DAE_OK) return 1;
  ci_b = (double *)calloc((size_t)blocos, sizeof(double));
  ef_b = (double *)calloc((size_t)blocos, sizeof(double));

  for (b = 0; b < blocos; ++b) {
    size_t k, tot = (size_t)na * (size_t)G.n * (size_t)G.n;
    int32_t s;
    /* Semente-base propria por bloco, derivada do indice: fluxo compartilhado
       faria o resultado depender da ordem de execucao. */
    bloco.proxima = 0; bloco.somadas = 0;
    memset(bloco.re, 0, tot * sizeof(double));
    memset(bloco.im, 0, tot * sizeof(double));
    memset(bloco.m2, 0, tot * sizeof(double));
    st = dae_traj_ensemble(&cfg, S->seed + 1000003ULL * (uint64_t)(b + 1),
                           &W, &H, &bloco, NULL, NULL);
    if (st != DAE_OK) { fprintf(stderr, "bloco %d: %s\n", b, dae_strerror(st)); return 1; }
    /* soma no total ANTES de normalizar o bloco */
    for (k = 0; k < tot; ++k) {
      cheio.re[k] += bloco.re[k]; cheio.im[k] += bloco.im[k]; cheio.m2[k] += bloco.m2[k];
    }
    cheio.somadas += bloco.somadas;
    dae_rho_acc_finalizar(&bloco);
    ci_b[b] = c_inter_wk(&bloco, &G, na - 1, por_bloco);
    ef_b[b] = 0.0;
    if (S->target >= 0 && S->target < G.n) {
      for (s = 0; s < na; ++s)
        ef_b[b] += bloco.re[(size_t)s * (size_t)G.n * (size_t)G.n +
                            (size_t)S->target * (size_t)G.n + (size_t)S->target];
      ef_b[b] /= (double)na;
    }
  }
  dae_rho_acc_finalizar(&cheio);
  ci_cheio = c_inter_wk(&cheio, &G, na - 1, n_total);
  ef_cheio = 0.0;
  if (S->target >= 0 && S->target < G.n) {
    int32_t s;
    for (s = 0; s < na; ++s)
      ef_cheio += cheio.re[(size_t)s * (size_t)G.n * (size_t)G.n +
                           (size_t)S->target * (size_t)G.n + (size_t)S->target];
    ef_cheio /= (double)na;
  }

  for (i = 0; i < blocos; ++i) { mci += ci_b[i]; mef += ef_b[i]; }
  mci /= blocos; mef /= blocos;
  for (i = 0; i < blocos; ++i) {
    sci += (ci_b[i] - mci) * (ci_b[i] - mci);
    sef += (ef_b[i] - mef) * (ef_b[i] - mef);
    cov += (ci_b[i] - mci) * (ef_b[i] - mef);
  }
  sci = sqrt(sci / (blocos - 1)); sef = sqrt(sef / (blocos - 1));
  cov /= (blocos - 1);

  printf("%.6g,%.6g,%.17g,%.17g,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.1f\n",
         S->gen.ws_p, S->gamma_deph, M.lambda2, M.modularity_Q,
         G.n, n_total, blocos,
         ci_cheio, sci / sqrt((double)blocos), (mci - ci_cheio) / (sqrt((double)blocos) - 1.0),
         ef_cheio, sef / sqrt((double)blocos), (mef - ef_cheio) / (sqrt((double)blocos) - 1.0),
         cov / (double)blocos,
         (sci > 0 && sef > 0) ? cov / (sci * sef) : 0.0,
         M.lambda2_converged, agora() - t0);
  fflush(stdout);

  free(ci_b); free(ef_b);
  dae_rho_acc_free(&cheio); dae_rho_acc_free(&bloco);
  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
  return 0;
}

/* SERIE TEMPORAL do emaranhamento da rede.
 *
 * No setor de uma excitacao a concurrence par a par vale C_ij = 2|rho_ij|
 * exatamente (estado X com |11> nunca populado), entao o emaranhamento TOTAL
 * da rede e sum_{i!=j} |rho_ij| — que e a mesma soma da coerencia l1. Aqui ela
 * sai decomposta em dentro e entre modulos, em cada instante amostrado.
 *
 * Com a correcao de Wardle-Kronberg, e ela e MAIS necessaria aqui que no plano:
 * em t = 0 o estado e uma delta e quase todo rho_ij verdadeiro vale zero, que e
 * exatamente onde o vies do |.| e maximo. Sem corrigir, a curva comecaria num
 * patamar alto que e puro ruido e desceria — o oposto da fisica. */
static int rodar_serie(dae_spec *S, int32_t n_total, int32_t blocos)
{
  dae_graph G; dae_csr H; dae_cheb W;
  dae_traj_cfg cfg; dae_rho_acc cheio, bloco;
  dae_status st;
  double *tot_b, *int_b, *dt_b;
  int32_t na, b, s2, i, j, por_bloco = n_total / blocos;

  { dae_gen_params gp = S->gen; gp.seed = S->seed; st = dae_graph_build(&G, &gp); }
  if (st != DAE_OK) { fprintf(stderr, "grafo: %s\n", dae_strerror(st)); return 1; }
  st = dae_hamiltonian(&H, &G.A, S->ham, S->gamma, S->norm, S->lanczos_steps, NULL);
  if (st != DAE_OK) return 1;
  st = dae_cheb_init(&W, &H, S->lanczos_steps);
  if (st != DAE_OK) return 1;

  /* A config tem ponteiros (init_re, gamma_tab): preencher campo a
     campo deixaria os demais com lixo da pilha. */
  dae_traj_cfg_default(&cfg);
  cfg.gamma_deph = S->gamma_deph; cfg.rho_stride = S->rho_stride;
  cfg.nt = S->nt; cfg.dt = S->t1 / (double)S->nt;
  cfg.sitio_inicial = S->init_site >= 0 ? S->init_site : 0;
  cfg.saida = DAE_SAIDA_ACUMULAR_RHO; cfg.n_traj = por_bloco;
  na = dae_traj_amostras(&cfg);

  if (dae_rho_acc_init_wk(&cheio, G.n, na) != DAE_OK) return 1;
  if (dae_rho_acc_init_wk(&bloco, G.n, na) != DAE_OK) return 1;
  tot_b = (double *)calloc((size_t)blocos * (size_t)na, sizeof(double));
  int_b = (double *)calloc((size_t)blocos * (size_t)na, sizeof(double));
  dt_b  = (double *)calloc((size_t)na, sizeof(double));

  for (b = 0; b < blocos; ++b) {
    size_t k, tot = (size_t)na * (size_t)G.n * (size_t)G.n;
    bloco.proxima = 0; bloco.somadas = 0;
    memset(bloco.re, 0, tot * sizeof(double));
    memset(bloco.im, 0, tot * sizeof(double));
    memset(bloco.m2, 0, tot * sizeof(double));
    st = dae_traj_ensemble(&cfg, S->seed + 1000003ULL * (uint64_t)(b + 1),
                           &W, &H, &bloco, NULL, NULL);
    if (st != DAE_OK) { fprintf(stderr, "bloco %d: %s\n", b, dae_strerror(st)); return 1; }
    for (k = 0; k < tot; ++k) {
      cheio.re[k] += bloco.re[k]; cheio.im[k] += bloco.im[k]; cheio.m2[k] += bloco.m2[k];
    }
    cheio.somadas += bloco.somadas;
    dae_rho_acc_finalizar(&bloco);
    for (s2 = 0; s2 < na; ++s2) {
      double t = 0.0, e = 0.0;
      for (i = 0; i < G.n; ++i)
        for (j = 0; j < G.n; ++j) {
          double m;
          if (i == j) continue;
          m = sqrt(dae_rho_acc_mod2_sem_vies(&bloco, s2, i, j, por_bloco));
          t += m;
          if (G.module_of[i] != G.module_of[j]) e += m;
        }
      tot_b[(size_t)b * (size_t)na + (size_t)s2] = t;
      int_b[(size_t)b * (size_t)na + (size_t)s2] = e;
    }
  }
  dae_rho_acc_finalizar(&cheio);

  printf("passo,gamma_t,c_total,sd_total,c_inter,sd_inter,c_intra\n");
  for (s2 = 0; s2 < na; ++s2) {
    double t = 0.0, e = 0.0, mt = 0.0, me = 0.0, st_ = 0.0, se = 0.0;
    for (i = 0; i < G.n; ++i)
      for (j = 0; j < G.n; ++j) {
        double m;
        if (i == j) continue;
        m = sqrt(dae_rho_acc_mod2_sem_vies(&cheio, s2, i, j, n_total));
        t += m;
        if (G.module_of[i] != G.module_of[j]) e += m;
      }
    for (b = 0; b < blocos; ++b) {
      mt += tot_b[(size_t)b * (size_t)na + (size_t)s2];
      me += int_b[(size_t)b * (size_t)na + (size_t)s2];
    }
    mt /= blocos; me /= blocos;
    for (b = 0; b < blocos; ++b) {
      const double dt2 = tot_b[(size_t)b * (size_t)na + (size_t)s2] - mt;
      const double de2 = int_b[(size_t)b * (size_t)na + (size_t)s2] - me;
      st_ += dt2 * dt2; se += de2 * de2;
    }
    st_ = sqrt(st_ / (blocos - 1)) / sqrt((double)blocos);
    se  = sqrt(se  / (blocos - 1)) / sqrt((double)blocos);
    printf("%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
           s2 * S->rho_stride,
           S->gamma * (double)(s2 * S->rho_stride) * cfg.dt,
           t, st_, e, se, t - e);
  }
  free(tot_b); free(int_b); free(dt_b);
  dae_rho_acc_free(&cheio); dae_rho_acc_free(&bloco);
  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
  return 0;
}

int main(int argc, char **argv)
{
  dae_spec S;
  dae_error err;
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  dae_traj_cfg cfg;
  dae_status st;
  FILE *fp;
  char *texto = NULL;
  long tam;
  int32_t niveis[16], nn = 0, replicas = 8, i, k, r, na;
  int quer_crel = 0, grafo_varia = 0, celula = 0, n_total = 63000, blocos = 16;

  if (argc < 2) {
    fprintf(stderr,
      "uso: sonda <spec.json> [--niveis 25,50,100,200,400] [--replicas 8] [--crel]\n");
    return 2;
  }
  niveis[nn++] = 25; niveis[nn++] = 50; niveis[nn++] = 100;
  niveis[nn++] = 200; niveis[nn++] = 400;
  for (i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "--replicas") == 0 && i + 1 < argc) replicas = atoi(argv[++i]);
    else if (strcmp(argv[i], "--crel") == 0) quer_crel = 1;
    /* Cada replica com OUTRO grafo. Sem isto a sonda mede so o ruido das
       trajetorias com o grafo fixo, e a varredura real tem tambem a variacao
       entre realizacoes da religacao — que nao esta naquele numero e
       provavelmente domina em p alto. */
    else if (strcmp(argv[i], "--grafo-varia") == 0) grafo_varia = 1;
    else if (strcmp(argv[i], "--celula") == 0) celula = 1;
    else if (strcmp(argv[i], "--serie") == 0) celula = 2;
    else if (strcmp(argv[i], "--n") == 0 && i + 1 < argc) n_total = atoi(argv[++i]);
    else if (strcmp(argv[i], "--blocos") == 0 && i + 1 < argc) blocos = atoi(argv[++i]);
    else if (strcmp(argv[i], "--niveis") == 0 && i + 1 < argc) {
      char *p = argv[++i]; nn = 0;
      while (*p && nn < 16) {
        niveis[nn++] = atoi(p);
        while (*p && *p != ',') ++p;
        if (*p == ',') ++p;
      }
    } else { fprintf(stderr, "argumento desconhecido: %s\n", argv[i]); return 2; }
  }

  fp = fopen(argv[1], "rb");
  if (!fp) { fprintf(stderr, "nao consegui ler %s\n", argv[1]); return 1; }
  fseek(fp, 0, SEEK_END); tam = ftell(fp); fseek(fp, 0, SEEK_SET);
  texto = (char *)malloc((size_t)tam + 1);
  if (fread(texto, 1, (size_t)tam, fp) != (size_t)tam) return 1;
  texto[tam] = '\0'; fclose(fp);

  st = dae_spec_parse(&S, texto, &err);
  free(texto);
  if (st != DAE_OK) {
    fprintf(stderr, "%s:%d:%d: %s\n", argv[1], (int)err.line, (int)err.col, err.msg);
    return 1;
  }
  if (celula == 2) { int r2 = rodar_serie(&S, n_total, blocos); dae_spec_free(&S); return r2; }
  if (celula) { int r2 = rodar_celula(&S, n_total, blocos); dae_spec_free(&S); return r2; }

  { dae_gen_params gp = S.gen; gp.seed = S.seed; st = dae_graph_build(&G, &gp); }
  if (st != DAE_OK) { fprintf(stderr, "grafo: %s\n", dae_strerror(st)); return 1; }
  st = dae_hamiltonian(&H, &G.A, S.ham, S.gamma, S.norm, S.lanczos_steps, NULL);
  if (st != DAE_OK) { fprintf(stderr, "H: %s\n", dae_strerror(st)); return 1; }
  st = dae_cheb_init(&W, &H, S.lanczos_steps);
  if (st != DAE_OK) { fprintf(stderr, "cheb: %s\n", dae_strerror(st)); return 1; }

  /* A config tem ponteiros (init_re, gamma_tab): preencher campo a
     campo deixaria os demais com lixo da pilha. */
  dae_traj_cfg_default(&cfg);
  cfg.gamma_deph = S.gamma_deph; cfg.rho_stride = S.rho_stride;
  cfg.nt = S.nt; cfg.dt = S.t1 / (double)S.nt;
  cfg.sitio_inicial = S.init_site >= 0 ? S.init_site : 0;
  cfg.saida = DAE_SAIDA_ACUMULAR_RHO;
  cfg.n_traj = 1;
  na = dae_traj_amostras(&cfg);

  printf("# sonda  n=%d  nmod=%d  |E|=%d  gamma_deph=%g  amostras=%d  alvo=%d\n",
         G.n, G.nmod, H.nnz / 2, S.gamma_deph, na, S.target);
  /* pmod0 e LINEAR em rho, limitado e bem condicionado. E o discriminador: se
     ele convergir como 1/sqrt(n) e C_inter nao, a anomalia e da NAO LINEARIDADE
     do |.| e nao do PRNG. Se nem ele convergir, o problema e correlacao entre
     trajetorias, e ai nada mais importa ate isso ser resolvido. */
  printf("n_traj,replica,c_inter,c_intra,c_l1,ef_transfer,c_rel,pmod0,c_inter_deb,c_inter_wk,segundos\n");

  for (k = 0; k < nn; ++k) {
    for (r = 0; r < replicas; ++r) {
      dae_rho_acc A, Aq;
      double ci, ca, cl, ciq, caq, clq, ci_wk, ef = 0.0, cr = NAN, t0;
      int32_t s;
      if (grafo_varia) {
        dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
        { dae_gen_params gp = S.gen;
          gp.seed = S.seed + 7919ULL * (uint64_t)(r + 1);
          st = dae_graph_build(&G, &gp); }
        if (st != DAE_OK) { fprintf(stderr, "grafo r=%d: %s\n", r, dae_strerror(st)); return 1; }
        st = dae_hamiltonian(&H, &G.A, S.ham, S.gamma, S.norm, S.lanczos_steps, NULL);
        if (st != DAE_OK) { fprintf(stderr, "H r=%d: %s\n", r, dae_strerror(st)); return 1; }
        st = dae_cheb_init(&W, &H, S.lanczos_steps);
        if (st != DAE_OK) { fprintf(stderr, "cheb r=%d: %s\n", r, dae_strerror(st)); return 1; }
      }
      cfg.n_traj = niveis[k];
      st = dae_rho_acc_init_wk(&A, G.n, na);
      if (st != DAE_OK) { fprintf(stderr, "acumulador: %s\n", dae_strerror(st)); return 1; }
      /* Segundo acumulador sobre o PRIMEIRO QUARTO das mesmas trajetorias.
         E o que torna a correcao de vies gratuita: nenhuma trajetoria a mais,
         so mais um acumulador. Ver o comentario de `dedois` abaixo. */
      st = dae_rho_acc_init(&Aq, G.n, na);
      if (st != DAE_OK) { fprintf(stderr, "acumulador q: %s\n", dae_strerror(st)); return 1; }
      Aq.n_traj_alvo = niveis[k] / 4;
      t0 = agora();
      /* Semente-base DIFERENTE por replica, derivada do indice: replicas que
         compartilhassem fluxo teriam correlacao e a dispersao sairia menor que
         a verdadeira — barra de erro otimista e invisivel. */
      st = dae_traj_ensemble_dupla(&cfg, S.seed + 1000003ULL * (uint64_t)(r + 1),
                                   &W, &H, &A, &Aq, niveis[k] / 4);
      if (st != DAE_OK) { fprintf(stderr, "ensemble: %s\n", dae_strerror(st)); return 1; }
      dae_rho_acc_finalizar(&A);
      dae_rho_acc_finalizar(&Aq);
      coerencias(&A, &G, na - 1, &ci, &ca, &cl);
      coerencias(&Aq, &G, na - 1, &ciq, &caq, &clq);
      ci_wk = c_inter_wk(&A, &G, na - 1, niveis[k]);
      /* Eficiencia de transferencia: media TEMPORAL de rho_alvo,alvo sobre as
         amostras. Sai do MESMO ensemble que C_inter — se viesse de outra
         rodada, a distancia entre as duas cristas mediria ruido
         estatistico independente em vez de fisica. */
      if (S.target >= 0 && S.target < G.n) {
        for (s = 0; s < na; ++s)
          ef += A.re[(size_t)s * (size_t)G.n * (size_t)G.n +
                     (size_t)S.target * (size_t)G.n + (size_t)S.target];
        ef /= (double)na;
      }
      if (quer_crel) cr = c_rel(&A, &G, na - 1);
      { double p0 = 0.0;
        for (s = 0; s < G.n; ++s)
          if (G.module_of[s] == 0)
            p0 += A.re[(size_t)(na - 1) * (size_t)G.n * (size_t)G.n +
                       (size_t)s * (size_t)G.n + (size_t)s];
        /* CORRECAO DE VIES POR RICHARDSON EM 1/sqrt(n).
           C_inter e Sigma|rho_ij| e |.| e convexo, entao por Jensen o estimador
           e VIESADO PARA CIMA: E[Sigma|rho^_ij|] >= Sigma|E rho^_ij|. Nas
           entradas onde rho_ij ~ 0 — a esmagadora maioria das N^2 — o vies vale
           ~sigma_ij sqrt(pi/2)/sqrt(n), e somado sobre N^2 entradas ele domina.
           Medido: C(n) = C_inf + B/sqrt(n).
           Com C(n) e C(n/4) do MESMO ensemble, 2C(n) - C(n/4) cancela B
           exatamente, e as duas estimativas sao correlacionadas, o que mantem a
           variancia menor do que se viessem de ensembles separados. */
        printf("%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.4f\n",
               niveis[k], r, ci, ca, cl, ef, cr, p0, 2.0 * ci - ciq, ci_wk,
               agora() - t0); }
      fflush(stdout);
      dae_rho_acc_free(&A); dae_rho_acc_free(&Aq);
    }
  }
  printf("# pico_rss_kb %ld\n", pico_rss_kb());

  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G); dae_spec_free(&S);
  return 0;
}
