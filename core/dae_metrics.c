#include "dae_metrics.h"
#include "dae_ham.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Teto de memória da base de Lanczos, em doubles (64 MB). A reortogonalização
   completa exige guardar a base inteira, e é ela que sustenta a precisão no
   regime modular; onde não couber, o teto de passos cai e a bandeira de
   convergência avisa. */
#define DAE_LANCZOS_BUDGET 8000000

void dae_metrics_cfg_default(dae_metrics_cfg *c)
{
  if (!c) return;
  c->lambda2_max_steps = 0;
  c->lambda2_tol = 0.0;
  c->bfs_sources = 0;
}

/* --------------------------------------------------------------- BFS ---- */

static void dae_bfs(const dae_csr *A, int32_t src, int32_t *dist, int32_t *q)
{
  int32_t head = 0, tail = 0, i;
  for (i = 0; i < A->n; ++i) dist[i] = -1;
  dist[src] = 0;
  q[tail++] = src;
  while (head < tail) {
    const int32_t u = q[head++];
    int32_t p;
    for (p = A->rowptr[u]; p < A->rowptr[u + 1]; ++p) {
      const int32_t v = A->colind[p];
      if (dist[v] < 0) { dist[v] = dist[u] + 1; q[tail++] = v; }
    }
  }
}

/* Vetor de Ritz do `alvo`-esimo menor autovalor da tridiagonal, projetado de
 * volta para o espaco cheio e deflacionado do vetor constante.
 *
 * Existe como funcao porque o layout espectral precisa do MESMO calculo para os
 * dois menores pares — e duplicar isso duplicaria tambem os pontos de deflacao,
 * que sao a parte que nao pode erodir. Devolve ||x|| antes de normalizar: se
 * cancelou tudo, o "vetor" nao e vetor.
 */
static double dae_ritz(const double *al, const double *be, int32_t m, int32_t alvo,
                       double glo, double ghi, const double *V, int32_t n,
                       double *y, double *dsc, double *x,
                       double *theta_out, double *vazamento_out)
{
  double theta, xn = 0.0;
  int32_t i;
  int it;
  theta = dae_tridiag_bisect(al, be, m, alvo, glo - 1.0, ghi + 1.0);
  for (i = 0; i < m; ++i) y[i] = 1.0 / sqrt((double)m);
  for (it = 0; it < 3; ++it) {
    double yn = 0.0;
    dae_tridiag_solve(al, be, m, theta - 1e-11 * (1.0 + fabs(theta)) - 1e-300, y, dsc);
    for (i = 0; i < m; ++i) yn += y[i] * y[i];
    yn = sqrt(yn);
    if (!(yn > 0.0)) break;
    for (i = 0; i < m; ++i) y[i] /= yn;
  }
  for (i = 0; i < n; ++i) x[i] = 0.0;
  for (i = 0; i < m; ++i) {
    const double *vi = V + (size_t)i * (size_t)n;
    const double c = y[i];
    int32_t k;
    for (k = 0; k < n; ++k) x[k] += c * vi[k];
  }
  { double c = 0.0, xnp = 0.0;            /* o vetor de Ritz VIVE em 1-perp */
    for (i = 0; i < n; ++i) c += x[i];
    for (i = 0; i < n; ++i) xnp += x[i] * x[i];
    xnp = sqrt(xnp);
    if (vazamento_out) {
      const double leak = (xnp > 0.0) ? fabs(c) / (sqrt((double)n) * xnp) : 0.0;
      if (leak > *vazamento_out) *vazamento_out = leak;
    }
    c /= (double)n;
    for (i = 0; i < n; ++i) x[i] -= c; }
  for (i = 0; i < n; ++i) xn += x[i] * x[i];
  xn = sqrt(xn);
  if (xn > 0.0) for (i = 0; i < n; ++i) x[i] /= xn;
  if (theta_out) *theta_out = theta;
  return xn;
}

/* ------------------------------------------------------------ lambda_2 --
 *
 * DEFLACAO DO VETOR CONSTANTE, em quatro pontos. L tem 1 como autovetor de
 * autovalor 0, e lambda_2 e o PROXIMO autovalor: qualquer vazamento da direcao
 * constante sequestra o valor de Ritz minimo e devolve ~0. Medido no P_64 com
 * deflacao incompleta: a componente constante do vetor de Ritz cresce de 1e-17
 * para 3e-14 ao longo de 60 passos e, em m = n, o vetor de Ritz VIRA o vetor
 * constante (media -1/sqrt(64)), com quociente de Rayleigh 5e-31 e residuo
 * 7e-16 — ou seja, "convergido" e completamente errado.
 *
 *   (A) o teto de passos e n-1, e nao n: o operador vive num espaco de
 *       dimensao n-1, e obrigar o Lanczos a produzir a n-esima direcao faz ele
 *       inventar ruido — que contem exatamente o vetor constante;
 *   (B) deflaciona antes de cada passada de reortogonalizacao;
 *   (C) deflaciona depois delas, porque reprojetar contra os v_q reintroduz a
 *       componente que eles carregam;
 *   (D) deflaciona o proprio vetor de Ritz antes do quociente de Rayleigh.
 *
 * (A), (B) e (C) sao MUTUAMENTE REDUNDANTES quanto ao lambda_2 final: remover
 * qualquer um sozinho nao muda o resultado, e so a remocao dos tres reproduz o
 * defeito. Comentario nao roda na CI, entao a defesa de verdade e
 * `lambda2_const_leak`, medido aqui e verificado em t94_metrics.c: ele ataca o
 * MECANISMO (a componente constante) em vez do sintoma (o lambda_2 errado), e
 * por isso morde quando qualquer um dos pontos some, muito antes de o
 * resultado final estragar.
 *
 * (D) tem outro papel e nao aparece nesse indicador: com (A)-(C) fora, ele
 * converte uma resposta silenciosamente errada num valor visivelmente alto com
 * a bandeira de convergencia em zero. Nenhum indicador o pega.
 */

static dae_status dae_fiedler(const dae_csr *L, const dae_metrics_cfg *cfg,
                              dae_metrics *out, int32_t nvec, double *vecs)
{
  const int32_t n = L->n;
  int32_t mmax, m, j, i, k;
  double *V = NULL, *w = NULL, *al = NULL, *be = NULL;
  double *y = NULL, *dsc = NULL, *x = NULL, *r = NULL;
  double glo, ghi, normL, tol, beta = 0.0;
  dae_status st = DAE_OK;

  dae_csr_gershgorin(L, &glo, &ghi);
  normL = (fabs(ghi) > fabs(glo)) ? fabs(ghi) : fabs(glo);
  if (!(normL > 0.0)) normL = 1.0;
  tol = (cfg->lambda2_tol > 0.0) ? cfg->lambda2_tol : 1e-8;

  mmax = (cfg->lambda2_max_steps > 0) ? cfg->lambda2_max_steps : 300;
  /* O operador vive no COMPLEMENTO do vetor constante, de dimensão n-1. Deixar
     m chegar a n obriga o Lanczos a inventar uma direção que não existe, e o
     que ele inventa é ruído contendo justamente o vetor constante — que tem
     autovalor 0 e sequestra o valor de Ritz mínimo. Medido: em P_64 com m = 64
     o quociente de Rayleigh cai para 5e-31 e o vetor de Ritz vira o próprio
     vetor constante (média -1/sqrt(64)). */
  if (mmax > n - 1) mmax = n - 1;
  if ((double)mmax * (double)n > (double)DAE_LANCZOS_BUDGET)
    mmax = (int32_t)((double)DAE_LANCZOS_BUDGET / (double)n);
  if (mmax < 2) mmax = 2;

  V   = (double *)malloc((size_t)mmax * (size_t)n * sizeof(double));
  w   = (double *)malloc((size_t)n * sizeof(double));
  x   = (double *)malloc((size_t)n * sizeof(double));
  r   = (double *)malloc((size_t)n * sizeof(double));
  al  = (double *)malloc((size_t)mmax * sizeof(double));
  be  = (double *)malloc((size_t)mmax * sizeof(double));
  y   = (double *)malloc((size_t)mmax * sizeof(double));
  dsc = (double *)malloc((size_t)mmax * sizeof(double));
  if (!V || !w || !x || !r || !al || !be || !y || !dsc) { st = DAE_ERR_ALLOC; goto done; }

  /* v_0 determinístico, deflacionado do vetor constante (que é o autovetor de
     autovalor 0 da laplaciana) e normalizado. */
  dae_csr_lanczos_start(V, n);
  {
    double s = 0.0, nn = 0.0;
    for (i = 0; i < n; ++i) s += V[i];
    s /= (double)n;
    for (i = 0; i < n; ++i) { V[i] -= s; nn += V[i] * V[i]; }
    nn = sqrt(nn);
    if (!(nn > 0.0)) { st = DAE_ERR_PARAM; goto done; }
    for (i = 0; i < n; ++i) V[i] /= nn;
  }
  out->lambda2_const_leak = 0.0;

  out->lambda2 = 0.0;
  out->lambda2_residual = normL;
  out->lambda2_converged = 0;
  out->lambda2_steps = 0;

  for (j = 0; j < mmax; ++j) {
    double a = 0.0, nn = 0.0, s = 0.0;
    const double *vj = V + (size_t)j * (size_t)n;

    dae_csr_spmv_real(L, vj, w);
    for (i = 0; i < n; ++i) s += w[i];               /* deflaciona 1 a cada passo */
    s /= (double)n;
    for (i = 0; i < n; ++i) w[i] -= s;

    for (i = 0; i < n; ++i) a += w[i] * vj[i];
    al[j] = a;
    for (i = 0; i < n; ++i) w[i] -= a * vj[i];
    if (j > 0) {
      const double *vp = V + (size_t)(j - 1) * (size_t)n;
      for (i = 0; i < n; ++i) w[i] -= beta * vp[i];
    }
    /* Reortogonalização COMPLETA, duas passadas, com o vetor constante tratado
       como mais um membro do conjunto — porque os v_q carregam componente
       constante residual e reprojetar contra eles a reintroduz. Deflacionar só
       uma vez, logo após o produto matriz-vetor, deixa essa componente crescer
       passo a passo (medido: de 1e-17 para 3e-14 em 60 passos) até dominar. */
    for (k = 0; k < 2; ++k) {
      int32_t q;
      double c = 0.0;
      for (i = 0; i < n; ++i) c += w[i];
      c /= (double)n;
      for (i = 0; i < n; ++i) w[i] -= c;
      for (q = 0; q <= j; ++q) {
        const double *vq = V + (size_t)q * (size_t)n;
        double d = 0.0;
        for (i = 0; i < n; ++i) d += w[i] * vq[i];
        for (i = 0; i < n; ++i) w[i] -= d * vq[i];
      }
    }
    { double c = 0.0;
      for (i = 0; i < n; ++i) c += w[i];
      c /= (double)n;
      for (i = 0; i < n; ++i) w[i] -= c; }
    for (i = 0; i < n; ++i) nn += w[i] * w[i];
    nn = sqrt(nn);
    beta = nn;
    if (j + 1 < mmax) be[j] = beta;
    m = j + 1;
    out->lambda2_steps = m;

    /* Avaliação de convergência: valor de Ritz mínimo por Sturm, autovetor da
       tridiagonal por iteração inversa, resíduo EXPLÍCITO em L. */
    if ((m % 5) == 0 || m == 1 || nn < 1e-13 || j + 1 == mmax) {
      double theta = 0.0, rq = 0.0, res = 0.0;
      const double xn = dae_ritz(al, be, m, 0, glo, ghi, V, n, y, dsc, x,
                                 &theta, &out->lambda2_const_leak);
      if (xn > 1e-8) {                        /* cancelou tudo: não é vetor */
        dae_csr_spmv_real(L, x, r);
        for (i = 0; i < n; ++i) rq += x[i] * r[i];      /* quociente de Rayleigh */
        for (i = 0; i < n; ++i) { const double d = r[i] - rq * x[i]; res += d * d; }
        res = sqrt(res);
        out->lambda2 = rq;
        out->lambda2_residual = res;
        out->lambda2_converged = (res < tol * normL) ? 1 : 0;
        if (out->lambda2_converged) break;
      }
    }
    if (nn < 1e-13) break;                             /* subespaço invariante */
    if (j + 1 < mmax) {
      double *vn = V + (size_t)(j + 1) * (size_t)n;
      double leak = 0.0;
      for (i = 0; i < n; ++i) vn[i] = w[i] / nn;
      for (i = 0; i < n; ++i) leak += vn[i];
      leak = fabs(leak) / sqrt((double)n);      /* |<v, 1/sqrt(n)>| */
      if (leak > out->lambda2_const_leak) out->lambda2_const_leak = leak;
    }
  }

  /* Os `nvec` menores pares de Ritz, extraidos da MESMA base — e por isso com
     as mesmas quatro deflacoes. O layout espectral usa os dois primeiros como
     coordenadas: no espaco ja deflacionado do vetor constante, eles sao os
     autovetores 2 e 3 da laplaciana. */
  if (nvec > 0 && vecs) {
    int32_t k;
    for (k = 0; k < nvec; ++k) {
      double theta = 0.0;
      if (k >= m) { for (i = 0; i < n; ++i) vecs[(size_t)k * (size_t)n + (size_t)i] = 0.0; continue; }
      dae_ritz(al, be, m, k, glo, ghi, V, n, y, dsc, x, &theta, &out->lambda2_const_leak);
      for (i = 0; i < n; ++i) vecs[(size_t)k * (size_t)n + (size_t)i] = x[i];
    }
  }

done:
  free(V); free(w); free(x); free(r); free(al); free(be); free(y); free(dsc);
  return st;
}

/* -------------------------------------------------------------- layout --
 *
 * Layout espectral: autovetores 2 e 3 da laplaciana como coordenadas (x, y).
 *
 * Nao e escolha estetica. Para grafo modular o embedding espectral SEPARA os
 * modulos no plano, que e exatamente o que se quer enxergar; e reaproveita o
 * maquinario que ja existe — e literalmente o vetor de Fiedler do lambda_2,
 * mais o proximo.
 *
 * Alternativa descartada: ordenar os sitios por indice e dobrar em linhas. Isso
 * desenha a ORDEM DE ROTULAGEM, nao a rede, e para SBM ou lista importada a
 * ordem nao significa nada.
 *
 * As coordenadas saem CRUAS. Normalizar para a caixa do canvas e trabalho de
 * quem desenha, que e quem sabe o tamanho da caixa.
 */
dae_status dae_layout_espectral(const dae_graph *G, float *xy)
{
  dae_csr L;
  dae_metrics M;
  dae_metrics_cfg cfg;
  double *vecs;
  int32_t i, n;
  dae_status st;

  if (!G || !xy || G->n <= 1) return DAE_ERR_PARAM;
  n = G->n;
  st = dae_hamiltonian(&L, &G->A, DAE_H_LAPLACIAN, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) return st;

  vecs = (double *)calloc(2u * (size_t)n, sizeof(double));
  if (!vecs) { dae_csr_free(&L); return DAE_ERR_ALLOC; }

  dae_metrics_cfg_default(&cfg);
  memset(&M, 0, sizeof(M));
  st = dae_fiedler(&L, &cfg, &M, 2, vecs);
  if (st == DAE_OK) {
    for (i = 0; i < n; ++i) {
      xy[2 * i]     = (float)vecs[i];
      xy[2 * i + 1] = (float)vecs[(size_t)n + (size_t)i];
    }
  }
  free(vecs);
  dae_csr_free(&L);
  return st;
}

/* ------------------------------------------------------------ principal - */

dae_status dae_metrics_compute(const dae_graph *G, const dae_metrics_cfg *cfg,
                               dae_metrics *out)
{
  dae_metrics_cfg def;
  const dae_csr *A;
  dae_csr L;
  int32_t n, i, p, nself = 0, noff = 0;
  int32_t *dist = NULL, *queue = NULL;
  double *krow = NULL, *WM = NULL, *KM = NULL;
  double Wtot = 0.0;
  dae_status st = DAE_OK;

  if (!G || !out) return DAE_ERR_PARAM;
  if (!cfg) { dae_metrics_cfg_default(&def); cfg = &def; }
  memset(out, 0, sizeof(*out));
  A = &G->A;
  n = A->n;
  if (n <= 0) return DAE_ERR_PARAM;

  for (i = 0; i < n; ++i)
    for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p) {
      if (A->colind[p] == i) ++nself;
      else                   ++noff;
    }
  out->n_edges = noff / 2 + nself;
  out->mean_degree = (double)noff / (double)n;

  /* --- modularidade de Newman, ponderada, sobre a partição do grafo --- */
  krow = (double *)calloc((size_t)n, sizeof(double));
  WM = (double *)calloc((size_t)(G->nmod > 0 ? G->nmod : 1), sizeof(double));
  KM = (double *)calloc((size_t)(G->nmod > 0 ? G->nmod : 1), sizeof(double));
  if (!krow || !WM || !KM) { st = DAE_ERR_ALLOC; goto done; }
  for (i = 0; i < n; ++i) {
    for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p) krow[i] += A->val[p];
    Wtot += krow[i];
  }
  Wtot *= 0.5;
  for (i = 0; i < n; ++i) {
    const int32_t mi = G->module_of ? G->module_of[i] : 0;
    KM[mi] += krow[i];
    for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p) {
      const int32_t mj = G->module_of ? G->module_of[A->colind[p]] : 0;
      if (mj == mi) WM[mi] += A->val[p];
    }
  }
  if (Wtot > 0.0) {
    double Q = 0.0;
    for (i = 0; i < G->nmod; ++i) {
      const double eM = 0.5 * WM[i] / Wtot;
      const double aM = KM[i] / (2.0 * Wtot);
      Q += eM - aM * aM;
    }
    out->modularity_Q = Q;
  }

  /* --- componentes e comprimento médio de caminho (BFS, sem pesos) --- */
  dist  = (int32_t *)malloc((size_t)n * sizeof(int32_t));
  queue = (int32_t *)malloc((size_t)n * sizeof(int32_t));
  if (!dist || !queue) { st = DAE_ERR_ALLOC; goto done; }
  {
    int32_t *seen = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    int32_t comp = 0;
    if (!seen) { st = DAE_ERR_ALLOC; goto done; }
    for (i = 0; i < n; ++i) {
      if (seen[i]) continue;
      ++comp;
      dae_bfs(A, i, dist, queue);
      { int32_t t; for (t = 0; t < n; ++t) if (dist[t] >= 0) seen[t] = 1; }
    }
    out->n_components = comp;
    free(seen);
  }
  {
    const int32_t want = (cfg->bfs_sources > 0 && cfg->bfs_sources < n)
                       ? cfg->bfs_sources : n;
    const int32_t nsrc = (n <= 5000 || cfg->bfs_sources > 0) ? want
                       : (256 < n ? 256 : n);
    double soma = 0.0;
    int64_t pares = 0;
    int32_t s;
    out->path_len_exact = (nsrc == n) ? 1 : 0;
    for (s = 0; s < nsrc; ++s) {
      /* fontes espaçadas de forma determinística: nada de PRNG aqui */
      const int32_t src = (int32_t)(((int64_t)s * (int64_t)n) / (int64_t)nsrc);
      int32_t t;
      dae_bfs(A, src, dist, queue);
      for (t = 0; t < n; ++t)
        if (t != src && dist[t] > 0) { soma += (double)dist[t]; ++pares; }
    }
    out->mean_path_len = (pares > 0) ? soma / (double)pares : 0.0;
  }

  /* --- lambda_2 --- */
  st = dae_hamiltonian(&L, A, DAE_H_LAPLACIAN, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) goto done;
  if (out->n_components > 1) {
    /* Grafo desconexo tem lambda_2 = 0 EXATO, com multiplicidade igual ao
       número de componentes. Não há o que iterar, e mandar Lanczos procurar
       num espectro com autovalor 0 degenerado só produziria ruído. */
    out->lambda2 = 0.0;
    out->lambda2_residual = 0.0;
    out->lambda2_converged = 1;
    out->lambda2_steps = 0;
  } else {
    st = dae_fiedler(&L, cfg, out, 0, NULL);
  }
  dae_csr_free(&L);

done:
  free(krow); free(WM); free(KM); free(dist); free(queue);
  return st;
}
