#include "dae_graph.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- arestas */

typedef struct { int32_t *i, *j; double *w; int32_t ne, cap; } dae_edges;

static dae_status dae_edges_init(dae_edges *E, int32_t cap)
{
  if (cap < 8) cap = 8;
  E->ne = 0; E->cap = cap;
  E->i = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  E->j = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  E->w = (double  *)malloc((size_t)cap * sizeof(double));
  if (!E->i || !E->j || !E->w) return DAE_ERR_ALLOC;
  return DAE_OK;
}

static void dae_edges_free(dae_edges *E)
{
  free(E->i); free(E->j); free(E->w);
  E->i = NULL; E->j = NULL; E->w = NULL; E->ne = 0; E->cap = 0;
}

static dae_status dae_edges_push(dae_edges *E, int32_t a, int32_t b, double w)
{
  if (E->ne == E->cap) {
    const int32_t nc = E->cap * 2;
    int32_t *ni = (int32_t *)realloc(E->i, (size_t)nc * sizeof(int32_t));
    int32_t *nj = (int32_t *)realloc(E->j, (size_t)nc * sizeof(int32_t));
    double  *nw = (double  *)realloc(E->w, (size_t)nc * sizeof(double));
    if (ni) E->i = ni;
    if (nj) E->j = nj;
    if (nw) E->w = nw;
    if (!ni || !nj || !nw) return DAE_ERR_ALLOC;
    E->cap = nc;
  }
  E->i[E->ne] = a; E->j[E->ne] = b; E->w[E->ne] = w;
  ++E->ne;
  return DAE_OK;
}

/* ------------------------------------------- conjunto de arestas existentes
 *
 * Endereçamento aberto com lápide. Serve SÓ para consulta: sem ele, uma
 * religação que caísse sobre uma aresta já existente seria descartada pela
 * deduplicação da CSR e |E| encolheria em silêncio — matando justamente a
 * disciplina de "religar mantendo |E| fixo". Nunca se itera sobre esta tabela:
 * a ordem de varredura de uma tabela de dispersão não é portável entre alvos,
 * e a ordem dos sorteios faz parte do contrato de reprodutibilidade. */

#define DAE_ESET_EMPTY 0ULL
#define DAE_ESET_DEAD  1ULL

typedef struct { uint64_t *key; uint64_t mask; } dae_eset;

static uint64_t dae_ekey(int32_t a, int32_t b, int32_t n)
{
  const int32_t lo = (a < b) ? a : b;
  const int32_t hi = (a < b) ? b : a;
  return (uint64_t)lo * (uint64_t)n + (uint64_t)hi + 2ULL;
}

static uint64_t dae_emix(uint64_t k)
{
  k ^= k >> 33; k *= 0xFF51AFD7ED558CCDULL;
  k ^= k >> 33; k *= 0xC4CEB9FE1A85EC53ULL;
  k ^= k >> 33;
  return k;
}

static dae_status dae_eset_init(dae_eset *S, int32_t expected)
{
  uint64_t cap = 16;
  while (cap < (uint64_t)expected * 4ULL) cap <<= 1;
  S->key = (uint64_t *)calloc((size_t)cap, sizeof(uint64_t));
  if (!S->key) return DAE_ERR_ALLOC;
  S->mask = cap - 1ULL;
  return DAE_OK;
}

static void dae_eset_free(dae_eset *S) { free(S->key); S->key = NULL; S->mask = 0; }

static int dae_eset_has(const dae_eset *S, uint64_t k)
{
  uint64_t h = dae_emix(k) & S->mask;
  for (;;) {
    const uint64_t v = S->key[h];
    if (v == DAE_ESET_EMPTY) return 0;
    if (v == k) return 1;
    h = (h + 1ULL) & S->mask;
  }
}

static void dae_eset_add(dae_eset *S, uint64_t k)
{
  uint64_t h = dae_emix(k) & S->mask;
  for (;;) {
    const uint64_t v = S->key[h];
    if (v == DAE_ESET_EMPTY || v == DAE_ESET_DEAD || v == k) { S->key[h] = k; return; }
    h = (h + 1ULL) & S->mask;
  }
}

static void dae_eset_del(dae_eset *S, uint64_t k)
{
  uint64_t h = dae_emix(k) & S->mask;
  for (;;) {
    const uint64_t v = S->key[h];
    if (v == DAE_ESET_EMPTY) return;
    if (v == k) { S->key[h] = DAE_ESET_DEAD; return; }
    h = (h + 1ULL) & S->mask;
  }
}

/* ------------------------------------------------------------- parâmetros */

uint64_t dae_graph_fingerprint(const dae_graph *G)
{
  uint64_t h = 1469598103934665603ULL;      /* FNV-1a, 64 bits */
  int32_t i;
  if (!G) return 0ULL;
  for (i = 0; i <= G->A.n; ++i) {
    h ^= (uint64_t)(uint32_t)G->A.rowptr[i]; h *= 1099511628211ULL;
  }
  for (i = 0; i < G->A.nnz; ++i) {
    h ^= (uint64_t)(uint32_t)G->A.colind[i]; h *= 1099511628211ULL;
  }
  /* Os pesos entram pelos bits, nao pelo valor: j_par trocado por j_perp muda a
     digital mesmo que a estrutura seja identica. */
  for (i = 0; i < G->A.nnz; ++i) {
    union { double d; uint64_t u; } c;
    c.d = G->A.val[i];
    h ^= c.u; h *= 1099511628211ULL;
  }
  for (i = 0; i < G->n; ++i) {
    h ^= (uint64_t)(uint32_t)(G->module_of ? G->module_of[i] : 0);
    h *= 1099511628211ULL;
  }
  return h;
}

void dae_gen_params_default(dae_gen_params *p)
{
  if (!p) return;
  memset(p, 0, sizeof(*p));
  p->kind = DAE_G_PATH;
  p->seed = 12345ULL;
  p->n_par = 200;
  p->n_perp = 13;
  p->seam_shift = 3;
  p->longitudinal_closed = 0;
  p->j_par = 1.0;
  p->j_perp = 1.0;
  p->n_modules = 1;
  p->conn_mode = DAE_REWIRE;
  p->n = 64;
  p->dim = 6;
  p->rows = 8;
  p->cols = 8;
}

void dae_graph_free(dae_graph *G)
{
  if (!G) return;
  dae_csr_free(&G->A);
  free(G->module_of);
  free(G->xy);
  memset(G, 0, sizeof(*G));
}

static dae_status dae_graph_alloc_node_arrays(dae_graph *G, int32_t n)
{
  G->n = n;
  G->module_of = (int32_t *)calloc((size_t)n, sizeof(int32_t));
  G->xy = (float *)calloc(2u * (size_t)n, sizeof(float));
  if (!G->module_of || !G->xy) return DAE_ERR_ALLOC;
  G->nmod = 1;
  G->geom = DAE_GEOM_NENHUMA;
  return DAE_OK;
}

static void dae_layout_circle(dae_graph *G)
{
  int32_t i;
  for (i = 0; i < G->n; ++i) {
    const double th = 6.283185307179586 * (double)i / (double)G->n;
    G->xy[2 * i]     = (float)cos(th);
    G->xy[2 * i + 1] = (float)sin(th);
  }
}

dae_status dae_graph_from_edges(dae_graph *G, int32_t n,
                                const int32_t *ei, const int32_t *ej,
                                const double *w, int32_t ne)
{
  dae_status st;
  if (!G || n <= 0) return DAE_ERR_PARAM;
  memset(G, 0, sizeof(*G));
  st = dae_csr_from_edges(&G->A, n, ei, ej, w, ne, &G->n_dropped);
  if (st != DAE_OK) { dae_graph_free(G); return st; }
  st = dae_graph_alloc_node_arrays(G, n);
  if (st != DAE_OK) { dae_graph_free(G); return st; }
  dae_layout_circle(G);
  return DAE_OK;
}

/* -------------------------------------------------------------- religação */

/* Watts-Strogatz. No modo DAE_REWIRE o número de arestas é INVARIANTE: cada
 * aresta sorteada troca um extremo, e a troca só é aceita se não colidir com
 * uma aresta existente. Uma colisão aceita seria descartada depois pela
 * deduplicação da CSR e |E| cairia sem aviso — o que faria o transporte parecer
 * pior por menos arestas, não por topologia.
 *
 * No modo DAE_ADD a aresta nova é uma CÓPIA da sorteada com um extremo trocado,
 * herdando o peso: a única diferença entre os dois modos passa a ser |E|.
 * Ver CONVENTIONS.md, parte 3.1. */
static dae_status dae_rewire(dae_edges *E, int32_t n, double p,
                             dae_conn_mode mode, dae_rng *g, int32_t *n_failed)
{
  dae_eset S;
  dae_status st;
  const int32_t base_ne = E->ne;
  int32_t e;

  *n_failed = 0;
  if (p <= 0.0 || base_ne == 0) return DAE_OK;
  if (p > 1.0) return DAE_ERR_PARAM;

  st = dae_eset_init(&S, base_ne * 2);
  if (st != DAE_OK) return st;
  for (e = 0; e < base_ne; ++e) dae_eset_add(&S, dae_ekey(E->i[e], E->j[e], n));

  for (e = 0; e < base_ne; ++e) {          /* ordem crescente: contrato */
    int attempt;
    if (dae_rng_uniform(g) >= p) continue;
    for (attempt = 0; attempt < 100; ++attempt) {
      const int32_t b = dae_rng_below(g, n);
      uint64_t k;
      if (b == E->i[e]) continue;
      k = dae_ekey(E->i[e], b, n);
      if (dae_eset_has(&S, k)) continue;
      if (mode == DAE_REWIRE) {
        dae_eset_del(&S, dae_ekey(E->i[e], E->j[e], n));
        E->j[e] = b;
      } else {
        st = dae_edges_push(E, E->i[e], b, E->w[e]);
        if (st != DAE_OK) { dae_eset_free(&S); return st; }
      }
      dae_eset_add(&S, k);
      break;
    }
    if (attempt == 100) ++(*n_failed);
  }
  dae_eset_free(&S);
  return DAE_OK;
}

/* ------------------------------------------------------------- geradores */

static dae_status dae_gen_microtubule(dae_edges *E, dae_graph *G,
                                      const dae_gen_params *p)
{
  const int32_t np = p->n_par, nq = p->n_perp;
  const int32_t n = np * nq;
  const int32_t nmod = (p->n_modules > 0) ? p->n_modules : 1;
  int32_t m, q;
  dae_status st;

  st = dae_graph_alloc_node_arrays(G, n);
  if (st != DAE_OK) return st;
  G->nmod = nmod;
  G->n_par = np;
  G->n_perp = nq;
  G->geom = DAE_GEOM_LATTICE;

  for (m = 0; m < np; ++m)
    for (q = 0; q < nq; ++q) {
      const int32_t j = m * nq + q;
      G->module_of[j] = (int32_t)(((int64_t)m * (int64_t)nmod) / (int64_t)np);
      G->xy[2 * j]     = (float)m;          /* rede desenrolada; o cilindro 3D */
      G->xy[2 * j + 1] = (float)q;          /* a interface deriva de (m, q)    */
    }

  for (m = 0; m < np; ++m) {
    for (q = 0; q < nq; ++q) {
      const int32_t j = m * nq + q;
      if (m + 1 < np) {
        st = dae_edges_push(E, j, (m + 1) * nq + q, p->j_par);
        if (st != DAE_OK) return st;
      } else if (p->longitudinal_closed) {
        st = dae_edges_push(E, j, q, p->j_par);
        if (st != DAE_OK) return st;
      }
      if (q + 1 < nq) {
        st = dae_edges_push(E, j, m * nq + q + 1, p->j_perp);
        if (st != DAE_OK) return st;
      } else if (nq > 2) {
        /* COSTURA: (m, nq-1) liga com (m + seam_shift, 0).
           Com pontas abertas as arestas que caem fora do cilindro somem, e
           somem TODAS numa ponta só — é a assimetria plus/minus end, física,
           não defeito de montagem. Ver CONVENTIONS.md, parte 5. */
        const int32_t ms = m + p->seam_shift;
        if (p->longitudinal_closed) {
          int32_t mm = ms % np;
          if (mm < 0) mm += np;
          st = dae_edges_push(E, j, mm * nq, p->j_perp);
          if (st != DAE_OK) return st;
        } else if (ms >= 0 && ms < np) {
          st = dae_edges_push(E, j, ms * nq, p->j_perp);
          if (st != DAE_OK) return st;
        }
      }
    }
  }
  return DAE_OK;
}

static dae_status dae_gen_sbm(dae_edges *E, dae_graph *G,
                              const dae_gen_params *p, dae_rng *g)
{
  const int32_t n = p->n;
  const int32_t nmod = (p->n_modules > 0) ? p->n_modules : 1;
  int32_t i, j;
  dae_status st;

  if (n < 2 || n > 20000) return DAE_ERR_TOOBIG;
  if (p->p_in < 0.0 || p->p_in > 1.0 || p->p_out < 0.0 || p->p_out > 1.0)
    return DAE_ERR_PARAM;

  st = dae_graph_alloc_node_arrays(G, n);
  if (st != DAE_OK) return st;
  G->nmod = nmod;
  for (i = 0; i < n; ++i)
    G->module_of[i] = (int32_t)(((int64_t)i * (int64_t)nmod) / (int64_t)n);

  /* Varredura i < j em ordem crescente: a ordem dos sorteios é contrato. */
  for (i = 0; i < n; ++i)
    for (j = i + 1; j < n; ++j) {
      const double pr = (G->module_of[i] == G->module_of[j]) ? p->p_in : p->p_out;
      if (dae_rng_uniform(g) < pr) {
        st = dae_edges_push(E, i, j, 1.0);
        if (st != DAE_OK) return st;
      }
    }
  return DAE_OK;
}

dae_status dae_graph_build(dae_graph *G, const dae_gen_params *p)
{
  dae_edges E;
  dae_rng g;
  int32_t n = 0, i;
  dae_status st;

  if (!G || !p) return DAE_ERR_PARAM;
  memset(G, 0, sizeof(*G));
  memset(&E, 0, sizeof(E));
  dae_rng_seed(&g, p->seed);

  switch (p->kind) {
    case DAE_G_PATH:
    case DAE_G_CYCLE:
      if (p->n < 2) return DAE_ERR_PARAM;
      n = p->n;
      break;
    case DAE_G_COMPLETE:
      if (p->n < 2) return DAE_ERR_PARAM;
      if (p->n > 20000) return DAE_ERR_TOOBIG;
      n = p->n;
      break;
    case DAE_G_HYPERCUBE:
      if (p->dim < 1 || p->dim > 24) return DAE_ERR_PARAM;
      n = (int32_t)1 << p->dim;
      break;
    case DAE_G_GRID2D:
      if (p->rows < 1 || p->cols < 1 || p->rows * p->cols < 2) return DAE_ERR_PARAM;
      n = p->rows * p->cols;
      break;
    case DAE_G_MICROTUBULE:
      if (p->n_par < 2 || p->n_perp < 2) return DAE_ERR_PARAM;
      if ((int64_t)p->n_par * (int64_t)p->n_perp > 20000000LL) return DAE_ERR_TOOBIG;
      n = p->n_par * p->n_perp;
      break;
    case DAE_G_SBM:
      n = p->n;
      break;
    case DAE_G_EDGELIST:
      return DAE_ERR_PARAM;                 /* use dae_graph_from_edges */
    default:
      return DAE_ERR_PARAM;
  }

  st = dae_edges_init(&E, 4 * n);
  if (st != DAE_OK) { dae_edges_free(&E); return st; }

  switch (p->kind) {
    case DAE_G_MICROTUBULE:
      st = dae_gen_microtubule(&E, G, p);
      break;
    case DAE_G_SBM:
      st = dae_gen_sbm(&E, G, p, &g);
      n = G->n;
      break;
    default:
      st = dae_graph_alloc_node_arrays(G, n);
      if (st != DAE_OK) break;
      dae_layout_circle(G);
      switch (p->kind) {
        case DAE_G_PATH:
        case DAE_G_CYCLE:
          for (i = 0; i + 1 < n && st == DAE_OK; ++i) st = dae_edges_push(&E, i, i + 1, 1.0);
          if (st == DAE_OK && p->kind == DAE_G_CYCLE) st = dae_edges_push(&E, n - 1, 0, 1.0);
          for (i = 0; i < n; ++i) { G->xy[2 * i] = (float)i; G->xy[2 * i + 1] = 0.0f; }
          if (p->kind == DAE_G_CYCLE) dae_layout_circle(G);
          break;
        case DAE_G_COMPLETE: {
          int32_t j;
          for (i = 0; i < n && st == DAE_OK; ++i)
            for (j = i + 1; j < n && st == DAE_OK; ++j) st = dae_edges_push(&E, i, j, 1.0);
          break;
        }
        case DAE_G_HYPERCUBE: {
          int32_t b;
          for (i = 0; i < n && st == DAE_OK; ++i)
            for (b = 0; b < p->dim && st == DAE_OK; ++b) {
              const int32_t j = i ^ ((int32_t)1 << b);
              if (j > i) st = dae_edges_push(&E, i, j, 1.0);
            }
          break;
        }
        case DAE_G_GRID2D: {
          int32_t r, c;
          for (r = 0; r < p->rows && st == DAE_OK; ++r)
            for (c = 0; c < p->cols && st == DAE_OK; ++c) {
              const int32_t j = r * p->cols + c;
              G->xy[2 * j] = (float)c; G->xy[2 * j + 1] = (float)r;
              if (c + 1 < p->cols) st = dae_edges_push(&E, j, j + 1, 1.0);
              if (st == DAE_OK && r + 1 < p->rows) st = dae_edges_push(&E, j, j + p->cols, 1.0);
            }
          break;
        }
        default: break;
      }
      break;
  }
  if (st != DAE_OK) { dae_edges_free(&E); dae_graph_free(G); return st; }

  st = dae_rewire(&E, n, p->ws_p, p->conn_mode, &g, &G->n_rewire_failed);
  if (st != DAE_OK) { dae_edges_free(&E); dae_graph_free(G); return st; }

  st = dae_csr_from_edges(&G->A, n, E.i, E.j, E.w, E.ne, &G->n_dropped);
  dae_edges_free(&E);
  if (st != DAE_OK) { dae_graph_free(G); return st; }
  return DAE_OK;
}
