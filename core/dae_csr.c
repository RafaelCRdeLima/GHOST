#include "dae_csr.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

dae_status dae_csr_alloc(dae_csr *A, int32_t n, int32_t nnz)
{
  size_t cap;
  if (!A || n <= 0 || nnz < 0) return DAE_ERR_PARAM;
  cap = (size_t)(nnz > 0 ? nnz : 1);
  A->n = n;
  A->nnz = nnz;
  A->rowptr = (int32_t *)malloc(((size_t)n + 1) * sizeof(int32_t));
  A->colind = (int32_t *)malloc(cap * sizeof(int32_t));
  A->val    = (double  *)malloc(cap * sizeof(double));
  if (!A->rowptr || !A->colind || !A->val) { dae_csr_free(A); return DAE_ERR_ALLOC; }
  memset(A->rowptr, 0, ((size_t)n + 1) * sizeof(int32_t));
  return DAE_OK;
}

void dae_csr_free(dae_csr *A)
{
  if (!A) return;
  free(A->rowptr); free(A->colind); free(A->val);
  A->rowptr = NULL; A->colind = NULL; A->val = NULL;
  A->n = 0; A->nnz = 0;
}

dae_status dae_csr_from_edges(dae_csr *A, int32_t n,
                              const int32_t *ei, const int32_t *ej,
                              const double *w, int32_t ne,
                              int32_t *n_dropped)
{
  int32_t i, e, p, cap, out, dropped, rowstart;
  int32_t *cursor;
  dae_status st;

  if (n_dropped) *n_dropped = 0;
  if (!A || n <= 0 || ne < 0 || (ne > 0 && (!ei || !ej))) return DAE_ERR_PARAM;

  cap = 0;
  for (e = 0; e < ne; ++e) {
    if (ei[e] < 0 || ei[e] >= n || ej[e] < 0 || ej[e] >= n) return DAE_ERR_PARAM;
    if (w && !(w[e] == w[e])) return DAE_ERR_PARAM;          /* NaN */
    cap += (ei[e] == ej[e]) ? 1 : 2;
  }

  st = dae_csr_alloc(A, n, cap);
  if (st != DAE_OK) return st;

  for (e = 0; e < ne; ++e) {
    A->rowptr[ei[e] + 1]++;
    if (ei[e] != ej[e]) A->rowptr[ej[e] + 1]++;
  }
  for (i = 0; i < n; ++i) A->rowptr[i + 1] += A->rowptr[i];

  cursor = (int32_t *)malloc((size_t)n * sizeof(int32_t));
  if (!cursor) { dae_csr_free(A); return DAE_ERR_ALLOC; }
  for (i = 0; i < n; ++i) cursor[i] = A->rowptr[i];

  for (e = 0; e < ne; ++e) {
    const double v = w ? w[e] : 1.0;
    p = cursor[ei[e]]++;
    A->colind[p] = ej[e]; A->val[p] = v;
    if (ei[e] != ej[e]) {
      p = cursor[ej[e]]++;
      A->colind[p] = ei[e]; A->val[p] = v;
    }
  }
  free(cursor);

  /* Ordena cada linha e compacta descartando duplicatas. A compactação escreve
     em `out` <= `p` sempre, então o in-place é seguro. */
  out = 0; dropped = 0; rowstart = 0;
  for (i = 0; i < n; ++i) {
    const int32_t beg = rowstart;
    const int32_t end = A->rowptr[i + 1];
    for (p = beg + 1; p < end; ++p) {                 /* inserção: linhas são curtas */
      const int32_t c = A->colind[p];
      const double  v = A->val[p];
      int32_t q = p - 1;
      while (q >= beg && A->colind[q] > c) {
        A->colind[q + 1] = A->colind[q];
        A->val[q + 1]    = A->val[q];
        --q;
      }
      A->colind[q + 1] = c;
      A->val[q + 1]    = v;
    }
    rowstart = end;
    A->rowptr[i] = out;
    for (p = beg; p < end; ++p) {
      if (p > beg && A->colind[p] == A->colind[p - 1]) { ++dropped; continue; }
      A->colind[out] = A->colind[p];
      A->val[out]    = A->val[p];
      ++out;
    }
  }
  A->rowptr[n] = out;
  A->nnz = out;
  if (n_dropped) *n_dropped = dropped;

  if (out > 0 && out < cap) {                         /* encolhe, sem falhar se não der */
    int32_t *ci = (int32_t *)realloc(A->colind, (size_t)out * sizeof(int32_t));
    double  *vv = (double  *)realloc(A->val,    (size_t)out * sizeof(double));
    if (ci) A->colind = ci;
    if (vv) A->val    = vv;
  }
  return DAE_OK;
}

void dae_csr_spmv(const dae_csr *A,
                  const double *xre, const double *xim,
                  double *yre, double *yim)
{
  int32_t i, p;
  for (i = 0; i < A->n; ++i) {
    const int32_t end = A->rowptr[i + 1];
    double sr = 0.0, si = 0.0;
    for (p = A->rowptr[i]; p < end; ++p) {
      const double  v = A->val[p];
      const int32_t c = A->colind[p];
      sr += v * xre[c];
      si += v * xim[c];
    }
    yre[i] = sr;
    yim[i] = si;
  }
}

void dae_csr_spmv_real(const dae_csr *A, const double *x, double *y)
{
  int32_t i, p;
  for (i = 0; i < A->n; ++i) {
    const int32_t end = A->rowptr[i + 1];
    double s = 0.0;
    for (p = A->rowptr[i]; p < end; ++p) s += A->val[p] * x[A->colind[p]];
    y[i] = s;
  }
}

void dae_csr_gershgorin(const dae_csr *A, double *lo, double *hi)
{
  int32_t i, p;
  double l = 0.0, h = 0.0;
  int first = 1;
  for (i = 0; i < A->n; ++i) {
    const int32_t end = A->rowptr[i + 1];
    double diag = 0.0, radius = 0.0;
    for (p = A->rowptr[i]; p < end; ++p) {
      if (A->colind[p] == i) diag += A->val[p];
      else                   radius += fabs(A->val[p]);
    }
    if (first) { l = diag - radius; h = diag + radius; first = 0; }
    else {
      if (diag - radius < l) l = diag - radius;
      if (diag + radius > h) h = diag + radius;
    }
  }
  *lo = l; *hi = h;
}

int32_t dae_tridiag_count(const double *al, const double *be,
                          int32_t m, double mu)
{
  int32_t i, cnt;
  double d = al[0] - mu;
  cnt = (d < 0.0) ? 1 : 0;
  for (i = 1; i < m; ++i) {
    if (d == 0.0) d = 1e-300;
    d = (al[i] - mu) - (be[i - 1] * be[i - 1]) / d;
    if (d < 0.0) ++cnt;
  }
  return cnt;
}

double dae_tridiag_bisect(const double *al, const double *be, int32_t m,
                          int32_t target, double lo, double hi)
{
  int it;
  for (it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) break;
    if (dae_tridiag_count(al, be, m, mid) <= target) lo = mid;
    else                                           hi = mid;
  }
  return 0.5 * (lo + hi);
}

void dae_csr_lanczos_start(double *v, int32_t n)
{
  /* Constante fixa, não dae_rng: ver a justificativa em dae_csr.h. */
  uint64_t z = 0xDAEDA105C0FFEE01ULL;
  int32_t i;
  for (i = 0; i < n; ++i) {
    uint64_t t;
    z += 0x9E3779B97F4A7C15ULL;
    t = z;
    t = (t ^ (t >> 30)) * 0xBF58476D1CE4E5B9ULL;
    t = (t ^ (t >> 27)) * 0x94D049BB133111EBULL;
    t ^= t >> 31;
    v[i] = (double)(t >> 11) / 9007199254740992.0 * 2.0 - 1.0;
  }
}

void dae_tridiag_solve(const double *al, const double *be, int32_t m,
                       double sigma, double *x, double *d)
{
  int32_t i;
  for (i = 0; i < m; ++i) d[i] = al[i] - sigma;
  for (i = 1; i < m; ++i) {
    double f;
    if (fabs(d[i - 1]) < 1e-300) d[i - 1] = 1e-300;   /* sigma quase autovalor */
    f = be[i - 1] / d[i - 1];
    d[i] -= f * be[i - 1];
    x[i] -= f * x[i - 1];
  }
  if (fabs(d[m - 1]) < 1e-300) d[m - 1] = 1e-300;
  x[m - 1] /= d[m - 1];
  for (i = m - 2; i >= 0; --i) x[i] = (x[i] - be[i] * x[i + 1]) / d[i];
}

dae_status dae_csr_lanczos_bounds(const dae_csr *A, int32_t m,
                                  double *lo, double *hi)
{
  const int32_t n = A->n;
  double *v = NULL, *vp = NULL, *wr = NULL, *al = NULL, *be = NULL;
  double glo, ghi, beta = 0.0, nrm;
  int32_t i, j, msteps;
  dae_status st = DAE_OK;

  if (!A || !lo || !hi || m <= 0) return DAE_ERR_PARAM;
  if (m > n) m = n;

  v    = (double *)malloc((size_t)n * sizeof(double));
  vp   = (double *)malloc((size_t)n * sizeof(double));
  wr   = (double *)malloc((size_t)n * sizeof(double));
  al   = (double *)malloc((size_t)m * sizeof(double));
  be   = (double *)malloc((size_t)m * sizeof(double));
  if (!v || !vp || !wr || !al || !be) { st = DAE_ERR_ALLOC; goto done; }

  /* Vetor inicial determinístico: mesmo limite espectral em qualquer alvo,
     sem tocar no fluxo compartilhado do PRNG. */
  dae_csr_lanczos_start(v, n);
  nrm = 0.0;
  for (i = 0; i < n; ++i) nrm += v[i] * v[i];
  nrm = sqrt(nrm);
  if (nrm == 0.0) { st = DAE_ERR_PARAM; goto done; }
  for (i = 0; i < n; ++i) v[i] /= nrm;
  memset(vp, 0, (size_t)n * sizeof(double));

  msteps = 0;
  for (j = 0; j < m; ++j) {
    double a = 0.0;
    dae_csr_spmv_real(A, v, wr);
    for (i = 0; i < n; ++i) a += wr[i] * v[i];
    for (i = 0; i < n; ++i) wr[i] -= a * v[i] + beta * vp[i];
    /* reortogonalização local contra v: barata e segura o suficiente aqui,
       porque a margem final vem de beta e do recorte por Gershgorin. */
    { double c = 0.0;
      for (i = 0; i < n; ++i) c += wr[i] * v[i];
      for (i = 0; i < n; ++i) wr[i] -= c * v[i]; }
    al[j] = a;
    msteps = j + 1;
    nrm = 0.0;
    for (i = 0; i < n; ++i) nrm += wr[i] * wr[i];
    nrm = sqrt(nrm);
    beta = nrm;
    if (j + 1 < m) be[j] = beta;
    if (nrm < 1e-13) break;                   /* subespaço invariante */
    for (i = 0; i < n; ++i) { vp[i] = v[i]; v[i] = wr[i] / nrm; }
  }

  dae_csr_gershgorin(A, &glo, &ghi);
  if (msteps == 1) { *lo = al[0] - beta; *hi = al[0] + beta; }
  else {
    const double tmin = dae_tridiag_bisect(al, be, msteps, 0,          glo, ghi);
    const double tmax = dae_tridiag_bisect(al, be, msteps, msteps - 1, glo, ghi);
    *lo = tmin - beta;
    *hi = tmax + beta;
  }

done:
  free(v); free(vp); free(wr); free(al); free(be);
  return st;
}
