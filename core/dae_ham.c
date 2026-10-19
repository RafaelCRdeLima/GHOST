#include "dae_ham.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

dae_status dae_hamiltonian(dae_csr *H, const dae_csr *A, dae_ham_kind kind,
                           double gamma, dae_norm_kind norm,
                           int32_t lanczos_steps, double *scale_out)
{
  const int32_t n = A ? A->n : 0;
  int32_t i, p, out, extra;
  double scale = 1.0;
  dae_status st;

  if (scale_out) *scale_out = 1.0;
  if (!H || !A || n <= 0) return DAE_ERR_PARAM;
  if (!(gamma == gamma)) return DAE_ERR_PARAM;

  /* Quantas linhas ainda não têm entrada diagonal (a laplaciana precisa dela). */
  extra = 0;
  if (kind == DAE_H_LAPLACIAN) {
    for (i = 0; i < n; ++i) {
      int has = 0;
      for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p)
        if (A->colind[p] == i) { has = 1; break; }
      if (!has) ++extra;
    }
  }

  st = dae_csr_alloc(H, n, A->nnz + extra);
  if (st != DAE_OK) return st;

  out = 0;
  for (i = 0; i < n; ++i) {
    H->rowptr[i] = out;
    if (kind == DAE_H_ADJACENCY) {
      for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p) {
        H->colind[out] = A->colind[p];
        H->val[out]    = -A->val[p];
        ++out;
      }
    } else {
      double deg = 0.0, self = 0.0;
      int placed = 0;
      for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p) {
        if (A->colind[p] == i) self += A->val[p];
        else                   deg  += A->val[p];
      }
      /* mantém a linha ordenada inserindo a diagonal na posição certa */
      for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p) {
        const int32_t c = A->colind[p];
        if (!placed && c >= i) {
          H->colind[out] = i;
          H->val[out]    = deg - self;
          ++out; placed = 1;
          if (c == i) continue;               /* a diagonal original foi absorvida */
        }
        H->colind[out] = c;
        H->val[out]    = -A->val[p];
        ++out;
      }
      if (!placed) {
        H->colind[out] = i;
        H->val[out]    = deg - self;
        ++out;
      }
    }
  }
  H->rowptr[n] = out;
  H->nnz = out;

  switch (norm) {
    case DAE_NORM_NONE:
      scale = 1.0;
      break;
    case DAE_NORM_MEAN_DEGREE: {
      double s = 0.0;
      for (i = 0; i < n; ++i)
        for (p = A->rowptr[i]; p < A->rowptr[i + 1]; ++p)
          if (A->colind[p] != i) s += fabs(A->val[p]);
      scale = s / (double)n;                    /* grau médio ponderado */
      break;
    }
    case DAE_NORM_SPECTRAL: {
      double lo, hi, r;
      const int32_t m = (lanczos_steps > 0) ? lanczos_steps : 40;
      st = dae_csr_lanczos_bounds(H, m, &lo, &hi);
      if (st != DAE_OK) { dae_csr_free(H); return st; }
      r = fabs(lo) > fabs(hi) ? fabs(lo) : fabs(hi);
      scale = r;
      break;
    }
    default:
      dae_csr_free(H);
      return DAE_ERR_PARAM;
  }
  if (!(scale > 0.0) || !(scale == scale)) scale = 1.0;   /* grafo sem arestas */

  {
    const double f = gamma / scale;
    for (p = 0; p < H->nnz; ++p) H->val[p] *= f;
  }
  if (scale_out) *scale_out = scale;
  return DAE_OK;
}
