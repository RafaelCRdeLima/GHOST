#include "dae_cheb.h"
#include "dae_bessel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static void dae_cheb_release(dae_cheb *W)
{
  free(W->t0re); free(W->t0im); free(W->t1re); free(W->t1im);
  free(W->t2re); free(W->t2im); free(W->acre); free(W->acim);
  free(W->ckre); free(W->ckim); free(W->jbuf);
  memset(W, 0, sizeof(*W));
}

void dae_cheb_free(dae_cheb *W) { if (W) dae_cheb_release(W); }

dae_status dae_cheb_init(dae_cheb *W, const dae_csr *H, int32_t lanczos_steps)
{
  double glo, ghi;
  size_t nb;
  if (!W || !H || H->n <= 0) return DAE_ERR_PARAM;
  memset(W, 0, sizeof(*W));
  W->n = H->n;

  dae_csr_gershgorin(H, &glo, &ghi);
  W->lo = glo; W->hi = ghi;

  if (lanczos_steps > 0) {
    double llo, lhi;
    if (dae_csr_lanczos_bounds(H, lanczos_steps, &llo, &lhi) == DAE_OK && lhi > llo) {
      const double mid  = 0.5 * (llo + lhi);
      const double half = 0.5 * (lhi - llo) * 1.05 + 1e-12 * (fabs(mid) + 1.0);
      double clo = mid - half, chi = mid + half;
      if (clo < glo) clo = glo;          /* nunca mais largo que a cota rigorosa */
      if (chi > ghi) chi = ghi;
      if (chi > clo) { W->lo = clo; W->hi = chi; W->lanczos_used = 1; }
    }
  }

  W->a = 0.5 * (W->hi - W->lo);
  W->b = 0.5 * (W->hi + W->lo);
  if (!(W->a > 0.0)) W->a = 1.0;         /* H múltiplo da identidade: H~ = 0 */

  nb = (size_t)W->n * sizeof(double);
  W->t0re = (double *)malloc(nb); W->t0im = (double *)malloc(nb);
  W->t1re = (double *)malloc(nb); W->t1im = (double *)malloc(nb);
  W->t2re = (double *)malloc(nb); W->t2im = (double *)malloc(nb);
  W->acre = (double *)malloc(nb); W->acim = (double *)malloc(nb);
  if (!W->t0re || !W->t0im || !W->t1re || !W->t1im ||
      !W->t2re || !W->t2im || !W->acre || !W->acim) {
    dae_cheb_release(W);
    return DAE_ERR_ALLOC;
  }
  W->coef_valid = 0;
  W->k_used = -1;
  return DAE_OK;
}

int32_t dae_cheb_order(const dae_cheb *W, double dt)
{
  double k = 1.2 * W->a * fabs(dt) + 20.0;
  if (k > 2.0e6) k = 2.0e6;
  return (k < 4.0) ? 4 : (int32_t)k + 1;
}

static dae_status dae_cheb_ensure(dae_cheb *W, int32_t k)
{
  double *a, *b, *c;
  if (k + 1 <= W->kcap) return DAE_OK;
  a = (double *)realloc(W->ckre, (size_t)(k + 1) * sizeof(double));
  if (a) W->ckre = a; else return DAE_ERR_ALLOC;
  b = (double *)realloc(W->ckim, (size_t)(k + 1) * sizeof(double));
  if (b) W->ckim = b; else return DAE_ERR_ALLOC;
  c = (double *)realloc(W->jbuf, (size_t)(k + 1) * sizeof(double));
  if (c) W->jbuf = c; else return DAE_ERR_ALLOC;
  W->kcap = k + 1;
  return DAE_OK;
}

static dae_status dae_cheb_coeffs(dae_cheb *W, double dt)
{
  const double alpha = W->a * dt;
  int32_t K = dae_cheb_order(W, dt), k;
  int guard;
  dae_status st;

  /* QUEM MANDA É A CAUDA, não a fórmula.
     K ~ 1.2 alpha + 20 é só o chute inicial. A margem que ela deixa acima do
     ponto de retorno, 0.2 alpha + 20, é assintoticamente suficiente mas NÃO na
     faixa alpha ~ 20 a 100: a transição de Airy de J_k(alpha) em k = alpha tem
     largura ~alpha^(1/3), e ali a fórmula corta a série com o coeficiente ainda
     grande. Medido: |2 J_45(20)| = 1.8e-12 e |2 J_81(50)| = 1.9e-11, o que
     tiraria a unitariedade do propagador justamente nessa faixa — 4.2e-13 de
     deriva da norma POR PASSO em alpha = 50, contra 5.5e-17 em alpha = 5.
     Em alpha = 300 a fórmula volta a bastar, e é por isso que o defeito é
     não-monotônico e passa despercebido num teste de um alpha só. */
  for (guard = 0; ; ++guard) {
    if (K > 1000000 || guard > 60) return DAE_ERR_TOOBIG;
    st = dae_cheb_ensure(W, K);
    if (st != DAE_OK) return st;
    st = dae_bessel_j_array(fabs(alpha), K, W->jbuf);
    if (st != DAE_OK) return st;
    if (fabs(2.0 * W->jbuf[K]) < DAE_CHEB_TAIL) break;
    K += (K / 4) + 8;
  }
  if (alpha < 0.0) for (k = 1; k <= K; k += 2) W->jbuf[k] = -W->jbuf[k];

  /* Agora apara o excesso: o laço acima garante que a cauda está abaixo do
     limiar, este remove os termos que sobraram por cima. */
  while (K > 2 && fabs(2.0 * W->jbuf[K]) < DAE_CHEB_TAIL) --K;

  for (k = 0; k <= K; ++k) {
    const double m = (k == 0) ? W->jbuf[0] : 2.0 * W->jbuf[k];
    switch (k & 3) {                                     /* (-i)^k */
      case 0:  W->ckre[k] =  m;  W->ckim[k] =  0.0; break;
      case 1:  W->ckre[k] = 0.0; W->ckim[k] = -m;   break;
      case 2:  W->ckre[k] = -m;  W->ckim[k] =  0.0; break;
      default: W->ckre[k] = 0.0; W->ckim[k] =  m;   break;
    }
  }
  W->k_used = K;
  W->cached_dt = dt;
  W->coef_valid = 1;
  return DAE_OK;
}

/* w <- H~ v = (H v - b v) / a */
static void dae_cheb_apply(const dae_cheb *W, const dae_csr *H,
                           const double *vre, const double *vim,
                           double *wre, double *wim)
{
  const int32_t n = H->n;
  const double b = W->b, inva = 1.0 / W->a;
  int32_t i;
  dae_csr_spmv(H, vre, vim, wre, wim);
  for (i = 0; i < n; ++i) {
    wre[i] = (wre[i] - b * vre[i]) * inva;
    wim[i] = (wim[i] - b * vim[i]) * inva;
  }
}

/* acumulador += c_k * T_k */
static void dae_cheb_accum(dae_cheb *W, int32_t k,
                           const double *tre, const double *tim, int32_t n)
{
  const double cr = W->ckre[k], ci = W->ckim[k];
  int32_t i;
  for (i = 0; i < n; ++i) {
    W->acre[i] += cr * tre[i] - ci * tim[i];
    W->acim[i] += cr * tim[i] + ci * tre[i];
  }
}

dae_status dae_cheb_step(dae_cheb *W, const dae_csr *H, double dt,
                         double *psire, double *psiim, dae_cheb_info *info)
{
  int32_t i, k, n;
  double n2in = 0.0, n2out = 0.0, cb, sb;
  double *p0re, *p0im, *p1re, *p1im, *p2re, *p2im, *tre, *tim;
  dae_status st;

  if (!W || !H || !psire || !psiim || H->n != W->n) return DAE_ERR_PARAM;
  n = W->n;

  if (!W->coef_valid || dt != W->cached_dt) {
    st = dae_cheb_coeffs(W, dt);
    if (st != DAE_OK) return st;
  }

  for (i = 0; i < n; ++i) n2in += psire[i] * psire[i] + psiim[i] * psiim[i];

  p0re = W->t0re; p0im = W->t0im;
  p1re = W->t1re; p1im = W->t1im;
  p2re = W->t2re; p2im = W->t2im;

  memcpy(p0re, psire, (size_t)n * sizeof(double));       /* T_0 = psi */
  memcpy(p0im, psiim, (size_t)n * sizeof(double));

  {
    const double cr = W->ckre[0], ci = W->ckim[0];
    for (i = 0; i < n; ++i) {
      W->acre[i] = cr * p0re[i] - ci * p0im[i];
      W->acim[i] = cr * p0im[i] + ci * p0re[i];
    }
  }

  if (W->k_used >= 1) {
    dae_cheb_apply(W, H, p0re, p0im, p1re, p1im);        /* T_1 = H~ psi */
    dae_cheb_accum(W, 1, p1re, p1im, n);
  }
  for (k = 2; k <= W->k_used; ++k) {
    dae_cheb_apply(W, H, p1re, p1im, p2re, p2im);
    for (i = 0; i < n; ++i) {                            /* T_{k+1} = 2 H~ T_k - T_{k-1} */
      p2re[i] = 2.0 * p2re[i] - p0re[i];
      p2im[i] = 2.0 * p2im[i] - p0im[i];
    }
    dae_cheb_accum(W, k, p2re, p2im, n);
    tre  = p0re; tim  = p0im;                            /* rotação dos três buffers */
    p0re = p1re; p0im = p1im;
    p1re = p2re; p1im = p2im;
    p2re = tre;  p2im = tim;
  }

  /* Fase global exp(-i b dt). Numa grade uniforme, cos e sin sao as MESMAS
     constantes em todo passo, e o defeito pitagorico delas — cb^2+sb^2 = 1 a
     menos de meia ulp — vira um fator multiplicativo FIXO sobre a norma. Isso
     produz deriva linear e sistematica: 5.5e-17 por passo, medidos, constantes
     ao longo de quatro decadas de numero de passos. Renormalizar o par custa
     uma raiz por passo (nao por elemento) e remove a deriva na origem. */
  cb =  cos(W->b * dt);
  sb = -sin(W->b * dt);
  {
    const double r = sqrt(cb * cb + sb * sb);
    if (r > 0.0) { cb /= r; sb /= r; }
  }
  for (i = 0; i < n; ++i) {
    const double ar = W->acre[i], ai = W->acim[i];
    psire[i] = cb * ar - sb * ai;
    psiim[i] = cb * ai + sb * ar;
    n2out += psire[i] * psire[i] + psiim[i] * psiim[i];
  }

  if (info) {
    info->lo = W->lo;
    info->hi = W->hi;
    info->k_used = W->k_used;
    info->c_tail = fabs(W->ckre[W->k_used]) + fabs(W->ckim[W->k_used]);
  }

  /* Rede de segurança. Um `a` pequeno demais põe o espectro de H~ fora de
     [-1,1]; os T_k crescem exponencialmente e o resultado sai errado SEM
     nenhum outro sinal. Aqui isso custa duas reduções O(N) por passo, contra
     K produtos matriz-vetor — ruído no orçamento. */
  if (n2in > 0.0) {
    const double ratio = sqrt(n2out / n2in);
    if (!(fabs(ratio - 1.0) <= DAE_NORM_TOL)) return DAE_ERR_NORM;
  }
  return DAE_OK;
}
