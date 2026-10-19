/* t100_mergulho.c — ACEITAÇÃO 100: três sabores COM fase CP, via mergulho real.
 *
 * `dae_csr.val` é `double`, e CONVENTIONS.md parte 1 proíbe `_Complex` no
 * núcleo. Com δ_CP ≠ 0, π a matriz PMNS torna H hermitiana COMPLEXA, e a
 * intuição diz que isso está fora. Não está:
 *
 *   H = X + iY   (X simétrica, Y ANTIssimétrica — é o que hermiticidade obriga)
 *   H~ = [ X  -Y ]   é REAL e SIMÉTRICA, porque (-Y)^T = Y
 *        [ Y   X ]
 *
 * O subespaço {(ψ, -iψ)} é invariante sob H~, e H~ age nele como H:
 *   H~ (ψ, -iψ) = (Xψ + iYψ, Yψ - iXψ) = (Hψ, -i Hψ)
 * logo exp(-iH~t)(ψ, -iψ) = (exp(-iHt)ψ, -i exp(-iHt)ψ).
 *
 * A referência não é outro propagador: é a solução espectral fechada
 *   ψ(t) = U diag(e^{-i λ_i t}) U† ψ(0),
 * montada aqui em pares de double. Comparar propagador com propagador esconde
 * erro comum aos dois; comparar com a forma fechada, não.
 *
 * ANTI-VACUIDADE — e este é o ponto do teste: um mergulho que descartasse Y em
 * silêncio passaria em TODOS os outros checks (a norma se conserva, o subespaço
 * é preservado, a energia bate). O que ele NÃO consegue fazer é reproduzir a
 * assimetria entre δ=+90° e δ=-90°, que é violação de CP. É o único check aqui
 * que enxerga a parte imaginária.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846
#define NF 3                            /* sabores */
#define NM (2 * NF)                     /* dimensão do mergulho */

typedef struct { double re[NF][NF], im[NF][NF]; } cmat3;

/* PMNS na convenção do PDG: U13 = s13 e^{-iδ}; os demais trazem e^{+iδ}. */
static void pmns(cmat3 *U, double t12, double t13, double t23, double dcp)
{
  const double c12 = cos(t12), s12 = sin(t12);
  const double c13 = cos(t13), s13 = sin(t13);
  const double c23 = cos(t23), s23 = sin(t23);
  const double cd = cos(dcp), sd = sin(dcp);   /* e^{+iδ} = cd + i sd */
  int i, j;
  for (i = 0; i < NF; ++i) for (j = 0; j < NF; ++j) { U->re[i][j] = 0.0; U->im[i][j] = 0.0; }

  U->re[0][0] = c12 * c13;
  U->re[0][1] = s12 * c13;
  U->re[0][2] =  s13 * cd;  U->im[0][2] = -s13 * sd;      /* s13 e^{-iδ} */

  U->re[1][0] = -s12 * c23 - c12 * s23 * s13 * cd;
  U->im[1][0] =            - c12 * s23 * s13 * sd;
  U->re[1][1] =  c12 * c23 - s12 * s23 * s13 * cd;
  U->im[1][1] =            - s12 * s23 * s13 * sd;
  U->re[1][2] =  s23 * c13;

  U->re[2][0] =  s12 * s23 - c12 * c23 * s13 * cd;
  U->im[2][0] =            - c12 * c23 * s13 * sd;
  U->re[2][1] = -c12 * s23 - s12 * c23 * s13 * cd;
  U->im[2][1] =            - s12 * c23 * s13 * sd;
  U->re[2][2] =  c23 * c13;
}

/* H = U diag(lam) U†,  H_ab = sum_i U_ai lam_i conj(U_bi) */
static void hamiltoniano(cmat3 *H, const cmat3 *U, const double *lam)
{
  int a, b, i;
  for (a = 0; a < NF; ++a) for (b = 0; b < NF; ++b) {
    double sr = 0.0, si = 0.0;
    for (i = 0; i < NF; ++i) {
      const double ar = U->re[a][i], ai = U->im[a][i];
      const double br = U->re[b][i], bi = -U->im[b][i];   /* conj */
      sr += lam[i] * (ar * br - ai * bi);
      si += lam[i] * (ar * bi + ai * br);
    }
    H->re[a][b] = sr; H->im[a][b] = si;
  }
}

/* Solução fechada: psi_a(t) = sum_i U_ai e^{-i lam_i t} conj(U_{b0 i}) */
static void exato(const cmat3 *U, const double *lam, int b0, double t,
                  double *pre, double *pim)
{
  int a, i;
  for (a = 0; a < NF; ++a) {
    double sr = 0.0, si = 0.0;
    for (i = 0; i < NF; ++i) {
      const double er = cos(lam[i] * t), ei = -sin(lam[i] * t);   /* e^{-i lam t} */
      const double cr = U->re[b0][i], ci = -U->im[b0][i];         /* conj(U_{b0 i}) */
      /* U_ai * e * conj(U_b0i) */
      double xr = U->re[a][i] * er - U->im[a][i] * ei;
      double xi = U->re[a][i] * ei + U->im[a][i] * er;
      sr += xr * cr - xi * ci;
      si += xr * ci + xi * cr;
    }
    pre[a] = sr; pim[a] = si;
  }
}

/* Propaga |b0> por `t` usando o MERGULHO 2n x 2n dentro do núcleo.
 * Devolve 0 em sucesso; escreve psi (3 componentes) e o desvio máximo do
 * estado em relação ao subespaço {(psi, -i psi)}. */
static int propaga(dae_test *T, const cmat3 *H, int b0, double t, int32_t nsteps,
                   double *pre, double *pim, double *fora_do_subespaco)
{
  int32_t ei[NM * NM], ej[NM * NM];
  double w[NM * NM];
  int32_t ne = 0, a, b;
  dae_csr A, Hc;
  dae_cheb W;
  dae_cheb_info info;
  double re[NM], im[NM];
  int32_t s;
  dae_status st;
  int rc = 1;

  /* H~ = [[X, -Y], [Y, X]]; o núcleo monta H = -A, entao A = -H~.
     Só o triângulo superior: dae_csr_from_edges simetriza. */
  for (a = 0; a < NM; ++a) for (b = a; b < NM; ++b) {
    double v;
    if      (a <  NF && b <  NF) v =  H->re[a][b];        /* X  */
    else if (a <  NF && b >= NF) v = -H->im[a][b - NF];   /* -Y */
    else                          v =  H->re[a - NF][b - NF];
    v = -v;                                              /* A = -H~ */
    if (fabs(v) > 1e-300) { ei[ne] = a; ej[ne] = b; w[ne] = v; ++ne; }
  }

  st = dae_csr_from_edges(&A, NM, ei, ej, w, ne, NULL);
  if (st != DAE_OK) { dae_test_ok(T, 0, "from_edges: %s", dae_strerror(st)); return 1; }
  st = dae_hamiltonian(&Hc, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) { dae_test_ok(T, 0, "H: %s", dae_strerror(st)); dae_csr_free(&A); return 1; }
  st = dae_cheb_init(&W, &Hc, 0);
  if (st != DAE_OK) { dae_test_ok(T, 0, "cheb: %s", dae_strerror(st));
                      dae_csr_free(&A); dae_csr_free(&Hc); return 1; }

  /* phi0 = (psi0, -i psi0)/sqrt(2), com psi0 = e_{b0} */
  for (a = 0; a < NM; ++a) { re[a] = 0.0; im[a] = 0.0; }
  re[b0] = 1.0 / sqrt(2.0);
  im[NF + b0] = -1.0 / sqrt(2.0);

  for (s = 0; s < nsteps; ++s) {
    st = dae_cheb_step(&W, &Hc, t / (double)nsteps, re, im, &info);
    if (st != DAE_OK) { dae_test_ok(T, 0, "passo: %s", dae_strerror(st)); goto fim; }
  }

  /* O estado ainda está em {(psi, -i psi)}?  phi_{n+a} tem de valer -i phi_a. */
  *fora_do_subespaco = 0.0;
  for (a = 0; a < NF; ++a) {
    const double dr = fabs(re[NF + a] - im[a]);      /* Re(-i z) =  Im(z) */
    const double di = fabs(im[NF + a] + re[a]);      /* Im(-i z) = -Re(z) */
    if (dr > *fora_do_subespaco) *fora_do_subespaco = dr;
    if (di > *fora_do_subespaco) *fora_do_subespaco = di;
  }

  /* projeta de volta: psi_a = (phi_a + i phi_{n+a})/sqrt(2) */
  for (a = 0; a < NF; ++a) {
    pre[a] = (re[a] - im[NF + a]) / sqrt(2.0);
    pim[a] = (im[a] + re[NF + a]) / sqrt(2.0);
  }
  rc = 0;

fim:
  dae_cheb_free(&W);
  dae_csr_free(&Hc);
  dae_csr_free(&A);
  return rc;
}

/* P(nu_mu -> nu_e) em t, para uma fase CP dada. */
static double p_mu_e(dae_test *T, double dcp, double t)
{
  static const double lam[NF] = { 0.0, 1.0, 33.8 };   /* dm21^2 : dm31^2 */
  cmat3 U, H;
  double pre[NF], pim[NF], fora = 0.0;
  pmns(&U, 33.4 * PI / 180.0, 8.6 * PI / 180.0, 49.0 * PI / 180.0, dcp);
  hamiltoniano(&H, &U, lam);
  if (propaga(T, &H, 1, t, 40, pre, pim, &fora) != 0) return -1.0;
  return pre[0] * pre[0] + pim[0] * pim[0];
}

int main(void)
{
  dae_test T;
  static const double lam[NF] = { 0.0, 1.0, 33.8 };
  const double dcp = -PI / 2.0;
  cmat3 U, H;
  double pior_amp = 0.0, pior_sub = 0.0;
  int32_t i;

  dae_test_begin(&T, "t100 tres sabores com fase CP, via mergulho real 2n x 2n");

  pmns(&U, 33.4 * PI / 180.0, 8.6 * PI / 180.0, 49.0 * PI / 180.0, dcp);
  hamiltoniano(&H, &U, lam);

  /* --- H e mesmo hermitiana, e Y e mesmo antissimetrica? ---------------- */
  {
    double herm = 0.0;
    int a, b;
    for (a = 0; a < NF; ++a) for (b = 0; b < NF; ++b) {
      const double dr = fabs(H.re[a][b] - H.re[b][a]);
      const double di = fabs(H.im[a][b] + H.im[b][a]);
      if (dr > herm) herm = dr;
      if (di > herm) herm = di;
    }
    dae_test_ok(&T, herm < 1e-15, "H hermitiana (X simetrica, Y antissimetrica): %.3e", herm);
  }

  /* --- o teste principal: amplitude contra a solucao espectral fechada --- */
  for (i = 1; i <= 20; ++i) {
    const double t = 0.15 * (double)i;
    double pre[NF], pim[NF], ere[NF], eim[NF], fora = 0.0;
    int a;
    if (propaga(&T, &H, 1, t, 40, pre, pim, &fora) != 0) break;
    exato(&U, lam, 1, t, ere, eim);
    if (fora > pior_sub) pior_sub = fora;
    for (a = 0; a < NF; ++a) {
      const double d = sqrt((pre[a] - ere[a]) * (pre[a] - ere[a])
                          + (pim[a] - eim[a]) * (pim[a] - eim[a]));
      if (d > pior_amp) pior_amp = d;
    }
  }
  dae_test_ok(&T, pior_amp < 1e-13,
              "amplitude contra U diag(e^-ilt) U^dag: pior %.3e (limite 1e-13)", pior_amp);
  dae_test_ok(&T, pior_sub < 1e-13,
              "o estado nao sai de {(psi,-i psi)}: desvio %.3e", pior_sub);
  dae_test_note("   . erro na amplitude %.3e   desvio do subespaco %.3e", pior_amp, pior_sub);

  /* --- ANTI-VACUIDADE: a fase CP e carregada, e viola CP ---------------- */
  {
    const double t = 3.0;
    const double p0   = p_mu_e(&T, 0.0,      t);
    const double pm90 = p_mu_e(&T, -PI / 2.0, t);
    const double pp90 = p_mu_e(&T,  PI / 2.0, t);
    const double p180 = p_mu_e(&T, PI,        t);

    dae_test_ok(&T, fabs(pm90 - pp90) > 1e-6,
                "P(delta=-90) != P(delta=+90): %.9f vs %.9f  (assimetria %.3e)",
                pm90, pp90, fabs(pm90 - pp90));
    dae_test_ok(&T, fabs(p0 - p180) > 1e-6,
                "P(delta=0) != P(delta=180): %.9f vs %.9f", p0, p180);
    dae_test_note("   . anti-vacuidade CP:  0=%.6f  -90=%.6f  +90=%.6f  180=%.6f",
                  p0, pm90, pp90, p180);
    dae_test_note("   . assimetria |P(-90)-P(+90)| = %.6e  (um mergulho sem Y daria 0)",
                  fabs(pm90 - pp90));
  }

  /* --- ANTI-VACUIDADE 2: com delta=0 a parte imaginaria some ------------ */
  {
    cmat3 U0, H0;
    double maxim = 0.0;
    int a, b;
    pmns(&U0, 33.4 * PI / 180.0, 8.6 * PI / 180.0, 49.0 * PI / 180.0, 0.0);
    hamiltoniano(&H0, &U0, lam);
    for (a = 0; a < NF; ++a) for (b = 0; b < NF; ++b)
      if (fabs(H0.im[a][b]) > maxim) maxim = fabs(H0.im[a][b]);
    dae_test_ok(&T, maxim < 1e-15, "delta=0 da H real: max|Im H| = %.3e", maxim);
    dae_test_note("   . com delta=0 o mergulho degenera para dois blocos reais iguais");
  }

  return dae_test_end(&T);
}
