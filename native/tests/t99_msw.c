/* t99_msw.c — ACEITAÇÃO 99: ressonância MSW em matéria de densidade constante.
 *
 * O potencial de Wolfenstein V = √2 G_F n_e afeta só ν_e. Com
 * A ≡ 2√2 G_F n_e E/Δm² e o traço removido:
 *
 *     H = [ -(cos2θ − A)      sin2θ     ]
 *         [   sin2θ        (cos2θ − A)  ]
 *
 * A oscilação vira P(t) = (sin²2θ/ω²) sin²(ω t), com ω = √[(cos2θ−A)² + sin²2θ].
 * A amplitude sin²2θ/ω² é MÁXIMA em A = cos2θ, onde a diagonal se anula: é a
 * ressonância, e o ponto físico é que ela chega perto de 1 com θ pequeno e
 * INALTERADO — a mistura grande vem da matéria, não do ângulo de vácuo.
 *
 * O teste avalia no máximo exato, t = π/(2ω), em vez de procurar o pico numa
 * grade. Procurar numa grade mediria a resolução da grade; avaliar no máximo
 * analítico mede o propagador.
 *
 * ANTI-VACUIDADE: a ressonância tem de SE MOVER com θ. Um código que
 * devolvesse "amplitude alta perto de A=1" para qualquer θ passaria num teste
 * de um ângulo só.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

/* Amplitude máxima de P(νe→νμ) para (θ, A), avaliada em t = π/(2ω). */
static double amplitude(dae_test *T, double theta, double a_msw, int32_t nsteps)
{
  const double c2 = cos(2.0 * theta), s2 = sin(2.0 * theta);
  const double d  = c2 - a_msw;
  const double w2 = d * d + s2 * s2;
  const double om = sqrt(w2);
  /* A = -H */
  const int32_t ei[3] = { 0, 1, 0 };
  const int32_t ej[3] = { 0, 1, 1 };
  const double  wt[3] = { d, -d, -s2 };
  dae_csr A, H;
  dae_cheb W;
  dae_cheb_info info;
  double re[2], im[2], p = -1.0, t;
  int32_t s;
  dae_status st;

  if (!(om > 0.0)) return -1.0;
  t = PI / (2.0 * om);                  /* o máximo exato */

  st = dae_csr_from_edges(&A, 2, ei, ej, wt, 3, NULL);
  if (st != DAE_OK) { dae_test_ok(T, 0, "from_edges: %s", dae_strerror(st)); return -1.0; }
  st = dae_hamiltonian(&H, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) { dae_test_ok(T, 0, "H: %s", dae_strerror(st)); dae_csr_free(&A); return -1.0; }
  st = dae_cheb_init(&W, &H, 0);
  if (st != DAE_OK) { dae_test_ok(T, 0, "cheb: %s", dae_strerror(st));
                      dae_csr_free(&A); dae_csr_free(&H); return -1.0; }

  re[0] = 1.0; re[1] = 0.0; im[0] = 0.0; im[1] = 0.0;
  for (s = 0; s < nsteps; ++s) {
    st = dae_cheb_step(&W, &H, t / (double)nsteps, re, im, &info);
    if (st != DAE_OK) { dae_test_ok(T, 0, "passo: %s", dae_strerror(st)); goto fim; }
  }
  p = re[1] * re[1] + im[1] * im[1];

fim:
  dae_cheb_free(&W);
  dae_csr_free(&H);
  dae_csr_free(&A);
  return p;
}

int main(void)
{
  dae_test T;
  const double theta = 10.0 * PI / 180.0;
  const double c2 = cos(2.0 * theta), s2sq = sin(2.0 * theta) * sin(2.0 * theta);
  double pior = 0.0;
  int32_t i;

  dae_test_begin(&T, "t99 ressonancia MSW em materia de densidade constante");

  /* --- o teste principal: varre A e compara com sin^2(2θ)/ω^2 ------------ */
  for (i = 0; i <= 40; ++i) {
    const double a = 2.0 * (double)i / 40.0;
    const double d = c2 - a;
    const double want = s2sq / (d * d + s2sq);
    const double got = amplitude(&T, theta, a, 12);
    double err;
    if (got < 0.0) break;
    err = fabs(got - want);
    if (err > pior) pior = err;
  }
  dae_test_ok(&T, pior < 1e-12,
              "amplitude contra sin^2(2t)/[(cos2t-A)^2+sin^2(2t)]: pior %.3e", pior);

  /* --- a ressonância existe: amplitude ~1 com theta pequeno -------------- */
  {
    const double na_res = amplitude(&T, theta, c2, 12);
    dae_test_near(&T, na_res, 1.0, 1e-12, "em A=cos2theta a amplitude e 1");
    dae_test_ok(&T, s2sq < 0.2, "e o angulo de vacuo continua pequeno (%.4f)", s2sq);
    dae_test_note("   . theta=10 graus: vacuo %.6f  ->  ressonancia %.12f", s2sq, na_res);
  }

  /* --- ANTI-VACUIDADE: a ressonancia SE MOVE com theta ------------------- */
  {
    static const double graus[3] = { 5.0, 15.0, 30.0 };
    int32_t k;
    for (k = 0; k < 3; ++k) {
      const double th = graus[k] * PI / 180.0;
      const double alvo = cos(2.0 * th);
      double melhor = -1.0, a_melhor = -1.0;
      for (i = 0; i <= 80; ++i) {
        const double a = 2.0 * (double)i / 80.0;
        const double p = amplitude(&T, th, a, 12);
        if (p > melhor) { melhor = p; a_melhor = a; }
      }
      /* a grade tem passo 0,025; o pico tem de cair na celula certa */
      dae_test_ok(&T, fabs(a_melhor - alvo) <= 0.025 + 1e-12,
                  "theta=%.0f: pico em A=%.3f, esperado cos2theta=%.3f",
                  graus[k], a_melhor, alvo);
      dae_test_note("   . theta=%2.0f graus  pico em A=%.3f  (cos2theta=%.3f)  amp=%.6f",
                    graus[k], a_melhor, alvo, melhor);
    }
  }

  /* --- ANTI-VACUIDADE 2: longe da ressonancia a materia SUPRIME ---------- */
  {
    const double longe = amplitude(&T, theta, 2.0, 12);
    dae_test_ok(&T, longe < s2sq,
                "A=2 suprime abaixo do vacuo: %.6f < %.6f", longe, s2sq);
    dae_test_note("   . anti-vacuidade: A=2 (longe) -> amplitude %.6f, abaixo do vacuo", longe);
  }

  return dae_test_end(&T);
}
