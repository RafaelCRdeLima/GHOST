/* t98_neutrino_vacuo.c — ACEITAÇÃO 98: oscilação de dois sabores no vácuo.
 *
 * Na base de sabor, com k = Δm²/4E = 1:
 *
 *     H = [ -cos2θ    sin2θ ]      P(νe→νμ)(t) = sin²(2θ) sin²(t)
 *         [  sin2θ    cos2θ ]
 *
 * Real e simétrica — a classe que o núcleo propaga. Os termos DIAGONAIS entram
 * como AUTO-LAÇOS na lista de arestas, e é isso que este teste trava: se
 * dae_csr_from_edges parar de aceitar laço, ou se dae_hamiltonian parar de
 * preservar a diagonal, toda a física de neutrinos deste programa cai junto.
 *
 * Como o núcleo monta H = -γA, a lista de arestas carrega -H.
 *
 * ANTI-VACUIDADE (t98b, no mesmo binário): θ=0 tem de dar P ≡ 0 e θ=45° tem de
 * dar amplitude 1. Sem esse par, um propagador que devolvesse o estado inicial
 * intacto passaria no teste principal em todo t onde sin²(t) fosse pequeno.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

/* Propaga |ν_e> = |0> até `t` e devolve P(νe→νμ) = |ψ_1|². */
static double p_transicao(dae_test *T, double theta, double t, int32_t nsteps)
{
  const double c2 = cos(2.0 * theta), s2 = sin(2.0 * theta);
  /* A = -H: o laço em 0 vale +cos2θ, o laço em 1 vale -cos2θ, a aresta -sin2θ */
  const int32_t ei[3] = { 0, 1, 0 };
  const int32_t ej[3] = { 0, 1, 1 };
  const double  w[3]  = { c2, -c2, -s2 };
  dae_csr A, H;
  dae_cheb W;
  dae_cheb_info info;
  double re[2], im[2], p = -1.0;
  int32_t s;
  dae_status st;

  st = dae_csr_from_edges(&A, 2, ei, ej, w, 3, NULL);
  dae_test_ok(T, st == DAE_OK, "from_edges: %s", dae_strerror(st));
  if (st != DAE_OK) return -1.0;

  st = dae_hamiltonian(&H, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  dae_test_ok(T, st == DAE_OK, "hamiltoniano: %s", dae_strerror(st));
  if (st != DAE_OK) { dae_csr_free(&A); return -1.0; }

  st = dae_cheb_init(&W, &H, 0);
  dae_test_ok(T, st == DAE_OK, "cheb_init: %s", dae_strerror(st));
  if (st != DAE_OK) { dae_csr_free(&A); dae_csr_free(&H); return -1.0; }

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
  const double theta = 10.0 * 3.14159265358979323846 / 180.0;
  const double s2sq = sin(2.0 * theta) * sin(2.0 * theta);
  double pior = 0.0;
  int32_t i;

  dae_test_begin(&T, "t98 oscilacao de dois sabores no vacuo");

  /* --- o teste principal: 24 pontos ao longo de dois periodos ------------ */
  for (i = 1; i <= 24; ++i) {
    const double t = (double)i * (2.0 * 3.14159265358979323846) / 12.0;
    const double got = p_transicao(&T, theta, t, 8);
    const double want = s2sq * sin(t) * sin(t);
    const double err = fabs(got - want);
    if (err > pior) pior = err;
    if (got < 0.0) break;
  }
  dae_test_ok(&T, pior < 1e-14,
              "P contra sin^2(2t) sin^2(t): pior erro %.3e (limite 1e-14)", pior);
  dae_test_note("   . theta=10 graus  sin^2(2theta)=%.12f  pior erro=%.3e", s2sq, pior);

  /* --- ANTI-VACUIDADE 1: sem mistura, nao ha transicao ------------------- */
  {
    double pmax = 0.0;
    for (i = 1; i <= 12; ++i) {
      const double p = p_transicao(&T, 0.0, (double)i * 0.7, 6);
      if (p > pmax) pmax = p;
    }
    dae_test_ok(&T, pmax < 1e-15, "theta=0 nao oscila: P_max=%.3e", pmax);
    dae_test_note("   . anti-vacuidade: theta=0  ->  P_max=%.3e", pmax);
  }

  /* --- ANTI-VACUIDADE 2: mistura maxima chega a 1 ------------------------ */
  {
    /* θ=45°: cos2θ=0, H=[[0,1],[1,0]], P = sin²(t), maximo em t=π/2 */
    const double p = p_transicao(&T, 3.14159265358979323846 / 4.0,
                                 3.14159265358979323846 / 2.0, 8);
    dae_test_near(&T, p, 1.0, 1e-14, "theta=45 graus atinge P=1 em t=pi/2");
    dae_test_note("   . anti-vacuidade: theta=45  ->  P(pi/2)=%.15f", p);
  }

  /* --- o laço diagonal REALMENTE entrou em H? ---------------------------- */
  {
    /* Sem o termo diagonal, H seria [[0,s2],[s2,0]] e a amplitude viraria 1,
       nao sin^2(2θ). Este check separa "o laço foi lido" de "o laço sumiu". */
    const double p = p_transicao(&T, theta, 3.14159265358979323846 / 2.0, 8);
    dae_test_near(&T, p, s2sq, 1e-14, "a amplitude e sin^2(2theta), nao 1");
    dae_test_ok(&T, fabs(p - 1.0) > 0.5, "o auto-laco entrou na diagonal de H");
  }

  return dae_test_end(&T);
}
