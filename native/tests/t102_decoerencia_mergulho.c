/* t102_decoerencia_mergulho.c — ACEITAÇÃO 102: defasagem COM o mergulho, e a
 * defasagem dependente do tempo.
 *
 * O IMPASSE QUE ISTO DESFAZ. Mergulho e trajetórias parecem incompatíveis: o
 * colapso |j><j| no espaço de 2n projeta numa componente só e destrói a forma
 * (psi, -i psi), jogando o estado para fora do subespaço invariante. O
 * resultado deixaria de significar coisa alguma, e com isso "três sabores com
 * fase CP" e "decoerência" seriam mutuamente exclusivos.
 *
 * A saída é geométrica, não numérica: o colapso não deve ser num sítio, e sim
 * no PAR {j, nf+j}. Então
 *     phi' = (psi_j e_j , -i psi_j e_j)/sqrt(2)
 * que É da forma certa, e a probabilidade sai de graça:
 *     |phi_j|^2 + |phi_{nf+j}|^2 = |psi_j|^2.
 *
 * O PORTÃO. Com delta_cp = 0 o Hamiltoniano é real e o mesmo problema pode ser
 * rodado dos DOIS jeitos: em 3 sítios com colapso em sítio (o caminho padrão
 * das trajetórias) e em 6 sítios mergulhados com colapso no par. Os dois
 * consomem o PRNG na mesma ordem e na mesma quantidade — um sorteio para o
 * tempo de espera, um para o colapso — então com a mesma semente as
 * trajetórias são as MESMAS, e a concordância tem de ser de máquina, não
 * estatística. É um teste muito mais duro que comparar médias de ensemble.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define NF 3

/* gamma vive na MESMA unidade de H — eV —, porque o tempo do propagador e
 * 1/eV. Com dt = 5e10 e nt = 200 o caminho dura T = 1e13, entao gamma = 5e-13
 * da ~5 saltos por trajetoria. Um gamma de aparencia inocente como 0.02 daria
 * 2e11 saltos e o teste nao terminaria.
 */
#define GAMA 5.0e-13

/* Roda o ensemble e devolve rho de sabor na ultima amostra.
 * `mergulhar` = 0: 3 sitios, colapso em sitio.
 * `mergulhar` = 1: 6 sitios (H~ = diag(X, X), porque Y = 0), colapso no par. */
static dae_status roda(dae_test *T, int mergulhar, double gamma_const,
                       const double *gamma_tab, int32_t n_traj,
                       double *rho_re, double *rho_im)
{
  dae_nu_params P;
  double hre[NF * NF], him[NF * NF];
  int32_t ei[8 * NF * NF], ej[8 * NF * NF];
  double w[8 * NF * NF];
  double init_re[2 * NF], init_im[2 * NF];
  int32_t ne = 0, n, a, b;
  dae_csr A, H;
  dae_cheb W;
  dae_traj_cfg cfg;
  dae_rho_acc acc;
  dae_status st;

  dae_nu_params_default(&P);
  P.delta_cp = 0.0;                 /* H real: o mesmo caso serve aos dois */
  P.energy_eV = 20.0e6;
  st = dae_nu_hamiltoniano(&P, 100.0, 0.5, hre, him);
  if (st != DAE_OK) return st;

  n = mergulhar ? 2 * NF : NF;
  for (a = 0; a < n; ++a) for (b = a; b < n; ++b) {
    double v;
    if (!mergulhar)                 v = hre[a * NF + b];
    else if (a < NF && b < NF)      v = hre[a * NF + b];
    else if (a < NF)                v = 0.0;            /* Y = 0 */
    else                            v = hre[(a - NF) * NF + (b - NF)];
    v = -v;                                             /* A = -H */
    if (v != 0.0) { ei[ne] = a; ej[ne] = b; w[ne] = v; ++ne; }
  }
  st = dae_csr_from_edges(&A, n, ei, ej, w, ne, NULL);
  if (st != DAE_OK) return st;
  st = dae_hamiltonian(&H, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) { dae_csr_free(&A); return st; }
  st = dae_cheb_init(&W, &H, 0);
  if (st != DAE_OK) { dae_csr_free(&A); dae_csr_free(&H); return st; }

  dae_nu_inicial(1, NF, mergulhar, init_re, init_im);   /* nu_mu */

  memset(&cfg, 0, sizeof(cfg));
  cfg.gamma_deph = gamma_const;
  cfg.gamma_tab = gamma_tab;
  cfg.n_traj = n_traj;
  cfg.rho_stride = 10;
  cfg.nt = 200;
  cfg.dt = 5.0e10;                  /* ~10 km em 1/eV */
  cfg.sitio_inicial = 0;
  cfg.init_re = init_re;
  cfg.init_im = init_im;
  cfg.nf_par = mergulhar ? NF : 0;
  cfg.saida = DAE_SAIDA_ACUMULAR_RHO;

  st = dae_rho_acc_init(&acc, n, dae_traj_amostras(&cfg));
  if (st != DAE_OK) { dae_cheb_free(&W); dae_csr_free(&A); dae_csr_free(&H); return st; }

  st = dae_traj_ensemble(&cfg, 4242ULL, &W, &H, &acc, NULL, NULL);
  if (st == DAE_OK) {
    const int32_t ult = acc.n_amostras - 1;
    dae_rho_acc_finalizar(&acc);
    for (a = 0; a < NF; ++a) for (b = 0; b < NF; ++b) {
      const size_t k = (size_t)ult * (size_t)n * (size_t)n
                     + (size_t)a * (size_t)n + (size_t)b;
      /* No mergulho, rho~_ab = rho_ab / 2 (porque phi_a = psi_a/sqrt(2)). */
      const double f = mergulhar ? 2.0 : 1.0;
      rho_re[a * NF + b] = acc.re[k] * f;
      rho_im[a * NF + b] = acc.im[k] * f;
    }
  } else {
    dae_test_ok(T, 0, "ensemble: %s", dae_strerror(st));
  }

  dae_rho_acc_free(&acc);
  dae_cheb_free(&W); dae_csr_free(&A); dae_csr_free(&H);
  return st;
}

int main(void)
{
  dae_test T;
  double r1[NF * NF], i1[NF * NF], r2[NF * NF], i2[NF * NF];
  double gtab[200];
  int32_t k;
  dae_status st;

  dae_test_begin(&T, "t102 defasagem com mergulho, e defasagem em gamma(t)");

  /* ---------------- O PORTÃO: os dois caminhos, bit a bit ---------------- */
  st = roda(&T, 0, GAMA, NULL, 400, r1, i1);
  dae_test_ok(&T, st == DAE_OK, "3 sitios, colapso em sitio: %s", dae_strerror(st));
  if (st == DAE_OK) {
    st = roda(&T, 1, GAMA, NULL, 400, r2, i2);
    dae_test_ok(&T, st == DAE_OK, "6 sitios mergulhados, colapso no par: %s",
                dae_strerror(st));
  }
  if (st == DAE_OK) {
    double pior = 0.0;
    int32_t a;
    for (a = 0; a < NF * NF; ++a) {
      const double d = fabs(r1[a] - r2[a]) + fabs(i1[a] - i2[a]);
      if (d > pior) pior = d;
    }
    dae_test_ok(&T, pior < 1e-14,
                "PORTAO: os dois caminhos dao o MESMO rho (%.3e)", pior);
    dae_test_note("   . portao sitio-vs-par: diferenca maxima em rho = %.3e", pior);
    dae_test_note("   . (mesma semente, mesmo consumo de PRNG -> concordancia de MAQUINA)");
  }

  /* ---- ANTI-VACUIDADE 1: o portao enxerga se o colapso for no lugar errado?
     Sem defasagem (gamma = 0) nao ha colapso nenhum, e os dois caminhos
     coincidiriam de qualquer jeito. Se o teste so rodasse assim, ele nao
     provaria nada sobre o colapso. Confere que gamma > 0 REALMENTE decoere. */
  {
    double rc[NF * NF], ic[NF * NF];
    st = roda(&T, 1, 0.0, NULL, 200, rc, ic);
    if (st == DAE_OK) {
      const double coe_com = sqrt(r2[1] * r2[1] + i2[1] * i2[1]);
      const double coe_sem = sqrt(rc[1] * rc[1] + ic[1] * ic[1]);
      dae_test_ok(&T, coe_sem > coe_com * 1.5,
                  "gamma>0 destroi coerencia: |rho_01| %.6f -> %.6f", coe_sem, coe_com);
      dae_test_note("   . anti-vacuidade: |rho_01| sem defasagem %.6f, com %.6f",
                    coe_sem, coe_com);
    }
  }

  /* ---- ANTI-VACUIDADE 2: o estado FICA no subespaco?
     A populacao total tem de somar 1 depois de centenas de colapsos. Se o
     colapso jogasse o estado para fora de {(psi,-i psi)}, o traco do rho
     PROJETADO deixaria de ser 1 — e essa e a assinatura direta do impasse
     que este arquivo desfaz. */
  {
    double tr = r2[0] + r2[4] + r2[8];
    dae_test_near(&T, tr, 1.0, 1e-12,
                  "traco do rho de sabor apos os colapsos (o subespaco sobreviveu)");
    dae_test_note("   . traco = %.15f  -> o colapso no par nao vazou do subespaco", tr);
  }

  /* ---------------------- gamma(t): constante por tabela == constante ---- */
  for (k = 0; k < 200; ++k) gtab[k] = GAMA;
  st = roda(&T, 1, 0.0, gtab, 400, r1, i1);
  if (st == DAE_OK) {
    double pior = 0.0;
    int32_t a;
    for (a = 0; a < NF * NF; ++a) {
      const double d = fabs(r1[a] - r2[a]) + fabs(i1[a] - i2[a]);
      if (d > pior) pior = d;
    }
    dae_test_ok(&T, pior < 1e-12,
                "gamma_tab constante reproduz gamma_deph constante (%.3e)", pior);
  }

  /* ---------------------- gamma(t) VARIAVEL faz diferenca ---------------- */
  {
    double rv[NF * NF], iv[NF * NF];
    /* mesma defasagem total, concentrada na segunda metade */
    for (k = 0; k < 200; ++k) gtab[k] = (k < 100) ? 0.0 : 2.0 * GAMA;
    st = roda(&T, 1, 0.0, gtab, 400, rv, iv);
    if (st == DAE_OK) {
      const double c_uni = sqrt(r2[1] * r2[1] + i2[1] * i2[1]);
      const double c_var = sqrt(rv[1] * rv[1] + iv[1] * iv[1]);
      dae_test_ok(&T, fabs(c_uni - c_var) > 1e-6,
                  "gamma(t) concentrada difere de uniforme: %.6f contra %.6f",
                  c_var, c_uni);
      dae_test_note("   . gamma(t) uniforme -> |rho_01| = %.6f", c_uni);
      dae_test_note("   . gamma(t) so na 2a metade, mesma dose -> %.6f", c_var);
      dae_test_note("   . (mesma defasagem TOTAL: se nao diferisse, a tabela seria ignorada)");
    }
  }

  return dae_test_end(&T);
}
