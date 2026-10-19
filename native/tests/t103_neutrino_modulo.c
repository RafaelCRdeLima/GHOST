/* t103_neutrino_modulo.c — PORTÃO DO MÓDULO DE NEUTRINOS: dae_neutrino contra
 * tudo o que já se sabe.
 *
 * A pergunta que este teste responde não é "o módulo roda", é "o módulo dá os
 * MESMOS números que a edgelist montada à mão". Se divergir, o módulo alterou
 * a física, e isso tem de aparecer aqui e não numa figura.
 *
 * Quatro frentes, e nenhuma delas confia nas outras:
 *   1. contra a forma fechada de dois sabores, sin^2(2t) sin^2(kL);
 *   2. contra a amplitude MSW, sin^2(2t)/[(cos2t - A)^2 + sin^2(2t)];
 *   3. contra a solução espectral U diag(e^{-i lam t}) U^dag, com fase CP;
 *   4. contra a EDGELIST montada à mão — o caminho que t98/t99/t100 travaram.
 *
 * A quarta é o portão de verdade: as outras três poderiam passar com um módulo
 * que reimplementasse a física de outro jeito e concordasse com a teoria mas
 * não com o resto do programa.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846
#define MAXN 3

/* Propaga |sabor> por `t` (em 1/eV) usando o MÓDULO, e devolve psi projetado. */
static int roda_modulo(dae_test *T, const dae_nu_params *P, double rho, double ye,
                       int32_t sabor, double t, int32_t nsteps,
                       double *pre, double *pim)
{
  double hre[MAXN * MAXN], him[MAXN * MAXN];
  double re[2 * MAXN], im[2 * MAXN];
  dae_csr A, H;
  dae_cheb W;
  dae_cheb_info info;
  int merg = 0;
  int32_t s;
  dae_status st;
  int rc = 1;

  st = dae_nu_hamiltoniano(P, rho, ye, hre, him);
  if (st != DAE_OK) { dae_test_ok(T, 0, "H: %s", dae_strerror(st)); return 1; }
  st = dae_nu_adjacencia(hre, him, P->n_flavors, &A, &merg);
  if (st != DAE_OK) { dae_test_ok(T, 0, "adjacencia: %s", dae_strerror(st)); return 1; }
  st = dae_hamiltonian(&H, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) { dae_test_ok(T, 0, "ham: %s", dae_strerror(st)); dae_csr_free(&A); return 1; }
  st = dae_cheb_init(&W, &H, 0);
  if (st != DAE_OK) { dae_test_ok(T, 0, "cheb: %s", dae_strerror(st));
                      dae_csr_free(&A); dae_csr_free(&H); return 1; }

  dae_nu_inicial(sabor, P->n_flavors, merg, re, im);
  for (s = 0; s < nsteps; ++s) {
    st = dae_cheb_step(&W, &H, t / (double)nsteps, re, im, &info);
    if (st != DAE_OK) { dae_test_ok(T, 0, "passo: %s", dae_strerror(st)); goto fim; }
  }
  dae_nu_projeta(re, im, P->n_flavors, merg, pre, pim);
  rc = 0;

fim:
  dae_cheb_free(&W); dae_csr_free(&H); dae_csr_free(&A);
  return rc;
}

/* O MESMO cenário, montado à mão como edgelist — o caminho de t98/t99/t100. */
static int roda_edgelist(dae_test *T, const double *hre, const double *him,
                         int32_t nf, int32_t sabor, double t, int32_t nsteps,
                         double *pre, double *pim)
{
  int32_t ei[4 * MAXN * MAXN], ej[4 * MAXN * MAXN];
  double w[4 * MAXN * MAXN];
  int32_t ne = 0, n, a, b, s;
  int merg = 0;
  double re[2 * MAXN], im[2 * MAXN];
  dae_csr A, H;
  dae_cheb W;
  dae_cheb_info info;
  dae_status st;
  int rc = 1;

  for (a = 0; a < nf * nf; ++a) if (him[a] != 0.0) merg = 1;
  n = merg ? 2 * nf : nf;

  for (a = 0; a < n; ++a) for (b = a; b < n; ++b) {
    double v;
    if (!merg)                    v = hre[a * nf + b];
    else if (a < nf && b < nf)    v = hre[a * nf + b];
    else if (a < nf)              v = -him[a * nf + (b - nf)];
    else                          v = hre[(a - nf) * nf + (b - nf)];
    v = -v;
    if (v != 0.0) { ei[ne] = a; ej[ne] = b; w[ne] = v; ++ne; }
  }

  st = dae_csr_from_edges(&A, n, ei, ej, w, ne, NULL);
  if (st != DAE_OK) { dae_test_ok(T, 0, "edgelist: %s", dae_strerror(st)); return 1; }
  st = dae_hamiltonian(&H, &A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, NULL);
  if (st != DAE_OK) { dae_csr_free(&A); return 1; }
  st = dae_cheb_init(&W, &H, 0);
  if (st != DAE_OK) { dae_csr_free(&A); dae_csr_free(&H); return 1; }

  for (a = 0; a < n; ++a) { re[a] = 0.0; im[a] = 0.0; }
  if (!merg) re[sabor] = 1.0;
  else { re[sabor] = 1.0 / sqrt(2.0); im[nf + sabor] = -1.0 / sqrt(2.0); }

  for (s = 0; s < nsteps; ++s) {
    st = dae_cheb_step(&W, &H, t / (double)nsteps, re, im, &info);
    if (st != DAE_OK) goto fim;
  }
  for (a = 0; a < nf; ++a) {
    if (!merg) { pre[a] = re[a]; pim[a] = im[a]; }
    else {
      pre[a] = (re[a] - im[nf + a]) / sqrt(2.0);
      pim[a] = (im[a] + re[nf + a]) / sqrt(2.0);
    }
  }
  rc = 0;

fim:
  dae_cheb_free(&W); dae_csr_free(&H); dae_csr_free(&A);
  return rc;
}

int main(void)
{
  dae_test T;
  dae_nu_params P;
  double pre[MAXN], pim[MAXN];
  double pior;
  int32_t i;

  dae_test_begin(&T, "t103 dae_neutrino: o modulo contra a teoria e contra a edgelist");

  /* ------------------------------------------------ 1. dois sabores, vacuo */
  {
    const double theta = 10.0 * PI / 180.0;
    double s2sq, k;
    dae_nu_params_default(&P);
    P.n_flavors = 2;
    P.theta12 = theta;
    P.dm21_sq = 7.42e-5;
    P.energy_eV = 20.0e6;
    s2sq = sin(2.0 * theta) * sin(2.0 * theta);
    /* Com o traço removido, H = k[[-cos2t, sin2t],[sin2t, cos2t]] e
       k = Dm2/4E — metade do autovalor Dm2/2E, porque o traço tirou a média. */
    k = P.dm21_sq / (4.0 * P.energy_eV);

    pior = 0.0;
    for (i = 1; i <= 16; ++i) {
      const double t = (double)i * PI / (8.0 * k);
      double want, err;
      if (roda_modulo(&T, &P, 0.0, 0.0, 0, t, 12, pre, pim) != 0) break;
      want = s2sq * sin(k * t) * sin(k * t);
      err = fabs(pre[1] * pre[1] + pim[1] * pim[1] - want);
      if (err > pior) pior = err;
    }
    dae_test_ok(&T, pior < 1e-13, "2 sabores no vacuo: pior erro %.3e", pior);
    dae_test_note("   . 2 sabores, vacuo: pior erro %.3e   sin^2(2theta)=%.9f", pior, s2sq);
  }

  /* ------------------------------------------- 2. dois sabores em materia */
  {
    const double theta = 10.0 * PI / 180.0;
    const double c2 = cos(2.0 * theta), s2 = sin(2.0 * theta);
    double k, rho_res, melhor = -1.0, rho_melhor = -1.0;
    dae_nu_params_default(&P);
    P.n_flavors = 2; P.theta12 = theta;
    P.dm21_sq = 7.42e-5; P.energy_eV = 20.0e6;
    k = P.dm21_sq / (4.0 * P.energy_eV);

    /* Ressonância: V = Dm2 cos2t / 2E  =>  rho = ... */
    rho_res = (P.dm21_sq * c2 / (2.0 * P.energy_eV)) / (DAE_NU_V_POR_RHO * 0.5);

    for (i = 0; i <= 60; ++i) {
      const double rho = rho_res * 2.0 * (double)i / 60.0;
      /* amplitude máxima: avalia no máximo exato t = pi/(2 w) */
      const double a = (DAE_NU_V_POR_RHO * 0.5 * rho) / (2.0 * k);
      const double d = c2 - a;
      const double om = k * sqrt(d * d + s2 * s2);
      double p;
      if (!(om > 0.0)) continue;
      if (roda_modulo(&T, &P, rho, 0.5, 0, PI / (2.0 * om), 16, pre, pim) != 0) break;
      p = pre[1] * pre[1] + pim[1] * pim[1];
      if (p > melhor) { melhor = p; rho_melhor = rho; }
    }
    dae_test_ok(&T, melhor > 0.99, "a ressonancia MSW aparece: amplitude %.6f", melhor);
    dae_test_ok(&T, fabs(rho_melhor - rho_res) < 0.05 * rho_res,
                "e esta em rho=%.1f g/cm3, esperado %.1f", rho_melhor, rho_res);
    dae_test_note("   . MSW: vacuo %.6f -> ressonancia %.6f em rho=%.1f g/cm3 (teoria %.1f)",
                  s2 * s2, melhor, rho_melhor, rho_res);
  }

  /* ------------------- 3. tres sabores com fase CP, contra a forma fechada */
  {
    double ure[9], uim[9], lam[3];
    double hre[9], him[9];
    int merg = 0;
    dae_nu_params_default(&P);
    P.delta_cp = -PI / 2.0;
    P.energy_eV = 20.0e6;

    dae_nu_pmns(&P, ure, uim);
    lam[0] = 0.0;
    lam[1] = P.dm21_sq / (2.0 * P.energy_eV);
    lam[2] = P.dm31_sq / (2.0 * P.energy_eV);
    dae_nu_hamiltoniano(&P, 0.0, 0.0, hre, him);
    dae_nu_precisa_mergulho(him, 3);
    merg = dae_nu_precisa_mergulho(him, 3);
    dae_test_ok(&T, merg == 1, "delta_cp=-90 exige mergulho");

    pior = 0.0;
    for (i = 1; i <= 12; ++i) {
      /* t em 1/eV, escalado pela oscilação atmosférica */
      const double t = (double)i * 0.25 / lam[2];
      double ere[3], eim[3];
      int32_t a, m;
      if (roda_modulo(&T, &P, 0.0, 0.0, 1, t, 24, pre, pim) != 0) break;
      /* solução fechada: psi_a = sum_m U_am e^{-i lam_m t} conj(U_1m) */
      for (a = 0; a < 3; ++a) {
        double sr = 0.0, si = 0.0;
        for (m = 0; m < 3; ++m) {
          /* o traço foi removido de H, entao a fase global difere; compara
             so o MODULO das amplitudes, que e o observavel */
          const double er = cos(lam[m] * t), eiv = -sin(lam[m] * t);
          const double cr = ure[1 * 3 + m], ci = -uim[1 * 3 + m];
          const double xr = ure[a * 3 + m] * er - uim[a * 3 + m] * eiv;
          const double xi = ure[a * 3 + m] * eiv + uim[a * 3 + m] * er;
          sr += xr * cr - xi * ci;
          si += xr * ci + xi * cr;
        }
        ere[a] = sr; eim[a] = si;
      }
      for (a = 0; a < 3; ++a) {
        const double pg = pre[a] * pre[a] + pim[a] * pim[a];
        const double pe = ere[a] * ere[a] + eim[a] * eim[a];
        if (fabs(pg - pe) > pior) pior = fabs(pg - pe);
      }
    }
    dae_test_ok(&T, pior < 1e-13,
                "3 sabores com CP contra U diag(e^-ilt) U^dag: pior %.3e", pior);
    dae_test_note("   . 3 sabores com delta_cp=-90: pior erro %.3e", pior);
  }

  /* --------------- 4. O PORTÃO: modulo contra a edgelist montada a mao ---- */
  {
    double hre[9], him[9];
    double ere[MAXN], eim[MAXN];
    double t;
    dae_nu_params_default(&P);
    P.delta_cp = -PI / 2.0;
    P.energy_eV = 20.0e6;
    dae_nu_hamiltoniano(&P, 1600.0, 0.5, hre, him);
    t = 2.0 / (P.dm31_sq / (2.0 * P.energy_eV));

    pior = 0.0;
    if (roda_modulo(&T, &P, 1600.0, 0.5, 1, t, 24, pre, pim) == 0 &&
        roda_edgelist(&T, hre, him, 3, 1, t, 24, ere, eim) == 0) {
      int32_t a;
      for (a = 0; a < 3; ++a) {
        const double d = sqrt((pre[a] - ere[a]) * (pre[a] - ere[a])
                            + (pim[a] - eim[a]) * (pim[a] - eim[a]));
        if (d > pior) pior = d;
      }
      dae_test_ok(&T, pior < 1e-15,
                  "PORTAO: modulo e edgelist dao a MESMA amplitude (%.3e)", pior);
      dae_test_note("   . portao modulo-vs-edgelist: diferenca %.3e (com materia e CP)", pior);
    } else {
      dae_test_ok(&T, 0, "PORTAO: uma das duas rotas falhou");
    }
  }

  /* -------------------------------- ANTI-VACUIDADE: o modulo enxerga algo? */
  {
    double hre[9], him[9];
    dae_nu_params_default(&P);

    /* delta_cp = 0 tem de dar H real — sem mergulho */
    P.delta_cp = 0.0;
    dae_nu_hamiltoniano(&P, 0.0, 0.0, hre, him);
    dae_test_ok(&T, dae_nu_precisa_mergulho(him, 3) == 0,
                "delta_cp=0 nao precisa de mergulho");

    /* dois sabores nunca precisam: fase CP nao existe em 2x2 */
    P.n_flavors = 2; P.delta_cp = -PI / 2.0;
    dae_nu_hamiltoniano(&P, 0.0, 0.0, hre, him);
    dae_test_ok(&T, dae_nu_precisa_mergulho(him, 2) == 0,
                "2 sabores nao tem fase CP, mesmo com delta_cp != 0");

    /* entrada invalida e RECUSADA, nao silenciosamente aceita */
    dae_nu_params_default(&P);
    P.energy_eV = 0.0;
    dae_test_ok(&T, dae_nu_hamiltoniano(&P, 0.0, 0.0, hre, him) == DAE_ERR_PARAM,
                "energia zero e recusada");
    dae_nu_params_default(&P);
    P.n_flavors = 4;
    dae_test_ok(&T, dae_nu_hamiltoniano(&P, 0.0, 0.0, hre, him) == DAE_ERR_PARAM,
                "4 sabores e recusado");
  }

  return dae_test_end(&T);
}
