/* t104_observaveis_nu.c — ACEITAÇÃO 104: as três famílias de observável, e o
 * fato de que elas são DIFERENTES.
 *
 * O CSV do GHOST emite três blocos de colunas, e a distinção entre eles é a
 * substância deste arquivo:
 *
 *   p_e, p_mu, p_tau      sabor COERENTE, no ponto. É o que o propagador dá.
 *   m_1, m_2, m_3         base de MASSA. É a métrica de convergência de um
 *                         perfil — P de sabor oscila e engana (ver t101).
 *   inc_e, inc_mu, inc_tau  INCOERENTE. Para fonte astrofísica é ESTE o
 *                         observável: os pacotes dos autoestados de massa se
 *                         separam em ~32 km, muito antes de sair da galáxia,
 *                         então a fase que define o P coerente já não existe
 *                         quando o neutrino chega ao detector.
 *
 * ANTI-VACUIDADE, e é a razão de o arquivo existir: as três famílias somam 1
 * cada. Um erro que copiasse uma coluna sobre a outra passaria em todo teste de
 * normalização, e o CSV sairia com nomes certos e números errados. O teste
 * exige que elas sejam MUTUAMENTE DIFERENTES.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static const char *SPEC_DUNE =
"{\"format_version\":3,\"seed\":1,"
"\"neutrino\":{\"n_flavors\":3,\"theta12_deg\":33.4,\"theta13_deg\":8.6,"
"\"theta23_deg\":49.0,\"delta_cp_deg\":-90.0,\"dm21_sq_eV2\":7.42e-5,"
"\"dm31_sq_eV2\":2.51e-3,\"energy_MeV\":2500.0,\"flavor\":\"mu\","
"\"rho_g_cm3\":2.8,\"Ye\":0.5},"
"\"time\":{\"t1\":1300.0,\"nt\":100},"
"\"observables\":{\"target\":0,\"population\":false},\"realizations\":1}";

/* Perfil de supernova: a conversao adiabatica manda `inc_e` para sin^2(th13). */
static const char *SPEC_SN =
"{\"format_version\":3,\"seed\":1,"
"\"neutrino\":{\"n_flavors\":3,\"theta12_deg\":33.4,\"theta13_deg\":8.6,"
"\"theta23_deg\":49.0,\"delta_cp_deg\":-90.0,\"dm21_sq_eV2\":7.42e-5,"
"\"dm31_sq_eV2\":2.51e-3,\"energy_MeV\":20.0,\"flavor\":\"e\"},"
"\"matter\":{\"profile\":\"power_law\",\"Ye\":0.5,\"rho_ref_g_cm3\":1.0e10,"
"\"r_ref_km\":100.0,\"exponent\":3.0,\"r_start_km\":5000.0,\"tolerance\":1e-6},"
"\"time\":{\"t1\":995000.0,\"nt\":300},"
"\"observables\":{\"target\":0,\"population\":false},\"realizations\":1}";

static int roda(dae_test *T, const char *texto, dae_spec *S,
                dae_graph *G, dae_series *R)
{
  dae_error err;
  dae_status st = dae_spec_parse(S, texto, &err);
  if (st != DAE_OK) { dae_test_ok(T, 0, "spec: %s", err.msg); return 1; }
  st = dae_run(S, S->seed, G, NULL, R, NULL, NULL);
  if (st != DAE_OK) { dae_test_ok(T, 0, "run: %s", dae_strerror(st));
                      dae_spec_free(S); return 1; }
  return 0;
}

int main(void)
{
  dae_test T;
  dae_spec S;
  dae_graph G;
  dae_series R;

  dae_test_begin(&T, "t104 sabor, massa e incoerente — e a diferenca entre eles");

  /* -------------------------------------------------- as tres normalizacoes */
  if (roda(&T, SPEC_DUNE, &S, &G, &R) == 0) {
    double pior_s = 0.0, pior_m = 0.0, pior_i = 0.0;
    double dif_sm = 0.0, dif_si = 0.0, dif_mi = 0.0;
    int32_t s, a;
    const int32_t nf = R.nu_nf;

    dae_test_ok(&T, R.nu_sabor && R.nu_massa && R.nu_inc && nf == 3,
                "as tres series existem, com 3 sabores");

    for (s = 0; s < R.nt; ++s) {
      const size_t o = (size_t)s * (size_t)nf;
      double ss = 0.0, sm = 0.0, si = 0.0;
      for (a = 0; a < nf; ++a) {
        ss += R.nu_sabor[o + (size_t)a];
        sm += R.nu_massa[o + (size_t)a];
        si += R.nu_inc[o + (size_t)a];
      }
      if (fabs(ss - 1.0) > pior_s) pior_s = fabs(ss - 1.0);
      if (fabs(sm - 1.0) > pior_m) pior_m = fabs(sm - 1.0);
      if (fabs(si - 1.0) > pior_i) pior_i = fabs(si - 1.0);
      for (a = 0; a < nf; ++a) {
        const double vs = R.nu_sabor[o + (size_t)a];
        const double vm = R.nu_massa[o + (size_t)a];
        const double vi = R.nu_inc[o + (size_t)a];
        if (fabs(vs - vm) > dif_sm) dif_sm = fabs(vs - vm);
        if (fabs(vs - vi) > dif_si) dif_si = fabs(vs - vi);
        if (fabs(vm - vi) > dif_mi) dif_mi = fabs(vm - vi);
      }
    }
    dae_test_ok(&T, pior_s < 1e-12, "sabor soma 1 (pior %.3e)", pior_s);
    dae_test_ok(&T, pior_m < 1e-12, "massa soma 1 (pior %.3e)", pior_m);
    dae_test_ok(&T, pior_i < 1e-12, "incoerente soma 1 (pior %.3e)", pior_i);

    /* ANTI-VACUIDADE: as tres tem de ser MUTUAMENTE diferentes. */
    dae_test_ok(&T, dif_sm > 1e-3, "sabor != massa (max %.6f)", dif_sm);
    dae_test_ok(&T, dif_si > 1e-3, "sabor != incoerente (max %.6f)", dif_si);
    dae_test_ok(&T, dif_mi > 1e-3, "massa != incoerente (max %.6f)", dif_mi);
    dae_test_note("   . as tres somam 1, e diferem entre si por ate "
                  "%.4f / %.4f / %.4f", dif_sm, dif_si, dif_mi);

    /* p_target tem de ser a coluna de sabor do alvo (target = 0 = nu_e). */
    {
      double pior = 0.0;
      for (s = 0; s < R.nt; ++s) {
        const double pa = R.scal[(size_t)s * DAE_S_NCOL + DAE_S_P_ALVO];
        const double pe = R.nu_sabor[(size_t)s * (size_t)nf];
        if (fabs(pa - pe) > pior) pior = fabs(pa - pe);
      }
      dae_test_ok(&T, pior < 1e-15, "p_target == p_e quando target=0 (%.3e)", pior);
    }
    dae_test_note("   . DUNE 1300 km: p_e=%.9f  m_3=%.9f  inc_e=%.9f",
                  R.nu_sabor[(size_t)(R.nt - 1) * 3],
                  R.nu_massa[(size_t)(R.nt - 1) * 3 + 2],
                  R.nu_inc[(size_t)(R.nt - 1) * 3]);
    dae_series_free(&R); dae_graph_free(&G); dae_spec_free(&S);
  }

  /* ------------- o incoerente e o que a fisica de supernova preve --------- */
  if (roda(&T, SPEC_SN, &S, &G, &R) == 0) {
    const int32_t u = R.nt - 1;
    const double sin2t13 = sin(S.nu.theta13) * sin(S.nu.theta13);
    const double m3 = R.nu_massa[(size_t)u * 3 + 2];
    const double inc_e = R.nu_inc[(size_t)u * 3];
    const double p_e = R.nu_sabor[(size_t)u * 3];

    dae_test_ok(&T, m3 > 1.0 - 1e-4,
                "a conversao foi adiabatica: m_3 = %.9f", m3);
    dae_test_near(&T, inc_e, sin2t13, 1e-3,
                  "inc_e contra sin^2(theta13) — a previsao de Dighe-Smirnov");
    dae_test_note("   . supernova: m_3=%.9f  inc_e=%.9f  sin^2(th13)=%.9f",
                  m3, inc_e, sin2t13);
    dae_test_note("   . e p_e (coerente) = %.9f — mesma corrida, outro numero:", p_e);
    dae_test_note("     e por isso que reportar o coerente seria reportar uma fase morta");
    dae_series_free(&R); dae_graph_free(&G); dae_spec_free(&S);
  }

  return dae_test_end(&T);
}
