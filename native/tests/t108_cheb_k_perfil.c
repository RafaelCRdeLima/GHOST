/* t108_cheb_k_perfil.c — ACEITAÇÃO 108: com perfil variável, `info.k_used`
 * (que vira `#! cheb_k` no CSV) é o MÁXIMO de subpassos por fatia que o
 * guarda de dae_nu_propaga escolheu.
 *
 * O DEFEITO. No caminho de perfil variável, dae_run guarda o máximo de
 * subpassos em `R->info.k_used` dentro do laço. Se, depois dele, o campo fosse
 * sobrescrito com `ci.k_used` — que só é preenchido no caminho de H fixo —,
 * todo CSV com perfil `table` ou `power_law` sairia com `#! cheb_k` de memória
 * não inicializada: comportamento indefinido, com valores como -64 ou
 * 808466230. Os observáveis não seriam afetados; o diagnóstico sim.
 *
 * O TESTE refaz as MESMAS fatias que dae_run anda (r_start + s * t1/nt),
 * chamando dae_nu_propaga diretamente, e exige que o máximo coincida.
 *
 * ANTI-VACUIDADE. O perfil é escolhido para exigir subdivisão de verdade (o
 * máximo esperado tem de ser > 1): um teste em que todo mundo dá 1 passaria
 * com qualquer valor fixo. E o caminho de H FIXO tem de continuar reportando
 * a ordem de Chebyshev (positiva), que é o outro significado do campo.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

/* Envelope r^-3 de supernova, nu_e de 20 MeV, 2 sabores: o perfil atravessa a
   ressonancia e o guarda subdivide. */
static const char *SPEC_PERFIL =
"{\"format_version\":3,\"seed\":1,"
"\"neutrino\":{\"n_flavors\":2,\"theta12_deg\":8.6,\"dm21_sq_eV2\":2.51e-3,"
"\"energy_MeV\":20.0,\"flavor\":\"e\"},"
"\"matter\":{\"profile\":\"table\",\"tolerance\":1e-6,"
"\"r_km\":[1000.0,3000.0,10000.0,30000.0,100000.0],"
"\"rho_table\":[1.0e7,3.7e5,1.0e4,3.7e2,1.0e1],"
"\"ye_table\":[0.5,0.49,0.48,0.47,0.46]},"
"\"time\":{\"t1\":99000.0,\"nt\":40},"
"\"observables\":{\"target\":0,\"population\":false},\"realizations\":1}";

static const char *SPEC_FIXO =
"{\"format_version\":3,\"seed\":1,"
"\"neutrino\":{\"n_flavors\":2,\"theta12_deg\":33.4,\"dm21_sq_eV2\":7.42e-5,"
"\"energy_MeV\":10.0,\"flavor\":\"e\",\"rho_g_cm3\":2.8,\"Ye\":0.5},"
"\"time\":{\"t1\":1000.0,\"nt\":10},"
"\"observables\":{\"target\":0,\"population\":false},\"realizations\":1}";

int main(void)
{
  dae_test T;
  dae_spec S;
  dae_graph G;
  dae_series R;
  dae_error err;
  dae_status st;

  dae_test_begin(&T, "t108 cheb_k com perfil variavel e o maximo de subpassos");

  st = dae_spec_parse(&S, SPEC_PERFIL, &err);
  if (st != DAE_OK) { dae_test_ok(&T, 0, "spec: %s", err.msg); return dae_test_end(&T); }

  {
    /* As fatias de dae_run, refeitas aqui. 2 sabores sem fase: sem mergulho. */
    const double passo = S.t1 / (double)S.nt;
    double re[2] = { 1.0, 0.0 }, im[2] = { 0.0, 0.0 };
    int32_t s, sub, maximo = 0;
    for (s = 0; s < S.nt; ++s) {
      const double ra = S.nu_perfil.r_start + (double)s * passo;
      const double rb = S.nu_perfil.r_start + (double)(s + 1) * passo;
      sub = 0;
      st = dae_nu_propaga(&S.nu, &S.nu_perfil, ra, rb, S.nu_tol, 0, re, im, &sub);
      if (st != DAE_OK) break;
      if (sub > maximo) maximo = sub;
    }
    dae_test_ok(&T, st == DAE_OK, "propagacao direta das fatias (%s)", dae_strerror(st));
    dae_test_ok(&T, maximo > 1,
                "o perfil exige subdivisao de verdade (maximo %d)", (int)maximo);

    st = dae_run(&S, S.seed, &G, NULL, &R, NULL, NULL);
    dae_test_ok(&T, st == DAE_OK, "dae_run com perfil (%s)", dae_strerror(st));
    if (st == DAE_OK) {
      dae_test_ok(&T, R.info.k_used == maximo,
                  "PORTAO: info.k_used = maximo de subpassos (%d contra %d)",
                  (int)R.info.k_used, (int)maximo);
      dae_series_free(&R); dae_graph_free(&G);
    }
  }
  dae_spec_free(&S);

  /* Companheira: H fixo continua reportando a ordem de Chebyshev. */
  st = dae_spec_parse(&S, SPEC_FIXO, &err);
  if (st != DAE_OK) { dae_test_ok(&T, 0, "spec fixo: %s", err.msg); return dae_test_end(&T); }
  st = dae_run(&S, S.seed, &G, NULL, &R, NULL, NULL);
  dae_test_ok(&T, st == DAE_OK, "dae_run com H fixo (%s)", dae_strerror(st));
  if (st == DAE_OK) {
    dae_test_ok(&T, R.info.k_used > 0 && R.info.k_used < 100000,
                "H fixo: k_used e a ordem de Chebyshev (%d)", (int)R.info.k_used);
    dae_series_free(&R); dae_graph_free(&G);
  }
  dae_spec_free(&S);

  return dae_test_end(&T);
}
