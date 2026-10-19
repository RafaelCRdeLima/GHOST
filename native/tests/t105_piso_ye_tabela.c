/* t105_piso_ye_tabela.c — ACEITAÇÃO 105: o piso de factibilidade de
 * dae_nu_propaga tem de usar o MESMO Y_e que a propagação usa.
 *
 * O DEFEITO. Se `escala_em` avaliasse H com `M->ye`, o Y_e escalar do perfil,
 * enquanto `passo_fatia` avalia com `dae_nu_ye(M, r)`, que com `ye_table`
 * interpola a tabela, os dois só coincidiriam quando a tabela é constante e
 * igual ao escalar. Quando o escalar é menor, o piso subestimaria o raio
 * espectral, escolheria poucos subpassos, e o Chebyshev receberia um alpha
 * grande demais: DAE_ERR_TOOBIG numa fatia que o método resolve. Quando é
 * maior, o piso só desperdiçaria subpassos — e esse lado não aparece em teste
 * nenhum, o que é o motivo de o caso abaixo ser o lado que falha.
 *
 * O CASO. Perfil tabelado com rho = 1e14 g/cm3 e ye_table = 0,5, mas Y_e
 * escalar 0 — o que um spec com `ye_table` pode muito bem trazer, já que o
 * escalar deixa de ser usado na propagação. Fatia de 1 m: o alpha verdadeiro
 * é ~1e7, que exige ~2^16 subpassos; um piso com Y_e = 0 só veria o termo de
 * vácuo, estimaria alpha < 1 e não subdividiria, e o Chebyshev recusaria a
 * fatia (K acima de 1e6).
 *
 * Com densidade menor o defeito NÃO aparece: a 1e10 g/cm3 o alpha verdadeiro
 * é ~1e3, que o Chebyshev resolve num passo só, e um teste nessa densidade
 * passaria com o defeito presente.
 *
 * O PORTÃO. Com o piso correto, o perfil tabelado tem de dar o MESMO estado
 * que um perfil constante com rho e Y_e iguais aos da tabela, dentro da
 * tolerância do guarda de passo.
 *
 * ANTI-VACUIDADE. O mesmo perfil com o Y_e escalar IGUAL ao da tabela tem de
 * passar qualquer que seja o Y_e que o piso usa. Sem isso, uma falha no caso
 * principal poderia ser de qualquer outra coisa no perfil, e não do Y_e do
 * piso.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

static const double R_TAB[3]   = { 1.0, 2.0, 3.0 };          /* km */
static const double RHO_TAB[3] = { 1.0e14, 1.0e14, 1.0e14 };  /* g/cm3 */
static const double YE_TAB[3]  = { 0.5, 0.5, 0.5 };

static void params(dae_nu_params *P)
{
  dae_nu_params_default(P);
  P->n_flavors = 2;          /* real: sem mergulho, o caso mais simples */
  P->energy_eV = 20.0e6;
}

static void perfil_tabela(dae_nu_perfil *M, double ye_escalar)
{
  dae_nu_perfil_default(M);
  M->kind = DAE_NU_PERFIL_TABELA;
  M->ye = ye_escalar;
  M->r_tab = R_TAB; M->rho_tab = RHO_TAB; M->ye_tab = YE_TAB;
  M->n_tab = 3;
  M->r_start = R_TAB[0];
}

/* Propaga nu_e numa fatia de 1 m e devolve o estado e o status. */
static dae_status fatia(const dae_nu_perfil *M, double *re, double *im,
                        int32_t *subdiv)
{
  dae_nu_params P;
  const double r0 = 1.5, r1 = 1.501, tol = 1e-8;
  params(&P);
  re[0] = 1.0; im[0] = 0.0; re[1] = 0.0; im[1] = 0.0;
  return dae_nu_propaga(&P, M, r0, r1, tol, 0, re, im, subdiv);
}

int main(void)
{
  dae_test T;
  dae_nu_perfil M;
  double re_ref[2], im_ref[2], re[2], im[2];
  int32_t sub_ref = 0, sub = 0, i;
  dae_status st;

  dae_test_begin(&T, "t105 piso de factibilidade com ye_table");

  /* Referencia: perfil constante com rho e Y_e da tabela. */
  dae_nu_perfil_default(&M);
  M.kind = DAE_NU_PERFIL_CONSTANTE;
  M.rho = RHO_TAB[0];
  M.ye = YE_TAB[0];
  st = fatia(&M, re_ref, im_ref, &sub_ref);
  dae_test_ok(&T, st == DAE_OK, "referencia constante propaga (%s)",
              dae_strerror(st));
  dae_test_ok(&T, sub_ref > 1000,
              "a fatia de referencia exige subdivisao de verdade (%d subpassos)",
              (int)sub_ref);

  /* Anti-vacuidade: Y_e escalar igual ao da tabela. */
  perfil_tabela(&M, YE_TAB[0]);
  st = fatia(&M, re, im, &sub);
  dae_test_ok(&T, st == DAE_OK, "tabela com Y_e escalar consistente propaga (%s)",
              dae_strerror(st));
  if (st == DAE_OK) {
    double pior = 0.0;
    for (i = 0; i < 2; ++i) {
      const double d = hypot(re[i] - re_ref[i], im[i] - im_ref[i]);
      if (d > pior) pior = d;
    }
    dae_test_ok(&T, pior < 1e-7,
                "tabela consistente concorda com o constante (%.3e)", pior);
  }

  /* O caso do defeito: Y_e escalar 0, tabela 0,5. */
  perfil_tabela(&M, 0.0);
  st = fatia(&M, re, im, &sub);
  dae_test_ok(&T, st == DAE_OK,
              "PORTAO: tabela com Y_e escalar divergente propaga (%s)",
              dae_strerror(st));
  if (st == DAE_OK) {
    double pior = 0.0;
    for (i = 0; i < 2; ++i) {
      const double d = hypot(re[i] - re_ref[i], im[i] - im_ref[i]);
      if (d > pior) pior = d;
    }
    dae_test_ok(&T, pior < 1e-7,
                "PORTAO: o piso usa o Y_e da tabela (diferenca %.3e, %d subpassos)",
                pior, (int)sub);
  }

  return dae_test_end(&T);
}
