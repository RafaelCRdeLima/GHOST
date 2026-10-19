/* t113_numag_zeros.c — ACEITAÇÃO 113: o localizador de zeros de Delta_ef do
 * integrador por janelas.
 *
 * O diagnóstico de cruzamentos alimenta larguras, gradientes e Gamma ao longo
 * das trajetórias; ele tem de distinguir um cruzamento isolado de um contato
 * sem troca de sinal e de um intervalo com Delta_ef identicamente nulo.
 *
 *  1. RAIZ ASCENDENTE NO INÍCIO DA CÉLULA: Y_e de 1/3 a 1/3 + 1e-7 em
 *     [0, 1e-4] km. O zero está em s = 0, não no fim da célula
 *     (s = 1,25e-5). Companheira: a rampa descendente também acha s = 0.
 *  2. ZERO NUM PONTO DE AMOSTRAGEM INTERIOR: rampa com a raiz exatamente no nó
 *     do meio da tabela — um cruzamento, no nó, contado UMA vez.
 *  3. CONTATO: Delta_ef = 0 num nó, com o mesmo sinal dos dois lados (|Y_e -
 *     1/3| em V): nenhum cruzamento, um contato.
 *  4. INTERVALO NULO: Y_e = 1/3 em toda a tabela (8 células): nenhum
 *     cruzamento, UM intervalo nulo com o comprimento inteiro, e não oito
 *     cruzamentos, um por célula.
 *  5. RAMPA COMUM (companheira): raiz fora dos pontos de amostragem, um
 *     cruzamento localizado por bisseção a 1e-12 do exato, zero contatos e
 *     zero intervalos.
 *  6. FATIAS: o CLI integra uma chamada por fatia de saída. Dividir a
 *     trajetória EXATAMENTE no zero não pode mudar a classificação: o contato
 *     em V continua contato (e não vira cruzamento), o intervalo nulo continua
 *     um (e não dois), o cruzamento no nó continua um; e P_R é a mesma, bit a
 *     bit, da chamada única.
 */
#include "harness.h"

#include <math.h>
#include <string.h>

#define RHO 1.0e14

/* Integra [s[0], s[n-1]] em nfat chamadas, cortando nos pontos cortes[]. */
static dae_status roda_fatias(const double *s, const double *ye, int32_t n,
                              const double *cortes, int nfat,
                              dae_numag_janelas_info *f, double *pr)
{
  static double rho[16], b1[16], b2[16];
  dae_numag_trajetoria t;
  dae_numag_params p;
  dae_numag_janelas_cfg c;
  double re[4], im[4];
  int32_t i;
  for (i = 0; i < n; ++i) { rho[i] = RHO; b1[i] = 1e15; b2[i] = 0.0; }
  memset(&t, 0, sizeof t);
  t.n = n; t.s_km = s; t.rho = rho; t.ye = ye; t.b1 = b1; t.b2 = b2;
  dae_numag_params_default(&p);
  p.mu_re[0] = 1e-12;
  dae_numag_janelas_default(&c);
  memset(f, 0, sizeof *f);
  dae_nu_inicial(0, 2, 1, re, im);
  {
    dae_status st = DAE_OK;
    int k;
    for (k = 0; k < nfat && st == DAE_OK; ++k) {
      const double a = (k == 0) ? s[0] : cortes[k - 1];
      const double b = (k == nfat - 1) ? s[n - 1] : cortes[k];
      st = dae_numag_propaga_janelas(&p, &t, a, b, &c, re, im, f);
    }
    if (pr) *pr = re[1] * re[1] + im[1] * im[1] + re[3] * re[3] + im[3] * im[3];
    return st;
  }
}

static dae_status roda(const double *s, const double *ye, int32_t n,
                       dae_numag_janelas_info *f)
{
  return roda_fatias(s, ye, n, NULL, 1, f, NULL);
}

int main(void)
{
  dae_test T;
  dae_numag_janelas_info f;
  dae_status st;

  dae_test_begin(&T, "t113 zeros de Delta_ef: inicio de celula, ponto de amostragem, contato, intervalo nulo");

  /* -------------------------------------------- 1. início da célula -- */
  {
    const double s[2] = { 0.0, 1e-4 };
    double ye[2] = { 1.0 / 3.0, 1.0 / 3.0 + 1e-7 };
    st = roda(s, ye, 2, &f);
    dae_test_ok(&T, st == DAE_OK && f.n_cruzamentos == 1 && f.cruzamento_s_km == 0.0,
                "1. raiz ascendente em s = 0: %d cruzamento(s) em s = %.3e",
                (int)f.n_cruzamentos, f.cruzamento_s_km);
    ye[1] = 1.0 / 3.0 - 1e-7;
    st = roda(s, ye, 2, &f);
    dae_test_ok(&T, st == DAE_OK && f.n_cruzamentos == 1 && f.cruzamento_s_km == 0.0,
                "1. companheira: descendente, %d cruzamento(s) em s = %.3e",
                (int)f.n_cruzamentos, f.cruzamento_s_km);
  }

  /* ------------------------------ 2. zero num ponto de amostragem -- */
  {
    const double s[3] = { 0.0, 1e-4, 2e-4 };
    const double ye[3] = { 1.0 / 3.0 - 1e-7, 1.0 / 3.0, 1.0 / 3.0 + 1e-7 };
    st = roda(s, ye, 3, &f);
    dae_test_ok(&T, st == DAE_OK && f.n_cruzamentos == 1 && f.cruzamento_s_km == 1e-4 &&
                f.n_contatos == 0,
                "2. raiz no no interior: %d cruzamento(s) em s = %.6e, %d contato(s)",
                (int)f.n_cruzamentos, f.cruzamento_s_km, (int)f.n_contatos);
  }

  /* ------------------------------------------------- 3. contato -- */
  {
    const double s[3] = { 0.0, 1e-4, 2e-4 };
    const double ye[3] = { 1.0 / 3.0 + 1e-7, 1.0 / 3.0, 1.0 / 3.0 + 1e-7 };
    st = roda(s, ye, 3, &f);
    dae_test_ok(&T, st == DAE_OK && f.n_cruzamentos == 0 && f.n_contatos == 1,
                "3. contato sem troca de sinal: %d cruzamento(s), %d contato(s)",
                (int)f.n_cruzamentos, (int)f.n_contatos);
  }

  /* -------------------------------------------- 4. intervalo nulo -- */
  {
    const double s[2] = { 0.0, 1e-4 };
    const double ye[2] = { 1.0 / 3.0, 1.0 / 3.0 };
    st = roda(s, ye, 2, &f);
    dae_test_ok(&T, st == DAE_OK && f.n_cruzamentos == 0 && f.n_intervalos_nulos == 1 &&
                fabs(f.comprimento_nulo_km - 1e-4) < 1e-18,
                "4. Delta_ef = 0 em toda a tabela: %d cruzamento(s), %d intervalo(s) nulo(s) "
                "de %.3e km",
                (int)f.n_cruzamentos, (int)f.n_intervalos_nulos, f.comprimento_nulo_km);
  }

  /* ---------------------------------------------- 5. rampa comum -- */
  {
    const double s[2] = { 0.0, 1e-4 };
    /* zero em 1/3 do intervalo: fora dos pontos de amostragem (h/8) */
    const double ye[2] = { 1.0 / 3.0 - 1e-7, 1.0 / 3.0 + 2e-7 };
    const double exato = 1e-4 / 3.0;
    st = roda(s, ye, 2, &f);
    dae_test_ok(&T, st == DAE_OK && f.n_cruzamentos == 1 &&
                fabs(f.cruzamento_s_km - exato) < 1e-12 &&
                f.n_contatos == 0 && f.n_intervalos_nulos == 0,
                "5. companheira: rampa comum, cruzamento em %.12e (exato %.12e)",
                f.cruzamento_s_km, exato);
  }

  /* ---------------------------------------------------- 6. fatias -- */
  {
    const double s[3] = { 0.0, 1e-4, 2e-4 };
    const double corte[1] = { 1e-4 };
    const double v[3] = { 1.0 / 3.0 + 1e-7, 1.0 / 3.0, 1.0 / 3.0 + 1e-7 };
    const double nulo[3] = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
    const double rampa[3] = { 1.0 / 3.0 - 1e-7, 1.0 / 3.0, 1.0 / 3.0 + 1e-7 };
    dae_numag_janelas_info g;
    dae_status s2;
    double p1, p2;

    st = roda_fatias(s, v, 3, NULL, 1, &f, &p1);
    s2 = roda_fatias(s, v, 3, corte, 2, &g, &p2);
    dae_test_ok(&T, st == DAE_OK && s2 == DAE_OK && g.n_cruzamentos == 0 &&
                g.n_contatos == 1 && f.n_contatos == 1 && p1 == p2,
                "6. contato no corte: uma chamada %d/%d, duas %d/%d (cruz./contatos); "
                "P_R %.17g / %.17g", (int)f.n_cruzamentos, (int)f.n_contatos,
                (int)g.n_cruzamentos, (int)g.n_contatos, p1, p2);

    st = roda_fatias(s, nulo, 3, NULL, 1, &f, &p1);
    s2 = roda_fatias(s, nulo, 3, corte, 2, &g, &p2);
    dae_test_ok(&T, st == DAE_OK && s2 == DAE_OK && g.n_intervalos_nulos == 1 &&
                g.n_cruzamentos == 0 && fabs(g.comprimento_nulo_km - 2e-4) < 1e-18 &&
                p1 == p2,
                "6. intervalo nulo cortado: %d intervalo(s) de %.3e km em duas chamadas", (int)g.n_intervalos_nulos, g.comprimento_nulo_km);

    st = roda_fatias(s, rampa, 3, NULL, 1, &f, &p1);
    s2 = roda_fatias(s, rampa, 3, corte, 2, &g, &p2);
    dae_test_ok(&T, st == DAE_OK && s2 == DAE_OK && g.n_cruzamentos == 1 &&
                g.cruzamento_s_km == 1e-4 && g.n_contatos == 0 && p1 == p2,
                "6. cruzamento no corte: %d cruzamento(s) em s = %.6e",
                (int)g.n_cruzamentos, g.cruzamento_s_km);

    /* corte fora de nó: 4 fatias iguais (o caso do CLI com nt = 4). Aqui os
       cortes mudam a partição dos passos do integrador, e P_R só concorda
       dentro da tolerância de força bruta (tol = 1e-10), não bit a bit. */
    {
      const double q[3] = { 5e-5, 1e-4, 1.5e-4 };
      s2 = roda_fatias(s, v, 3, q, 4, &g, &p2);
      st = roda_fatias(s, v, 3, NULL, 1, &f, &p1);
      dae_test_ok(&T, st == DAE_OK && s2 == DAE_OK && g.n_cruzamentos == 0 &&
                  g.n_contatos == 1 && fabs(p1 - p2) < 1e-10,
                  "6. contato com nt = 4: %d/%d, |dP_R| = %.2e",
                  (int)g.n_cruzamentos, (int)g.n_contatos, fabs(p1 - p2));
    }
  }

  return dae_test_end(&T);
}
