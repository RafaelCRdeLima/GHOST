/* t112_numag_fundo_neutrinos.c — ACEITAÇÃO 112: fundo de nu_e aprisionados no
 * desajuste do setor magnético.
 *
 *     Delta = sqrt(2) G_F n_b [(3 Y_e - 1)/2 + 2 Y_nue]   (Nötzold e Raffelt)
 *
 *  1. O VALOR: dae_numag_desajuste contra a expressão, e o zero deslocado para
 *     Y_e + (4/3) Y_nue = 1/3. Companheira: Y_nue = 0 dá o valor sem fundo
 *     bit a bit, e o hamiltoniano recusa Y_nue NaN ou |Y_nue| > 1.
 *  2. A TABELA: Y_nue linear entre nós, exato nos nós; ynu = NULL vale 0, e uma
 *     tabela de zeros propaga BIT A BIT igual a NULL.
 *  3. CRUZAMENTO PELO FUNDO: Y_e fixo em 0,3 (longe de 1/3) e rampa de Y_nue
 *     através de 0,025. Força bruta contra Landau-Zener, e contra a tabela
 *     EQUIVALENTE Y_e' = Y_e + (4/3) Y_nue sem fundo (mesmo Delta nos nós e
 *     linear entre eles: tem de coincidir a arredondamento). Companheira: a
 *     mesma tabela sem o fundo não converte (P_R ~ (Omega/Delta)^2).
 *  4. O INTEGRADOR POR JANELAS lê a DERIVADA de Delta com o fundo: na rampa
 *     longa (janela K = 50 dentro de um domínio de K = 500), estado final,
 *     cota adiabática, eps fora da janela e posição do cruzamento coincidem
 *     com os da tabela equivalente. Um dY_nue/ds esquecido em jan_ponto zera
 *     theta~' e muda a cota; o teste falharia por ordens de grandeza.
 *  5. SPEC: trajectory.ynu_table é lido, vai ao canônico (e só quando existe),
 *     tamanho e |Y_nu| <= 1 conferidos; ponta a ponta, Y_e = 0,3 com
 *     Y_nue = 0,025 uniforme está na ressonância: P_R = sin^2(Omega s).
 */
#include "harness.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define TK DAE_NU_T_POR_KM
#define RHO 1.0e14
#define VR (DAE_NU_V_POR_RHO * RHO)
#define YE0 0.3
#define YNU_RES 0.025             /* (1/3 - 0,3) * 3/4 */

static dae_numag_trajetoria traj(const double *s, const double *rho, const double *ye,
                                 const double *ynu, const double *b1, const double *b2,
                                 int32_t n)
{
  dae_numag_trajetoria T;
  T.s_km = s; T.rho = rho; T.ye = ye; T.ynu = ynu; T.b1 = b1; T.b2 = b2;
  T.razao_energia = NULL; T.n = n;
  return T;
}

/* P_R partindo de nu_L, força bruta em nf fatias */
static double pr_bruta(const dae_numag_params *P, dae_numag_trajetoria *T, double s0,
                       double s1, int32_t nf, double re[4], double im[4])
{
  double pr[2], pi[2];
  int32_t k;
  dae_nu_inicial(0, 2, 1, re, im);
  for (k = 0; k < nf; ++k)
    if (dae_numag_propaga_tabela(P, T, s0 + (s1 - s0) * k / nf, s0 + (s1 - s0) * (k + 1) / nf,
                                 1e-11, re, im, NULL) != DAE_OK)
      return -1.0;
  dae_nu_projeta(re, im, 2, 1, pr, pi);
  return pr[1] * pr[1] + pi[1] * pi[1];
}

static double pr_jan(const dae_numag_params *P, const dae_numag_trajetoria *T, double s0,
                     double s1, const dae_numag_janelas_cfg *c, dae_numag_janelas_info *inf,
                     double re[4], double im[4])
{
  double pr[2], pi[2];
  memset(inf, 0, sizeof(*inf));
  dae_nu_inicial(0, 2, 1, re, im);
  {
    const dae_status st = dae_numag_propaga_janelas(P, T, s0, s1, c, re, im, inf);
    if (st != DAE_OK) { printf("   ! propaga_janelas: %s\n", dae_strerror(st)); return -1.0; }
  }
  dae_nu_projeta(re, im, 2, 1, pr, pi);
  return pr[1] * pr[1] + pi[1] * pi[1];
}

static double dmax(const double *a, const double *b, int n)
{
  double m = 0.0;
  int k;
  for (k = 0; k < n; ++k) if (fabs(a[k] - b[k]) > m) m = fabs(a[k] - b[k]);
  return m;
}

static char *spec_json(double L, int32_t nt, const char *ynu)
{
  char *b = (char *)malloc(1024);
  snprintf(b, 1024,
           "{\"format_version\":4,\"seed\":1,"
           "\"magnetic\":{\"nature\":\"dirac\",\"mu_nu_bohr\":1e-12,\"energy_MeV\":10.0,"
           "\"initial_helicity\":\"left\",\"tolerance\":1e-11},"
           "\"trajectory\":{\"transverse_basis\":\"parallel_transported\",\"s_km\":[0,%.17g],"
           "\"rho_table\":[1e14,1e14],\"ye_table\":[0.3,0.3]%s,"
           "\"b1_gauss\":[1e15,1e15],\"b2_gauss\":[0,0]},"
           "\"time\":{\"t1\":%.17g,\"nt\":%d},\"realizations\":1}",
           L, ynu ? ynu : "", L, (int)nt);
  return b;
}

int main(void)
{
  dae_test T;
  dae_numag_params P;
  const double om = 1e-12 * DAE_NUMAG_MU_B_EV_POR_T * 1e15 * DAE_NUMAG_T_POR_GAUSS;

  dae_test_begin(&T, "t112 fundo de nu_e no desajuste: valor, tabela, cruzamento, janelas, spec");
  dae_numag_params_default(&P);
  P.mu_re[0] = 1e-12;

  /* ------------------------------------------------------- 1. o valor -- */
  {
    dae_numag_meio m;
    double hre[4], him[4], d, ref;
    memset(&m, 0, sizeof m);
    m.rho = RHO; m.ye = 0.25; m.ynu = 0.04;
    d = dae_numag_desajuste(&m);
    ref = VR * ((3.0 * 0.25 - 1.0) / 2.0 + 2.0 * 0.04);
    dae_test_ok(&T, fabs(d - ref) <= 1e-15 * fabs(ref),
                "1. Delta = V rho [(3Ye-1)/2 + 2Ynu] (%.12e eV, ref %.12e)", d, ref);
    m.ye = YE0; m.ynu = YNU_RES;
    d = dae_numag_desajuste(&m);
    dae_test_ok(&T, fabs(d) < 1e-15 * VR,
                "1. zero em Y_e + (4/3) Y_nue = 1/3 (Ye = 0,3, Ynu = 0,025: |Delta| = %.1e eV)", fabs(d));
    m.ynu = 0.0;
    d = dae_numag_desajuste(&m);
    dae_test_ok(&T, d == DAE_NU_V_POR_RHO * RHO * (3.0 * YE0 - 1.0) * 0.5,
                "1. companheira: Y_nue = 0 da a expressao sem fundo bit a bit");
    m.ynu = NAN;
    dae_test_ok(&T, dae_numag_hamiltoniano(&P, &m, hre, him) == DAE_ERR_PARAM, "1. Y_nue NaN recusado");
    m.ynu = 1.5;
    dae_test_ok(&T, dae_numag_hamiltoniano(&P, &m, hre, him) == DAE_ERR_PARAM, "1. |Y_nue| > 1 recusado");
    m.ynu = -0.2;
    dae_test_ok(&T, dae_numag_hamiltoniano(&P, &m, hre, him) == DAE_OK,
                "1. Y_nue negativo (excesso de antineutrinos) aceito");
  }

  /* ------------------------------------------------------ 2. a tabela -- */
  {
    /* microquilômetros: com Delta ~ 0,2 eV, 1 km já são ~1e9 rad de fase, e a
       igualdade bit a bit não precisa disso */
    static const double s[3] = { 0.0, 1.0e-6, 3.0e-6 }, rho[3] = { RHO, RHO, RHO };
    static const double ye[3] = { 0.3, 0.3, 0.3 }, ynu[3] = { 0.01, 0.03, -0.01 };
    static const double z[3] = { 0.0, 0.0, 0.0 }, b1[3] = { 1e15, 1e15, 1e15 };
    dae_numag_trajetoria A = traj(s, rho, ye, ynu, b1, z, 3);
    dae_numag_trajetoria N = traj(s, rho, ye, NULL, b1, z, 3);
    dae_numag_trajetoria Z = traj(s, rho, ye, z, b1, z, 3);
    dae_numag_meio m;
    double re1[4], im1[4], re2[4], im2[4];
    dae_numag_meio_tabela(&A, 1.0e-6, &m);
    dae_test_ok(&T, m.ynu == 0.03, "2. Y_nue exato no no (%.17g)", m.ynu);
    dae_numag_meio_tabela(&A, 2.0e-6, &m);
    dae_test_near(&T, m.ynu, 0.01, 1e-16, "2. Y_nue linear entre nos");
    m.ynu = 7.0;
    dae_numag_meio_tabela(&N, 2.0e-6, &m);
    dae_test_ok(&T, m.ynu == 0.0, "2. ynu = NULL vale 0");
    pr_bruta(&P, &N, 0.0, 3.0e-6, 4, re1, im1);
    pr_bruta(&P, &Z, 0.0, 3.0e-6, 4, re2, im2);
    dae_test_ok(&T, memcmp(re1, re2, sizeof re1) == 0 && memcmp(im1, im2, sizeof im1) == 0,
                "2. tabela de zeros propaga bit a bit igual a NULL");
  }

  /* ------------------------------------- 3. cruzamento pelo fundo, LZ -- */
  {
    const double gamma = 1.0, K = 200.0;
    const double dD = 2.0 * PI * om * om * TK / gamma;        /* eV/km */
    const double meia = K * 2.0 * om / dD;
    double s[2] = { 0.0, 2.0 * meia }, rho[2] = { RHO, RHO }, ye[2] = { YE0, YE0 };
    double ynu[2], yeq[2], z[2] = { 0.0, 0.0 }, b1[2] = { 1e15, 1e15 };
    double re1[4], im1[4], re2[4], im2[4], re3[4], im3[4], p1, p2, p3;
    dae_numag_trajetoria A, E, S0;
    const double lz = 1.0 - exp(-gamma);
    int k;
    for (k = 0; k < 2; ++k) {
      ynu[k] = YNU_RES + dD * (s[k] - meia) / (2.0 * VR);
      yeq[k] = YE0 + 4.0 / 3.0 * ynu[k];
    }
    A = traj(s, rho, ye, ynu, b1, z, 2);
    E = traj(s, rho, yeq, NULL, b1, z, 2);
    S0 = traj(s, rho, ye, NULL, b1, z, 2);
    p1 = pr_bruta(&P, &A, 0.0, 2.0 * meia, 64, re1, im1);
    p2 = pr_bruta(&P, &E, 0.0, 2.0 * meia, 64, re2, im2);
    p3 = pr_bruta(&P, &S0, 0.0, 2.0 * meia, 64, re3, im3);
    dae_test_ok(&T, fabs(p1 - lz) < 2.0 / K,
                "3. rampa de Y_nue: P_R = %.5f contra LZ %.5f (dif %.1e, resolucao 2/K = %.0e)",
                p1, lz, fabs(p1 - lz), 2.0 / K);
    dae_test_ok(&T, dmax(re1, re2, 4) < 1e-9 && dmax(im1, im2, 4) < 1e-9,
                "3. estado final igual ao da tabela equivalente Y_e' = Y_e + (4/3) Y_nue "
                "(P_R %.12f x %.12f, dif de amplitude %.1e)", p1, p2,
                fmax(dmax(re1, re2, 4), dmax(im1, im2, 4)));
    dae_test_ok(&T, p3 < 1e-8,
                "3. companheira: sem o fundo, Y_e = 0,3 nao cruza (P_R = %.1e)", p3);
  }

  /* -------------------------------- 4. janelas: a derivada com o fundo -- */
  {
    enum { NN = 401 };
    const double gamma = 1.0, Kd = 500.0;
    const double dD = 2.0 * PI * om * om * TK / gamma;
    const double meia = Kd * 2.0 * om / dD;
    static double s[NN], rho[NN], ye[NN], ynu[NN], yeq[NN], z[NN], b1[NN];
    double re1[4], im1[4], re2[4], im2[4], reb[4], imb[4], p1, p2, pb;
    dae_numag_trajetoria A, E;
    dae_numag_janelas_cfg c;
    dae_numag_janelas_info ia, ie;
    int k;
    for (k = 0; k < NN; ++k) {
      s[k] = 2.0 * meia * k / (NN - 1);
      rho[k] = RHO; ye[k] = YE0; z[k] = 0.0; b1[k] = 1e15;
      ynu[k] = YNU_RES + dD * (s[k] - meia) / (2.0 * VR);
      yeq[k] = YE0 + 4.0 / 3.0 * ynu[k];
    }
    A = traj(s, rho, ye, ynu, b1, z, NN);
    E = traj(s, rho, yeq, NULL, b1, z, NN);
    dae_numag_janelas_default(&c);
    c.K = 50.0;
    p1 = pr_jan(&P, &A, 0.0, 2.0 * meia, &c, &ia, re1, im1);
    p2 = pr_jan(&P, &E, 0.0, 2.0 * meia, &c, &ie, re2, im2);
    pb = pr_bruta(&P, &A, 0.0, 2.0 * meia, 64, reb, imb);
    dae_test_ok(&T, ia.n_janelas == 1 && ia.n_trechos >= 2 && ia.cota_adiabatica > 0.0,
                "4. ha trechos adiabaticos fora da janela (%d janela, %d trechos, cota %.1e)",
                (int)ia.n_janelas, (int)ia.n_trechos, ia.cota_adiabatica);
    /* P a arredondamento; as amplitudes, até a tolerância da quadratura
       adaptativa da fase (tol_fase), que pode refinar diferente */
    dae_test_ok(&T, fabs(p1 - p2) < 1e-9 &&
                fmax(dmax(re1, re2, 4), dmax(im1, im2, 4)) < 10.0 * c.tol_fase,
                "4. janelas: estado final igual ao da tabela equivalente (dif P %.1e, amplitude %.1e, "
                "tol_fase %.0e)", fabs(p1 - p2), fmax(dmax(re1, re2, 4), dmax(im1, im2, 4)), c.tol_fase);
    dae_test_ok(&T, fabs(ia.cota_adiabatica - ie.cota_adiabatica) <= 1e-9 * ie.cota_adiabatica &&
                fabs(ia.eps_fora_max - ie.eps_fora_max) <= 1e-9 * ie.eps_fora_max,
                "4. cota (%.6e x %.6e) e eps fora (%.6e x %.6e) iguais: a derivada inclui dY_nue/ds",
                ia.cota_adiabatica, ie.cota_adiabatica, ia.eps_fora_max, ie.eps_fora_max);
    /* o cruzamento é dado no meio da célula de classificação que o contém:
       resolução h/sub */
    {
      const double cel = 2.0 * meia / (NN - 1) / c.sub;
      dae_test_ok(&T, fabs(ia.cruzamento_s_km - meia) <= cel &&
                  fabs(ia.cruzamento_s_km - ie.cruzamento_s_km) <= cel,
                  "4. cruzamento em s = %.8f km (exato %.8f, celula %.1e km), o mesmo da equivalente",
                  ia.cruzamento_s_km, meia, cel);
    }
    dae_test_ok(&T, fabs(p1 - pb) < 1e-4,
                "4. janelas contra forca bruta na rampa de Y_nue: %.6f x %.6f (dif %.1e)",
                p1, pb, fabs(p1 - pb));
  }

  /* ---------------------------------------------------------- 5. spec -- */
  {
    const double L = PI / om / TK;
    dae_spec S, S2;
    dae_error err;
    char *js;
    char ynu[128];
    snprintf(ynu, sizeof ynu, ",\"ynu_table\":[%.17g,%.17g]", YNU_RES, YNU_RES);

    js = spec_json(L, 16, ynu);
    if (dae_spec_parse(&S, js, &err) != DAE_OK) {
      dae_test_ok(&T, 0, "5. spec com ynu_table: %s", err.msg);
    } else {
      const int32_t n1 = dae_spec_canonical(&S, NULL, 0) + 1;
      char *c1 = (char *)malloc((size_t)n1), *c2;
      dae_numag_serie R;
      dae_status st;
      dae_spec_canonical(&S, c1, n1);
      dae_test_ok(&T, strstr(c1, "\"ynu_table\":[") != NULL, "5. o canonico carrega ynu_table");
      if (dae_spec_parse(&S2, c1, &err) == DAE_OK) {
        const int32_t n2 = dae_spec_canonical(&S2, NULL, 0) + 1;
        c2 = (char *)malloc((size_t)n2);
        dae_spec_canonical(&S2, c2, n2);
        dae_test_ok(&T, strcmp(c1, c2) == 0, "5. canonico idempotente com ynu_table");
        free(c2);
        dae_spec_free(&S2);
      } else {
        dae_test_ok(&T, 0, "5. canonico recusado: %s", err.msg);
      }
      free(c1);
      st = dae_numag_run(&S, &R, NULL, NULL);
      if (st == DAE_OK) {
        double pior = 0.0, x;
        int32_t i;
        for (i = 0; i < R.nt; ++i) {
          x = sin(om * R.s_km[i] * TK);
          if (fabs(R.p_r[i] - x * x) > pior) pior = fabs(R.p_r[i] - x * x);
        }
        dae_test_ok(&T, pior < 1e-10,
                    "5. Y_e = 0,3 e Y_nue = 0,025 uniformes: P_R = sin^2(Omega s) (pior %.1e)", pior);
        dae_numag_serie_free(&R);
      } else {
        dae_test_ok(&T, 0, "5. dae_numag_run: %s", dae_strerror(st));
      }
      dae_spec_free(&S);
    }
    free(js);

    js = spec_json(L, 16, NULL);
    if (dae_spec_parse(&S, js, &err) == DAE_OK) {
      const int32_t n1 = dae_spec_canonical(&S, NULL, 0) + 1;
      char *c1 = (char *)malloc((size_t)n1);
      dae_numag_serie R;
      dae_spec_canonical(&S, c1, n1);
      dae_test_ok(&T, strstr(c1, "ynu_table") == NULL,
                  "5. sem fundo, o canonico nao ganha ynu_table");
      if (dae_numag_run(&S, &R, NULL, NULL) == DAE_OK) {
        dae_test_ok(&T, R.p_r[R.nt - 1] < 1e-8,
                    "5. companheira: sem o fundo, Y_e = 0,3 fica em nu_L (P_R = %.1e)", R.p_r[R.nt - 1]);
        dae_numag_serie_free(&R);
      }
      free(c1);
      dae_spec_free(&S);
    } else {
      dae_test_ok(&T, 0, "5. spec sem ynu_table: %s", err.msg);
    }
    free(js);

    js = spec_json(L, 16, ",\"ynu_table\":[0.01]");
    dae_test_ok(&T, dae_spec_parse(&S, js, &err) != DAE_OK && strstr(err.msg, "ynu_table"),
                "5. ynu_table de tamanho errado recusado (%s)", err.msg);
    free(js);
    js = spec_json(L, 16, ",\"ynu_table\":[0.01,1.5]");
    dae_test_ok(&T, dae_spec_parse(&S, js, &err) != DAE_OK && strstr(err.msg, "Y_nu"),
                "5. |Y_nu| > 1 recusado (%s)", err.msg);
    free(js);
  }

  return dae_test_end(&T);
}
