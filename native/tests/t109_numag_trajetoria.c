/* t109_numag_trajetoria.c — ACEITAÇÃO 109: a trajetória tabelada do setor
 * magnético e os blocos "magnetic"/"trajectory" do spec.
 *
 *  1. INTERPOLAÇÃO da tabela: exata nos nós; rho geométrico entre nós (ln rho
 *     linear); Y_e e b_t lineares; segura a borda fora da tabela.
 *  2. O PARSER RECUSA, com mensagem, o que define o problema errado: versão
 *     < 4, um bloco sem o outro, mistura com "neutrino", base transversal não
 *     declarada, s_km não crescente, tamanhos diferentes, t1 além da tabela,
 *     Majorana. Companheira: o spec correto é aceito.
 *  3. CANÔNICO: parse -> canônico -> parse -> canônico é idempotente, com as
 *     tabelas inteiras dentro.
 *  4. PONTA A PONTA: spec -> dae_numag_run em meio uniforme na ressonância
 *     (Y_e = 1/3) reproduz sin^2(Omega s) em cada linha; dae_run recusa o
 *     spec magnético (senão rodaria o grafo padrão em silêncio).
 *  5. CAMPO GIRANTE TABELADO por componentes (b1 = B cos phi, b2 = B sin phi,
 *     2001 nós): em Delta = -phi'/T converte tudo, em Delta = +phi'/T fica em
 *     0,2. O sinal de t107 sobrevive à interpolação de componentes, que é a
 *     escolha que evita o salto de 2 pi do ângulo.
 */
#include "harness.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define RHO 1.0e14

/* ------------------------------------------------- montagem do spec ---- */

typedef struct { char *b; size_t n, cap; } texto;

static void put(texto *t, const char *s)
{
  const size_t m = strlen(s);
  if (t->n + m + 1 > t->cap) {
    t->cap = 2 * (t->n + m + 1);
    t->b = (char *)realloc(t->b, t->cap);
  }
  memcpy(t->b + t->n, s, m + 1);
  t->n += m;
}

static void putv(texto *t, const char *nome, const double *v, int32_t n)
{
  char tmp[64];
  int32_t i;
  snprintf(tmp, sizeof tmp, ",\"%s\":[", nome); put(t, tmp);
  for (i = 0; i < n; ++i) { snprintf(tmp, sizeof tmp, i ? ",%.17g" : "%.17g", v[i]); put(t, tmp); }
  put(t, "]");
}

/* Monta um spec magnético. `extra_mag` e `extra_traj` entram crus dentro dos
   blocos (para os casos de recusa); `basis` NULL omite a declaração. */
static char *spec(int versao, double mu, double t1, int32_t nt,
                  const double *s, const double *rho, const double *ye,
                  const double *b1, const double *b2, int32_t n,
                  const char *basis, const char *extra_mag, const char *extra_raiz)
{
  texto t = { NULL, 0, 0 };
  char tmp[256];
  snprintf(tmp, sizeof tmp, "{\"format_version\":%d,\"seed\":1", versao); put(&t, tmp);
  snprintf(tmp, sizeof tmp, ",\"magnetic\":{\"nature\":\"dirac\",\"mu_nu_bohr\":%.17g,"
           "\"energy_MeV\":10.0,\"initial_helicity\":\"left\",\"tolerance\":1e-11%s}",
           mu, extra_mag ? extra_mag : ""); put(&t, tmp);
  put(&t, ",\"trajectory\":{");
  if (basis) { snprintf(tmp, sizeof tmp, "\"transverse_basis\":\"%s\"", basis); put(&t, tmp); }
  else put(&t, "\"_\":0");          /* marcador removido abaixo */
  putv(&t, "s_km", s, n);
  putv(&t, "rho_table", rho, n);
  putv(&t, "ye_table", ye, n);
  putv(&t, "b1_gauss", b1, n);
  putv(&t, "b2_gauss", b2, n);
  put(&t, "}");
  snprintf(tmp, sizeof tmp, ",\"time\":{\"t1\":%.17g,\"nt\":%d}", t1, (int)nt); put(&t, tmp);
  if (extra_raiz) put(&t, extra_raiz);
  put(&t, ",\"realizations\":1}");
  if (!basis) {
    /* tira o marcador e a vírgula seguinte: "{"_":0,"s_km" -> "{"s_km" */
    char *p = strstr(t.b, "\"_\":0,");
    memmove(p, p + 6, strlen(p + 6) + 1);
  }
  return t.b;
}

static int aceita(const char *json, dae_spec *S, char *msg, size_t cap)
{
  dae_error err;
  const dae_status st = dae_spec_parse(S, json, &err);
  if (st != DAE_OK) { snprintf(msg, cap, "%s", err.msg); return 0; }
  msg[0] = '\0';
  return 1;
}

static double exata(double delta, double omega, double t)
{
  const double w = sqrt(delta * delta + 4.0 * omega * omega);
  const double s = sin(0.5 * t * w);
  return (w == 0.0) ? 0.0 : 4.0 * omega * omega / (w * w) * s * s;
}

int main(void)
{
  dae_test T;
  const double TK = DAE_NU_T_POR_KM;
  const double omega = 1e-12 * DAE_NUMAG_MU_B_EV_POR_T * 1e15 * DAE_NUMAG_T_POR_GAUSS;
  char msg[256];
  dae_spec S;
  int32_t i;

  dae_test_begin(&T, "t109 trajetoria tabelada e os blocos magnetic/trajectory");

  /* ------------------------------------------------- 1. interpolação ---- */
  {
    static const double s[3] = { 0.0, 1.0, 3.0 };
    static const double rho[3] = { 1.0e14, 1.0e12, 1.0e10 };
    static const double ye[3] = { 0.1, 0.3, 0.5 };
    static const double b1[3] = { 1.0, -1.0, 2.0 };
    static const double b2[3] = { 0.0, 2.0, 4.0 };
    dae_numag_trajetoria tr;
    dae_numag_meio m;
    tr.s_km = s; tr.rho = rho; tr.ye = ye; tr.b1 = b1; tr.b2 = b2;
    tr.razao_energia = NULL; tr.ynu = NULL; tr.n = 3;
    dae_numag_meio_tabela(&tr, 1.0, &m);
    dae_test_ok(&T, m.rho == 1.0e12 && m.ye == 0.3 && m.b_t[0] == -1.0 && m.b_t[1] == 2.0,
                "1. exata no no");
    dae_numag_meio_tabela(&tr, 2.0, &m);
    dae_test_near(&T, m.rho / 1.0e11, 1.0, 1e-12, "1. rho geometrico entre nos (1e12, 1e10 -> 1e11)");
    dae_test_near(&T, m.ye, 0.4, 1e-15, "1. Ye linear");
    dae_test_ok(&T, fabs(m.b_t[0] - 0.5) < 1e-15 && fabs(m.b_t[1] - 3.0) < 1e-15, "1. b_t linear");
    dae_test_ok(&T, m.razao_energia == 1.0, "1. razao_energia ausente vale 1");
    dae_numag_meio_tabela(&tr, -5.0, &m);
    dae_test_ok(&T, m.rho == 1.0e14 && m.ye == 0.1, "1. antes da tabela segura a borda");
    dae_numag_meio_tabela(&tr, 9.0, &m);
    dae_test_ok(&T, m.rho == 1.0e10 && m.ye == 0.5 && m.b_t[1] == 4.0, "1. depois da tabela segura a borda");
  }

  /* --------------------------------------- 2. recusas e 3. canônico ---- */
  {
    static const double s[3] = { 0.0, 1.0, 2.0 }, rho[3] = { RHO, RHO, RHO };
    static const double ye[3] = { 0.3, 0.3, 0.3 }, b1[3] = { 1e15, 1e15, 1e15 };
    static const double b2[3] = { 0.0, 0.0, 0.0 }, sruim[3] = { 0.0, 2.0, 1.0 };
    static const char *NEU = ",\"neutrino\":{\"n_flavors\":2,\"energy_MeV\":10,\"flavor\":\"e\"}";
    struct { const char *nome; char *json; const char *espera; } caso[8];
    char *bom = spec(4, 1e-12, 1.5, 3, s, rho, ye, b1, b2, 3, "parallel_transported", NULL, NULL);
    int32_t c;
    caso[0].nome = "versao 3";            caso[0].json = spec(3, 1e-12, 1.5, 3, s, rho, ye, b1, b2, 3, "parallel_transported", NULL, NULL); caso[0].espera = "format_version 4";
    caso[1].nome = "com neutrino";        caso[1].json = spec(4, 1e-12, 1.5, 3, s, rho, ye, b1, b2, 3, "parallel_transported", NULL, NEU); caso[1].espera = "exclusivo";
    caso[2].nome = "base nao declarada";  caso[2].json = spec(4, 1e-12, 1.5, 3, s, rho, ye, b1, b2, 3, NULL, NULL, NULL); caso[2].espera = "transverse_basis";
    caso[3].nome = "base desconhecida";   caso[3].json = spec(4, 1e-12, 1.5, 3, s, rho, ye, b1, b2, 3, "local", NULL, NULL); caso[3].espera = "desconhecido";
    caso[4].nome = "s nao crescente";     caso[4].json = spec(4, 1e-12, 0.5, 3, sruim, rho, ye, b1, b2, 3, "parallel_transported", NULL, NULL); caso[4].espera = "crescente";
    caso[5].nome = "tamanho diferente";   caso[5].json = NULL;
    caso[6].nome = "t1 alem da tabela";   caso[6].json = spec(4, 1e-12, 2.5, 3, s, rho, ye, b1, b2, 3, "parallel_transported", NULL, NULL); caso[6].espera = "caber";
    caso[7].nome = "majorana";            caso[7].json = NULL;
    /* caso 5: ye_table com um valor a menos que s_km */
    caso[5].json = spec(4, 1e-12, 1.5, 3, s, rho, ye, b1, b2, 3, "parallel_transported", NULL, NULL);
    { char *p = strstr(caso[5].json, "\"ye_table\":[");
      char *v = strchr(p, ',');                    /* remove o primeiro valor */
      memmove(p + 12, v + 1, strlen(v + 1) + 1); }
    caso[5].espera = "tamanho";
    free(caso[7].json);
    {
      /* troca "dirac" por "majorana" no spec correto */
      texto t = { NULL, 0, 0 };
      char *base = spec(4, 1e-12, 1.5, 3, s, rho, ye, b1, b2, 3, "parallel_transported", NULL, NULL);
      char *p = strstr(base, "\"dirac\"");
      *p = '\0';
      put(&t, base); put(&t, "\"majorana\""); put(&t, p + 7);
      free(base);
      caso[7].json = t.b;
      caso[7].espera = "Dirac";
    }

    for (c = 0; c < 8; ++c) {
      const int ok = aceita(caso[c].json, &S, msg, sizeof msg);
      dae_test_ok(&T, !ok && (!caso[c].espera || strstr(msg, caso[c].espera)),
                  "2. recusa: %s (\"%s\")", caso[c].nome, msg);
      if (ok) dae_spec_free(&S);
      free(caso[c].json);
    }

    dae_test_ok(&T, aceita(bom, &S, msg, sizeof msg), "2. companheira: o spec correto e aceito (%s)", msg);
    {
      const int32_t n1 = dae_spec_canonical(&S, NULL, 0) + 1;
      char *c1 = (char *)malloc((size_t)n1), *c2;
      dae_spec S2;
      int32_t n2;
      dae_spec_canonical(&S, c1, n1);
      dae_test_ok(&T, strstr(c1, "\"transverse_basis\":\"parallel_transported\"") != NULL &&
                      strstr(c1, "\"b1_gauss\":[") != NULL,
                  "3. o canonico carrega a base e as tabelas");
      dae_test_ok(&T, aceita(c1, &S2, msg, sizeof msg), "3. o canonico e aceito (%s)", msg);
      n2 = dae_spec_canonical(&S2, NULL, 0) + 1;
      c2 = (char *)malloc((size_t)n2);
      dae_spec_canonical(&S2, c2, n2);
      dae_test_ok(&T, n1 == n2 && strcmp(c1, c2) == 0, "3. canonico idempotente (%d bytes)", (int)n1);
      dae_spec_free(&S2);
      free(c1); free(c2);
    }
    dae_spec_free(&S);
    free(bom);
  }

  /* ------------------------------------------------ 4. ponta a ponta ---- */
  {
    /* meio uniforme em Y_e = 1/3: Delta = 0, P_R = sin^2(Omega s) */
    const double L = PI / omega / TK;                /* um periodo completo */
    double s[2], rho[2] = { RHO, RHO }, ye[2] = { 1.0 / 3.0, 1.0 / 3.0 };
    double b1[2] = { 1e15, 1e15 }, b2[2] = { 0.0, 0.0 };
    char *json;
    dae_numag_serie R;
    dae_graph G;
    dae_series Rs;
    dae_status st;
    double pior = 0.0;
    s[0] = 0.0; s[1] = L;
    json = spec(4, 1e-12, L, 16, s, rho, ye, b1, b2, 2, "parallel_transported", NULL, NULL);
    if (!aceita(json, &S, msg, sizeof msg)) { dae_test_ok(&T, 0, "4. spec: %s", msg); free(json); return dae_test_end(&T); }
    st = dae_numag_run(&S, &R, NULL, NULL);
    dae_test_ok(&T, st == DAE_OK && R.nt == 16, "4. dae_numag_run (%s)", dae_strerror(st));
    if (st == DAE_OK) {
      for (i = 0; i < R.nt; ++i) {
        const double d = fabs(R.p_r[i] - exata(0.0, omega, R.s_km[i] * TK));
        if (d > pior) pior = d;
      }
      dae_test_ok(&T, pior < 1e-10, "4. P_R = sin^2(Omega s) nas 16 linhas (pior %.2e)", pior);
      dae_test_ok(&T, fabs(R.s_km[R.nt - 1] - L) < 1e-15 * L + 1e-18, "4. a ultima linha esta em s = t1");
      {
        const int32_t n = dae_numag_csv(&S, &R, NULL, 0) + 1;
        char *csv = (char *)malloc((size_t)n);
        dae_numag_csv(&S, &R, csv, n);
        dae_test_ok(&T, strstr(csv, "#! setor magnetico") && strstr(csv, "s_km,norm,p_L,p_R,subdiv\n") &&
                        strstr(csv, "#! spec {\"format_version\":4"),
                    "4. CSV com procedencia e colunas");
        free(csv);
      }
      dae_numag_serie_free(&R);
    }
    st = dae_run(&S, S.seed, &G, NULL, &Rs, NULL, NULL);
    dae_test_ok(&T, st == DAE_ERR_PARAM, "4. dae_run recusa o spec magnetico (%s)", dae_strerror(st));
    dae_spec_free(&S);
    free(json);
  }

  /* --------------------------------- 5. campo girante tabelado ---- */
  {
    enum { N = 2001 };
    const double dphi = 2.0 * omega * TK;              /* rad/km: phi' = 2 Omega */
    const double L = 0.5 * PI / omega / TK;            /* primeiro maximo        */
    static double s[N], rho[N], ye_m[N], ye_p[N], b1[N], b2[N];
    const double ye_menos = (2.0 * (-2.0 * omega) / (DAE_NU_V_POR_RHO * RHO) + 1.0) / 3.0;
    const double ye_mais = (2.0 * (+2.0 * omega) / (DAE_NU_V_POR_RHO * RHO) + 1.0) / 3.0;
    double pmenos = -1.0, pmais_max = 0.0;
    int lado;
    for (i = 0; i < N; ++i) {
      s[i] = L * (double)i / (double)(N - 1);
      rho[i] = RHO; ye_m[i] = ye_menos; ye_p[i] = ye_mais;
      b1[i] = 1e15 * cos(dphi * s[i]);
      b2[i] = 1e15 * sin(dphi * s[i]);
    }
    for (lado = 0; lado < 2; ++lado) {
      char *json = spec(4, 1e-12, L, 8, s, rho, lado ? ye_p : ye_m, b1, b2, N,
                        "parallel_transported", NULL, NULL);
      dae_numag_serie R;
      if (!aceita(json, &S, msg, sizeof msg)) { dae_test_ok(&T, 0, "5. spec: %s", msg); free(json); break; }
      if (dae_numag_run(&S, &R, NULL, NULL) == DAE_OK) {
        if (lado == 0) pmenos = R.p_r[R.nt - 1];
        else for (i = 0; i < R.nt; ++i) if (R.p_r[i] > pmais_max) pmais_max = R.p_r[i];
        dae_numag_serie_free(&R);
      }
      dae_spec_free(&S);
      free(json);
    }
    /* Resolução: entre nós, a corda encurta |B| por ~(dphi ds)^2/8, com
       dphi ds = pi/2000 -> ~3e-7. */
    dae_test_near(&T, pmenos, 1.0, 1e-5,
                  "5. tabelado: Delta = -phi'/T converte tudo (sinal de t107, resolucao 1e-5)");
    dae_test_ok(&T, pmais_max < 0.21, "5. tabelado: Delta = +phi'/T fica em 0,2 (max %.4f)", pmais_max);
  }

  return dae_test_end(&T);
}
