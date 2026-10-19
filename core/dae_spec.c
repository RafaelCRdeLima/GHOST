#include "dae_spec.h"

/* EMENDA A REGRA 3 DO NUCLEO (CONVENTIONS.md, parte 1): este e o unico arquivo
 * do nucleo que inclui <stdio.h>, e so por snprintf. Formatar para um buffer
 * nao e entrada e saida — nenhum printf, nenhum arquivo. A alternativa seria
 * escrever a serializacao canonica fora do nucleo, e ai ela existiria duas
 * vezes: uma no TypeScript e outra no main do .cpp exportado. Duas versoes do
 * JSON canonico e exatamente a divergencia silenciosa que a amalgamacao existe
 * para impedir. */
#include <stdio.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------ varredura -- */

typedef struct {
  const char *s;
  int32_t     i, linha, coluna;
  dae_error  *err;
  int         ok;
} dae_js;

static void js_erro(dae_js *j, const char *msg)
{
  if (!j->ok) return;                       /* guarda o primeiro erro */
  j->ok = 0;
  if (!j->err) return;
  j->err->line = j->linha;
  j->err->col = j->coluna;
  { size_t n = strlen(msg);
    if (n > sizeof(j->err->msg) - 1) n = sizeof(j->err->msg) - 1;
    memcpy(j->err->msg, msg, n);
    j->err->msg[n] = '\0'; }
}

static void js_avanca(dae_js *j)
{
  if (j->s[j->i] == '\n') { ++j->linha; j->coluna = 1; }
  else ++j->coluna;
  ++j->i;
}

static void js_ws(dae_js *j)
{
  for (;;) {
    const char c = j->s[j->i];
    if (c == ' ' || c == '\t' || c == '\n' || c == '\r') js_avanca(j);
    else return;
  }
}

static int js_ch(dae_js *j, char c, const char *msg)
{
  js_ws(j);
  if (j->s[j->i] != c) { js_erro(j, msg); return 0; }
  js_avanca(j);
  return 1;
}

static int js_texto(dae_js *j, char *out, int32_t cap)
{
  int32_t k = 0;
  js_ws(j);
  if (!js_ch(j, '"', "esperava string")) return 0;
  for (;;) {
    const char c = j->s[j->i];
    if (c == '\0') { js_erro(j, "string sem fechamento"); return 0; }
    if (c == '"') { js_avanca(j); break; }
    if (c == '\\') {
      js_avanca(j);
      { const char e = j->s[j->i];
        char d = e;
        if (e == 'n') d = '\n'; else if (e == 't') d = '\t';
        else if (e == 'r') d = '\r'; else if (e == 'b') d = '\b';
        else if (e == 'f') d = '\f';
        else if (e != '"' && e != '\\' && e != '/') { js_erro(j, "escape nao suportado"); return 0; }
        if (k + 1 < cap) out[k++] = d;
        js_avanca(j); }
      continue;
    }
    if (k + 1 < cap) out[k++] = c;
    js_avanca(j);
  }
  out[k] = '\0';
  return 1;
}

static int js_numero(dae_js *j, double *out)
{
  char buf[64];
  int32_t k = 0;
  js_ws(j);
  while (k + 1 < (int32_t)sizeof(buf)) {
    const char c = j->s[j->i];
    if ((c >= '0' && c <= '9') || c == '+' || c == '-' || c == '.' ||
        c == 'e' || c == 'E') { buf[k++] = c; js_avanca(j); }
    else break;
  }
  buf[k] = '\0';
  if (k == 0) { js_erro(j, "esperava numero"); return 0; }
  *out = strtod(buf, NULL);
  if (!(*out == *out)) { js_erro(j, "numero invalido (NaN)"); return 0; }
  return 1;
}

static int js_bool(dae_js *j, int *out)
{
  js_ws(j);
  if (strncmp(j->s + j->i, "true", 4) == 0) {
    int k; for (k = 0; k < 4; ++k) js_avanca(j);
    *out = 1; return 1;
  }
  if (strncmp(j->s + j->i, "false", 5) == 0) {
    int k; for (k = 0; k < 5; ++k) js_avanca(j);
    *out = 0; return 1;
  }
  js_erro(j, "esperava true ou false");
  return 0;
}

/* Iteração de objeto: devolve 0 quando acabou (ou em erro). */
static int js_obj_inicio(dae_js *j) { return js_ch(j, '{', "esperava '{'"); }

static int js_obj_chave(dae_js *j, char *k, int32_t cap, int *primeiro)
{
  js_ws(j);
  if (j->s[j->i] == '}') { js_avanca(j); return 0; }
  if (!*primeiro) { if (!js_ch(j, ',', "esperava ',' ou '}'")) return 0; }
  *primeiro = 0;
  js_ws(j);
  if (j->s[j->i] == '}') { js_avanca(j); return 0; }
  if (!js_texto(j, k, cap)) return 0;
  if (!js_ch(j, ':', "esperava ':'")) return 0;
  return 1;
}

static int js_arr_inicio(dae_js *j) { return js_ch(j, '[', "esperava '['"); }

static int js_arr_item(dae_js *j, int *primeiro)
{
  js_ws(j);
  if (j->s[j->i] == ']') { js_avanca(j); return 0; }
  if (!*primeiro) { if (!js_ch(j, ',', "esperava ',' ou ']'")) return 0; }
  *primeiro = 0;
  js_ws(j);
  if (j->s[j->i] == ']') { js_avanca(j); return 0; }
  return 1;
}

/* --------------------------------------------------------------- enums -- */

typedef struct { const char *nome; int valor; } dae_mapa;

static const dae_mapa MAPA_GER[] = {
  { "microtubule", DAE_G_MICROTUBULE }, { "path", DAE_G_PATH },
  { "cycle", DAE_G_CYCLE }, { "complete", DAE_G_COMPLETE },
  { "hypercube", DAE_G_HYPERCUBE }, { "grid2d", DAE_G_GRID2D },
  { "sbm", DAE_G_SBM }, { "edgelist", DAE_G_EDGELIST }, { 0, 0 }
};
static const dae_mapa MAPA_CONN[] = {
  { "rewire", DAE_REWIRE }, { "add", DAE_ADD }, { 0, 0 }
};
static const dae_mapa MAPA_HAM[] = {
  { "adjacency", DAE_H_ADJACENCY }, { "laplacian", DAE_H_LAPLACIAN }, { 0, 0 }
};
/* Sem entrada para INDEFINIDA: ela nao e escrevivel no arquivo, e sim o estado
   de "o bloco existe e nao disse". */
static const dae_mapa MAPA_SAIDA[] = {
  { "accumulate_rho", DAE_SAIDA_ACUMULAR_RHO },
  { "archive_psi",    DAE_SAIDA_ARQUIVAR_PSI }, { 0, 0 }
};
static const dae_mapa MAPA_PERFIL[] = {
  { "constant", DAE_NU_PERFIL_CONSTANTE }, { "power_law", DAE_NU_PERFIL_POTENCIA },
  { "table", DAE_NU_PERFIL_TABELA },
  { 0, 0 }
};
static const dae_mapa MAPA_SABOR[] = {
  { "e", 0 }, { "mu", 1 }, { "tau", 2 }, { 0, 0 }
};
static const dae_mapa MAPA_NATUREZA[] = {
  { "dirac", DAE_NUMAG_DIRAC }, { "majorana", DAE_NUMAG_MAJORANA }, { 0, 0 }
};
static const dae_mapa MAPA_INTEGRADOR[] = {
  { "brute_force", 0 }, { "windows", 1 }, { 0, 0 }
};
static const dae_mapa MAPA_COERENCIA[] = {
  { "coherent", 0 }, { "incoherent", 1 }, { 0, 0 }
};
static const dae_mapa MAPA_HELICIDADE[] = {
  { "left", 0 }, { "right", 1 }, { 0, 0 }
};
/* Um valor só, e é o ponto: quem escreve o spec DECLARA que b1/b2 estão na
   base transportada ao longo do raio. Uma base recalculada ponto a ponto por
   regra local cria ressonância artificial (CONVENTIONS.md 13.2). */
static const dae_mapa MAPA_BASE[] = {
  { "parallel_transported", 1 }, { 0, 0 }
};
static const dae_mapa MAPA_NORM[] = {
  { "none", DAE_NORM_NONE }, { "spectral", DAE_NORM_SPECTRAL },
  { "mean_degree", DAE_NORM_MEAN_DEGREE }, { 0, 0 }
};

static const char *nome_de(const dae_mapa *m, int v)
{
  for (; m->nome; ++m) if (m->valor == v) return m->nome;
  return "?";
}

static int js_enum(dae_js *j, const dae_mapa *m, int *out)
{
  char buf[40];
  if (!js_texto(j, buf, (int32_t)sizeof(buf))) return 0;
  for (; m->nome; ++m) {
    if (strcmp(m->nome, buf) == 0) { *out = m->valor; return 1; }
  }
  js_erro(j, "valor desconhecido para este campo");
  return 0;
}

/* --------------------------------------------------------------- spec --- */

void dae_spec_default(dae_spec *S)
{
  if (!S) return;
  memset(S, 0, sizeof(*S));
  S->format_version = DAE_FORMAT_VERSION;
  S->seed = 12345ULL;
  dae_gen_params_default(&S->gen);
  S->ham = DAE_H_ADJACENCY;
  S->gamma = 1.0;
  S->norm = DAE_NORM_SPECTRAL;      /* ligada por padrao: CONVENTIONS 3.2 */
  S->lanczos_steps = 40;
  S->init_site = 0;
  S->t1 = 50.0;
  S->nt = 500;
  S->target = -1;
  S->want_pop = 1;
  S->want_conc_mod = 1;
  S->pop_stride = 1;
  S->nu_ativo = 0;
  dae_nu_params_default(&S->nu);
  S->nu_rho = 0.0;
  S->nu_ye = 0.5;
  S->nu_sabor = 0;
  S->nu_perfil_ativo = 0;
  dae_nu_perfil_default(&S->nu_perfil);
  S->perf_r = NULL; S->perf_rho = NULL; S->perf_ye = NULL; S->perf_n = 0;
  S->ruido_amp = NULL; S->ruido_fase = NULL;
  S->nu_tol = 1e-6;
  S->scan_pontos = 0;
  S->scan_min = 0.0; S->scan_max = 0.0;
  S->scan_log = 0; S->scan_ambos = 0;
  S->realizations = 1;
  S->n_traj = 0;
  S->gamma_deph = 0.0;
  S->rho_stride = 1;
  S->saida_traj = DAE_SAIDA_INDEFINIDA;
  S->numag_ativo = 0;
  S->numag_tem_traj = 0;
  dae_numag_params_default(&S->numag);
  S->numag_tol = 1e-10;
  S->numag_helicidade = 0;
  S->numag_base_declarada = 0;
  S->numag_integrador = 0;
  S->numag_coerencia = -1;
  dae_numag_janelas_default(&S->numag_jan);
  S->numag_gamma_eV = 0.0;
  S->traj_gamma = NULL;
  S->numag_n_traj = 0;
  S->traj_s = NULL; S->traj_rho = NULL; S->traj_ye = NULL;
  S->traj_b1 = NULL; S->traj_b2 = NULL; S->traj_e = NULL;
  S->traj_ynu = NULL;
  S->traj_n = 0;
  memset(&S->numag_traj, 0, sizeof(S->numag_traj));
}

void dae_spec_free(dae_spec *S)
{
  if (!S) return;
  free(S->ei); free(S->ej); free(S->w);
  free(S->init_re); free(S->init_im);
  free(S->perf_r); free(S->perf_rho); free(S->perf_ye);
  free(S->ruido_amp); free(S->ruido_fase);
  S->ei = NULL; S->ej = NULL; S->w = NULL;
  S->init_re = NULL; S->init_im = NULL;
  S->perf_r = NULL; S->perf_rho = NULL; S->perf_ye = NULL;
  S->ruido_amp = NULL; S->ruido_fase = NULL;
  S->nu_perfil.ruido_amp = NULL; S->nu_perfil.ruido_fase = NULL;
  /* os ponteiros do perfil apontavam para os arrays acima: zera junto, senao
     ficam pendurados e o proximo uso le memoria liberada. */
  S->nu_perfil.r_tab = NULL; S->nu_perfil.rho_tab = NULL; S->nu_perfil.ye_tab = NULL;
  S->nu_perfil.n_tab = 0;
  S->n_edges = 0; S->init_n = 0; S->perf_n = 0;
  free(S->traj_s); free(S->traj_rho); free(S->traj_ye);
  free(S->traj_b1); free(S->traj_b2); free(S->traj_e); free(S->traj_gamma);
  free(S->traj_ynu);
  S->traj_gamma = NULL;
  S->traj_s = NULL; S->traj_rho = NULL; S->traj_ye = NULL;
  S->traj_b1 = NULL; S->traj_b2 = NULL; S->traj_e = NULL;
  S->traj_ynu = NULL;
  S->traj_n = 0;
  memset(&S->numag_traj, 0, sizeof(S->numag_traj));
}

int dae_spec_hash_confere(const dae_spec *S)
{
  return S && S->core_hash[0] && strcmp(S->core_hash, DAE_CORE_HASH) == 0;
}

/* Le um array de numeros para *out (realocado). Devolve 0 em erro. */
static int le_vetor_num(dae_js *j, double **out, int32_t *n_out)
{
  int32_t cap = 64, n = 0;
  int primeiro = 1;
  double *v = (double *)malloc((size_t)cap * sizeof(double));
  if (!v) { js_erro(j, "memoria"); return 0; }
  if (!js_arr_inicio(j)) { free(v); return 0; }
  for (;;) {
    double x;
    if (!js_arr_item(j, &primeiro)) break;
    if (!js_numero(j, &x)) { free(v); return 0; }
    if (n == cap) {
      double *nv = (double *)realloc(v, (size_t)(cap * 2) * sizeof(double));
      if (!nv) { free(v); js_erro(j, "memoria"); return 0; }
      v = nv; cap *= 2;
    }
    v[n++] = x;
  }
  free(*out);
  *out = v; *n_out = n;
  return j->ok;
}

static int le_arestas(dae_js *j, dae_spec *S)
{
  int32_t cap = 64;
  int primeiro = 1;
  S->ei = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  S->ej = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
  S->w  = (double  *)malloc((size_t)cap * sizeof(double));
  if (!S->ei || !S->ej || !S->w) { js_erro(j, "memoria"); return 0; }
  if (!js_arr_inicio(j)) return 0;
  for (;;) {
    double a, b, p = 1.0;
    int interno = 1;
    if (!js_arr_item(j, &primeiro)) break;
    if (!js_arr_inicio(j)) return 0;
    if (!js_arr_item(j, &interno) || !js_numero(j, &a)) return 0;
    if (!js_arr_item(j, &interno) || !js_numero(j, &b)) return 0;
    if (js_arr_item(j, &interno)) {              /* peso opcional */
      if (!js_numero(j, &p)) return 0;
      if (js_arr_item(j, &interno)) { js_erro(j, "aresta com mais de 3 campos"); return 0; }
    }
    if (S->n_edges == cap) {
      int32_t nc = cap * 2;
      int32_t *ni = (int32_t *)realloc(S->ei, (size_t)nc * sizeof(int32_t));
      int32_t *nj = (int32_t *)realloc(S->ej, (size_t)nc * sizeof(int32_t));
      double  *nw = (double  *)realloc(S->w,  (size_t)nc * sizeof(double));
      if (!ni || !nj || !nw) { js_erro(j, "memoria"); return 0; }
      S->ei = ni; S->ej = nj; S->w = nw; cap = nc;
    }
    S->ei[S->n_edges] = (int32_t)a;
    S->ej[S->n_edges] = (int32_t)b;
    S->w[S->n_edges] = p;
    ++S->n_edges;
  }
  return j->ok;
}

static int le_vetor_inicial(dae_js *j, dae_spec *S)
{
  int32_t cap = 64;
  int primeiro = 1;
  S->init_re = (double *)malloc((size_t)cap * sizeof(double));
  S->init_im = (double *)malloc((size_t)cap * sizeof(double));
  if (!S->init_re || !S->init_im) { js_erro(j, "memoria"); return 0; }
  if (!js_arr_inicio(j)) return 0;
  for (;;) {
    double re, im;
    int interno = 1;
    if (!js_arr_item(j, &primeiro)) break;
    if (!js_arr_inicio(j)) return 0;
    if (!js_arr_item(j, &interno) || !js_numero(j, &re)) return 0;
    if (!js_arr_item(j, &interno) || !js_numero(j, &im)) return 0;
    if (js_arr_item(j, &interno)) { js_erro(j, "amplitude com mais de 2 campos"); return 0; }
    if (S->init_n == cap) {
      int32_t nc = cap * 2;
      double *nr = (double *)realloc(S->init_re, (size_t)nc * sizeof(double));
      double *ni = (double *)realloc(S->init_im, (size_t)nc * sizeof(double));
      if (!nr || !ni) { js_erro(j, "memoria"); return 0; }
      S->init_re = nr; S->init_im = ni; cap = nc;
    }
    S->init_re[S->init_n] = re;
    S->init_im[S->init_n] = im;
    ++S->init_n;
  }
  S->init_site = -1;
  return j->ok;
}

#define CHAVE(x) (strcmp(k, x) == 0)

static int le_params_grafo(dae_js *j, dae_spec *S)
{
  char k[48];
  int primeiro = 1;
  if (!js_obj_inicio(j)) return 0;
  while (js_obj_chave(j, k, (int32_t)sizeof(k), &primeiro)) {
    double v; int b, e;
    if      (CHAVE("n_par"))      { if (!js_numero(j, &v)) return 0; S->gen.n_par = (int32_t)v; }
    else if (CHAVE("n_perp"))     { if (!js_numero(j, &v)) return 0; S->gen.n_perp = (int32_t)v; }
    else if (CHAVE("seam_shift")) { if (!js_numero(j, &v)) return 0; S->gen.seam_shift = (int32_t)v; }
    else if (CHAVE("longitudinal_closed")) { if (!js_bool(j, &b)) return 0; S->gen.longitudinal_closed = b; }
    else if (CHAVE("j_par"))      { if (!js_numero(j, &v)) return 0; S->gen.j_par = v; }
    else if (CHAVE("j_perp"))     { if (!js_numero(j, &v)) return 0; S->gen.j_perp = v; }
    else if (CHAVE("n_modules"))  { if (!js_numero(j, &v)) return 0; S->gen.n_modules = (int32_t)v; }
    else if (CHAVE("ws_p"))       { if (!js_numero(j, &v)) return 0; S->gen.ws_p = v; }
    else if (CHAVE("conn_mode"))  { if (!js_enum(j, MAPA_CONN, &e)) return 0; S->gen.conn_mode = (dae_conn_mode)e; }
    else if (CHAVE("p_in"))       { if (!js_numero(j, &v)) return 0; S->gen.p_in = v; }
    else if (CHAVE("p_out"))      { if (!js_numero(j, &v)) return 0; S->gen.p_out = v; }
    else if (CHAVE("n"))          { if (!js_numero(j, &v)) return 0; S->gen.n = (int32_t)v; }
    else if (CHAVE("dim"))        { if (!js_numero(j, &v)) return 0; S->gen.dim = (int32_t)v; }
    else if (CHAVE("rows"))       { if (!js_numero(j, &v)) return 0; S->gen.rows = (int32_t)v; }
    else if (CHAVE("cols"))       { if (!js_numero(j, &v)) return 0; S->gen.cols = (int32_t)v; }
    else { js_erro(j, "chave desconhecida em graph.params"); return 0; }
  }
  return j->ok;
}

static int le_grafo(dae_js *j, dae_spec *S)
{
  char k[48];
  int primeiro = 1;
  if (!js_obj_inicio(j)) return 0;
  while (js_obj_chave(j, k, (int32_t)sizeof(k), &primeiro)) {
    double v; int e;
    if      (CHAVE("generator")) { if (!js_enum(j, MAPA_GER, &e)) return 0; S->gen.kind = (dae_gen_kind)e; }
    else if (CHAVE("params"))    { if (!le_params_grafo(j, S)) return 0; }
    else if (CHAVE("n"))         { if (!js_numero(j, &v)) return 0; S->gen.n = (int32_t)v; }
    else if (CHAVE("edges"))     { if (!le_arestas(j, S)) return 0; }
    else { js_erro(j, "chave desconhecida em graph"); return 0; }
  }
  return j->ok;
}

dae_status dae_spec_parse(dae_spec *S, const char *json, dae_error *err)
{
  dae_js j;
  char k[48];
  int primeiro = 1;
  int tem_traj = 0;

  if (!S || !json) return DAE_ERR_PARAM;
  dae_spec_default(S);
  if (err) { err->line = 0; err->col = 0; err->msg[0] = '\0'; }

  j.s = json; j.i = 0; j.linha = 1; j.coluna = 1; j.err = err; j.ok = 1;

  if (!js_obj_inicio(&j)) return DAE_ERR_JSON;
  while (js_obj_chave(&j, k, (int32_t)sizeof(k), &primeiro)) {
    double v; int b, e;
    if      (CHAVE("format_version")) { if (!js_numero(&j, &v)) break; S->format_version = (int32_t)v; }
    else if (CHAVE("core_hash"))      { if (!js_texto(&j, S->core_hash, (int32_t)sizeof(S->core_hash))) break; }
    else if (CHAVE("seed"))           { if (!js_numero(&j, &v)) break; S->seed = (uint64_t)v; }
    else if (CHAVE("graph"))          { if (!le_grafo(&j, S)) break; }
    else if (CHAVE("hamiltonian")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("kind"))          { if (!js_enum(&j, MAPA_HAM, &e)) break; S->ham = (dae_ham_kind)e; }
        else if (CHAVE("gamma"))         { if (!js_numero(&j, &v)) break; S->gamma = v; }
        else if (CHAVE("normalization")) { if (!js_enum(&j, MAPA_NORM, &e)) break; S->norm = (dae_norm_kind)e; }
        else if (CHAVE("lanczos_steps")) { if (!js_numero(&j, &v)) break; S->lanczos_steps = (int32_t)v; }
        else { js_erro(&j, "chave desconhecida em hamiltonian"); break; }
      }
    }
    else if (CHAVE("initial")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("site"))   { if (!js_numero(&j, &v)) break; S->init_site = (int32_t)v; }
        else if (CHAVE("vector")) { if (!le_vetor_inicial(&j, S)) break; }
        else { js_erro(&j, "chave desconhecida em initial"); break; }
      }
    }
    else if (CHAVE("time")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("t1")) { if (!js_numero(&j, &v)) break; S->t1 = v; }
        else if (CHAVE("nt")) { if (!js_numero(&j, &v)) break; S->nt = (int32_t)v; }
        else { js_erro(&j, "chave desconhecida em time"); break; }
      }
    }
    else if (CHAVE("observables")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("target"))            { if (!js_numero(&j, &v)) break; S->target = (int32_t)v; }
        else if (CHAVE("population"))        { if (!js_bool(&j, &b)) break; S->want_pop = b; }
        else if (CHAVE("module_concurrence")){ if (!js_bool(&j, &b)) break; S->want_conc_mod = b; }
        else if (CHAVE("full_concurrence"))  { if (!js_bool(&j, &b)) break; S->want_conc_full = b; }
        else if (CHAVE("pop_stride"))        { if (!js_numero(&j, &v)) break; S->pop_stride = (int32_t)v; }
        else { js_erro(&j, "chave desconhecida em observables"); break; }
      }
    }
    else if (CHAVE("neutrino")) {
      int p2 = 1;
      const double GRAU = 3.14159265358979323846 / 180.0;
      if (!js_obj_inicio(&j)) break;
      S->nu_ativo = 1;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("n_flavors"))    { if (!js_numero(&j, &v)) break; S->nu.n_flavors = (int32_t)v; }
        else if (CHAVE("theta12_deg"))  { if (!js_numero(&j, &v)) break; S->nu.theta12 = v * GRAU; }
        else if (CHAVE("theta13_deg"))  { if (!js_numero(&j, &v)) break; S->nu.theta13 = v * GRAU; }
        else if (CHAVE("theta23_deg"))  { if (!js_numero(&j, &v)) break; S->nu.theta23 = v * GRAU; }
        else if (CHAVE("delta_cp_deg")) { if (!js_numero(&j, &v)) break; S->nu.delta_cp = v * GRAU; }
        else if (CHAVE("dm21_sq_eV2"))  { if (!js_numero(&j, &v)) break; S->nu.dm21_sq = v; }
        else if (CHAVE("dm31_sq_eV2"))  { if (!js_numero(&j, &v)) break; S->nu.dm31_sq = v; }
        else if (CHAVE("energy_MeV"))   { if (!js_numero(&j, &v)) break; S->nu.energy_eV = v * 1.0e6; }
        else if (CHAVE("antineutrino")) { if (!js_bool(&j, &b)) break; S->nu.antineutrino = b; }
        else if (CHAVE("flavor"))       { if (!js_enum(&j, MAPA_SABOR, &e)) break; S->nu_sabor = (int32_t)e; }
        else if (CHAVE("rho_g_cm3"))    { if (!js_numero(&j, &v)) break; S->nu_rho = v; }
        else if (CHAVE("Ye"))           { if (!js_numero(&j, &v)) break; S->nu_ye = v; }
        else { js_erro(&j, "chave desconhecida em neutrino"); break; }
      }
    }
    else if (CHAVE("scan")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("energy_min_MeV")) { if (!js_numero(&j, &v)) break; S->scan_min = v; }
        else if (CHAVE("energy_max_MeV")) { if (!js_numero(&j, &v)) break; S->scan_max = v; }
        else if (CHAVE("points"))         { if (!js_numero(&j, &v)) break; S->scan_pontos = (int32_t)v; }
        else if (CHAVE("log"))            { if (!js_bool(&j, &b)) break; S->scan_log = b; }
        else if (CHAVE("both_signs"))     { if (!js_bool(&j, &b)) break; S->scan_ambos = b; }
        else { js_erro(&j, "chave desconhecida em scan"); break; }
      }
    }
    else if (CHAVE("matter")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      S->nu_perfil_ativo = 1;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("profile"))        { if (!js_enum(&j, MAPA_PERFIL, &e)) break;
                                            S->nu_perfil.kind = (dae_nu_perfil_kind)e; }
        else if (CHAVE("Ye"))             { if (!js_numero(&j, &v)) break; S->nu_perfil.ye = v; }
        else if (CHAVE("rho_g_cm3"))      { if (!js_numero(&j, &v)) break; S->nu_perfil.rho = v; }
        else if (CHAVE("rho_ref_g_cm3"))  { if (!js_numero(&j, &v)) break; S->nu_perfil.rho_ref = v; }
        else if (CHAVE("r_ref_km"))       { if (!js_numero(&j, &v)) break; S->nu_perfil.r_ref = v; }
        else if (CHAVE("exponent"))       { if (!js_numero(&j, &v)) break; S->nu_perfil.expoente = v; }
        else if (CHAVE("r_start_km"))     { if (!js_numero(&j, &v)) break; S->nu_perfil.r_start = v; }
        else if (CHAVE("tolerance"))      { if (!js_numero(&j, &v)) break; S->nu_tol = v; }
        else if (CHAVE("r_km"))           { if (!le_vetor_num(&j, &S->perf_r, &S->perf_n)) break; }
        else if (CHAVE("rho_table"))      { int32_t q; if (!le_vetor_num(&j, &S->perf_rho, &q)) break;
                                            if (q != S->perf_n) { js_erro(&j, "rho_table tem tamanho diferente de r_km"); break; } }
        else if (CHAVE("ye_table"))       { int32_t q; if (!le_vetor_num(&j, &S->perf_ye, &q)) break;
                                            if (q != S->perf_n) { js_erro(&j, "ye_table tem tamanho diferente de r_km"); break; } }
        else if (CHAVE("noise_rms"))      { if (!js_numero(&j, &v)) break; S->nu_perfil.ruido_rms = v; }
        else if (CHAVE("noise_exponent")) { if (!js_numero(&j, &v)) break; S->nu_perfil.ruido_p = v; }
        else if (CHAVE("noise_n_min"))    { if (!js_numero(&j, &v)) break; S->nu_perfil.ruido_n_min = (int32_t)v; }
        else if (CHAVE("noise_n_max"))    { if (!js_numero(&j, &v)) break; S->nu_perfil.ruido_n_max = (int32_t)v; }
        else if (CHAVE("noise_seed"))     { if (!js_numero(&j, &v)) break; S->nu_perfil.ruido_semente = (uint64_t)v; }
        else { js_erro(&j, "chave desconhecida em matter"); break; }
      }
    }
    else if (CHAVE("magnetic")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      S->numag_ativo = 1;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        if      (CHAVE("nature"))           { if (!js_enum(&j, MAPA_NATUREZA, &e)) break;
                                              S->numag.natureza = (dae_numag_natureza)e; }
        else if (CHAVE("mu_nu_bohr"))       { if (!js_numero(&j, &v)) break; S->numag.mu_re[0] = v; }
        else if (CHAVE("energy_MeV"))       { if (!js_numero(&j, &v)) break; S->numag.sabor.energy_eV = v * 1.0e6; }
        else if (CHAVE("antineutrino"))     { if (!js_bool(&j, &b)) break; S->numag.sabor.antineutrino = b; }
        else if (CHAVE("initial_helicity")) { if (!js_enum(&j, MAPA_HELICIDADE, &e)) break; S->numag_helicidade = (int32_t)e; }
        else if (CHAVE("tolerance"))        { if (!js_numero(&j, &v)) break; S->numag_tol = v; }
        else if (CHAVE("integrator"))       { if (!js_enum(&j, MAPA_INTEGRADOR, &e)) break; S->numag_integrador = e; }
        else if (CHAVE("coherence"))        { if (!js_enum(&j, MAPA_COERENCIA, &e)) break; S->numag_coerencia = e; }
        else if (CHAVE("window_K"))         { if (!js_numero(&j, &v)) break; S->numag_jan.K = v; }
        else if (CHAVE("window_eps_max"))   { if (!js_numero(&j, &v)) break; S->numag_jan.eps_max = v; }
        else if (CHAVE("window_eps_tol"))   { if (!js_numero(&j, &v)) break; S->numag_jan.eps_tol = v; }
        else if (CHAVE("phase_tolerance_rad")) { if (!js_numero(&j, &v)) break; S->numag_jan.tol_fase = v; }
        else if (CHAVE("samples_per_interval")) { if (!js_numero(&j, &v)) break; S->numag_jan.sub = (int32_t)v; }
        else if (CHAVE("collision_rate_eV")) { if (!js_numero(&j, &v)) break; S->numag_gamma_eV = v; }
        else if (CHAVE("n_trajectories"))   { if (!js_numero(&j, &v)) break; S->numag_n_traj = (int32_t)v; }
        else { js_erro(&j, "chave desconhecida em magnetic"); break; }
      }
    }
    else if (CHAVE("trajectory")) {
      int p2 = 1;
      if (!js_obj_inicio(&j)) break;
      S->numag_tem_traj = 1;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p2)) {
        int32_t q = 0;
        if      (CHAVE("transverse_basis")) { if (!js_enum(&j, MAPA_BASE, &e)) break; S->numag_base_declarada = e; }
        else if (CHAVE("s_km"))             { if (!le_vetor_num(&j, &S->traj_s, &S->traj_n)) break; }
        else if (CHAVE("rho_table"))        { if (!le_vetor_num(&j, &S->traj_rho, &q)) break;
                                              if (q != S->traj_n) { js_erro(&j, "trajectory.rho_table tem tamanho diferente de s_km"); break; } }
        else if (CHAVE("ye_table"))         { if (!le_vetor_num(&j, &S->traj_ye, &q)) break;
                                              if (q != S->traj_n) { js_erro(&j, "trajectory.ye_table tem tamanho diferente de s_km"); break; } }
        else if (CHAVE("b1_gauss"))         { if (!le_vetor_num(&j, &S->traj_b1, &q)) break;
                                              if (q != S->traj_n) { js_erro(&j, "trajectory.b1_gauss tem tamanho diferente de s_km"); break; } }
        else if (CHAVE("b2_gauss"))         { if (!le_vetor_num(&j, &S->traj_b2, &q)) break;
                                              if (q != S->traj_n) { js_erro(&j, "trajectory.b2_gauss tem tamanho diferente de s_km"); break; } }
        else if (CHAVE("ynu_table"))        { if (!le_vetor_num(&j, &S->traj_ynu, &q)) break;
                                              if (q != S->traj_n) { js_erro(&j, "trajectory.ynu_table tem tamanho diferente de s_km"); break; } }
        else if (CHAVE("energy_ratio"))     { if (!le_vetor_num(&j, &S->traj_e, &q)) break;
                                              if (q != S->traj_n) { js_erro(&j, "trajectory.energy_ratio tem tamanho diferente de s_km"); break; } }
        else if (CHAVE("collision_rate_table_eV")) { if (!le_vetor_num(&j, &S->traj_gamma, &q)) break;
                                              if (q != S->traj_n) { js_erro(&j, "trajectory.collision_rate_table_eV tem tamanho diferente de s_km"); break; } }
        else { js_erro(&j, "chave desconhecida em trajectory"); break; }
      }
    }
    else if (CHAVE("realizations")) { if (!js_numero(&j, &v)) break; S->realizations = (int32_t)v; }
    else if (CHAVE("trajectories")) {
      int p3 = 1;
      tem_traj = 1;
      if (!js_obj_inicio(&j)) break;
      while (js_obj_chave(&j, k, (int32_t)sizeof(k), &p3)) {
        if      (CHAVE("n_traj"))     { if (!js_numero(&j, &v)) break; S->n_traj = (int32_t)v; }
        else if (CHAVE("gamma_deph")) { if (!js_numero(&j, &v)) break; S->gamma_deph = v; }
        else if (CHAVE("rho_stride")) { if (!js_numero(&j, &v)) break; S->rho_stride = (int32_t)v; }
        else if (CHAVE("output_mode")) {
          if (!js_enum(&j, MAPA_SAIDA, &e)) break;
          S->saida_traj = (dae_saida_traj)e;
        }
        else { js_erro(&j, "chave desconhecida em trajectories"); break; }
      }
    }
    else { js_erro(&j, "chave desconhecida na raiz"); break; }
  }
  if (!j.ok) { dae_spec_free(S); return DAE_ERR_JSON; }

  js_ws(&j);
  if (j.s[j.i] != '\0') { js_erro(&j, "lixo depois do objeto"); dae_spec_free(S); return DAE_ERR_JSON; }

  if (S->format_version < DAE_FORMAT_VERSION_MIN ||
      S->format_version > DAE_FORMAT_VERSION) {
    js_erro(&j, "format_version fora da faixa que este nucleo le");
    dae_spec_free(S);
    return DAE_ERR_VERSION;
  }
  S->gen.seed = S->seed;
  if (S->n_edges > 0) S->gen.kind = DAE_G_EDGELIST;
  /* SEM PADRAO IMPLICITO. O bloco existir e nao dizer o modo e erro, e a
     mensagem diz o que fazer — a escolha e do usuario porque as duas saidas
     resolvem problemas diferentes, e escolher por ele so apareceria meses
     depois, quando o observavel mudasse e os psi nao existissem mais. */
  if (tem_traj && S->saida_traj == DAE_SAIDA_INDEFINIDA) {
    js_erro(&j, "trajectories exige output_mode: \"accumulate_rho\" ou \"archive_psi\"");
    dae_spec_free(S);
    return DAE_ERR_JSON;
  }
  if (S->n_traj < 0 || S->gamma_deph < 0.0 || S->rho_stride <= 0) {
    js_erro(&j, "n_traj e gamma_deph nao podem ser negativos, rho_stride precisa ser positivo");
    dae_spec_free(S);
    return DAE_ERR_JSON;
  }
  if (S->nt <= 0 || S->pop_stride <= 0 || S->realizations <= 0) {
    js_erro(&j, "nt, pop_stride e realizations precisam ser positivos");
    dae_spec_free(S);
    return DAE_ERR_PARAM;
  }
  /* O bloco neutrino e campo novo, e o parser e estrito: campo novo exige
     subir format_version (ver o cabecalho de dae_spec.h). Sem esta checagem,
     um spec de versao 2 com bloco neutrino seria aceito por um binario novo e
     recusado por um antigo — divergencia silenciosa entre versoes. */
  if (S->nu_ativo) {
    if (S->format_version < 3) {
      js_erro(&j, "o bloco neutrino exige format_version 3");
      dae_spec_free(S);
      return DAE_ERR_JSON;
    }
    if (S->nu.n_flavors < 2 || S->nu.n_flavors > DAE_NU_MAX_SABORES) {
      js_erro(&j, "neutrino.n_flavors precisa ser 2 ou 3");
      dae_spec_free(S);
      return DAE_ERR_JSON;
    }
    if (!(S->nu.energy_eV > 0.0)) {
      js_erro(&j, "neutrino.energy_MeV precisa ser positiva");
      dae_spec_free(S);
      return DAE_ERR_JSON;
    }
    if (S->nu_sabor >= S->nu.n_flavors) {
      js_erro(&j, "neutrino.flavor fora do numero de sabores");
      dae_spec_free(S);
      return DAE_ERR_JSON;
    }
    if (S->nu_rho < 0.0 || S->nu_ye < 0.0) {
      js_erro(&j, "neutrino.rho_g_cm3 e Ye nao podem ser negativos");
      dae_spec_free(S);
      return DAE_ERR_JSON;
    }
  }
  if (S->scan_pontos > 0) {
    if (!S->nu_ativo) {
      js_erro(&j, "o bloco scan exige o bloco neutrino");
      dae_spec_free(S); return DAE_ERR_JSON;
    }
    if (S->scan_pontos < 2) {
      js_erro(&j, "scan.points precisa ser >= 2");
      dae_spec_free(S); return DAE_ERR_JSON;
    }
    if (!(S->scan_min > 0.0) || !(S->scan_max > S->scan_min)) {
      js_erro(&j, "scan exige 0 < energy_min_MeV < energy_max_MeV");
      dae_spec_free(S); return DAE_ERR_JSON;
    }
  }
  if (S->nu_perfil_ativo) {
    if (!S->nu_ativo) {
      js_erro(&j, "o bloco matter exige o bloco neutrino");
      dae_spec_free(S); return DAE_ERR_JSON;
    }
    if (!(S->nu_tol > 0.0)) {
      js_erro(&j, "matter.tolerance precisa ser positiva");
      dae_spec_free(S); return DAE_ERR_JSON;
    }
    if (S->nu_perfil.ruido_rms > 0.0 &&
        (S->nu_perfil.ruido_n_min < 1 ||
         S->nu_perfil.ruido_n_max < S->nu_perfil.ruido_n_min)) {
      js_erro(&j, "noise exige 1 <= noise_n_min <= noise_n_max");
      dae_spec_free(S); return DAE_ERR_JSON;
    }
    if (S->nu_perfil.ruido_rms > 0.0) {
      const int32_t nm = S->nu_perfil.ruido_n_max - S->nu_perfil.ruido_n_min + 1;
      S->ruido_amp = (double *)malloc((size_t)nm * sizeof(double));
      S->ruido_fase = (double *)malloc((size_t)nm * sizeof(double));
      if (!S->ruido_amp || !S->ruido_fase) { dae_spec_free(S); return DAE_ERR_ALLOC; }
      dae_nu_ruido_prepara(&S->nu_perfil, S->ruido_amp, S->ruido_fase);
    }
    if (S->nu_perfil.kind == DAE_NU_PERFIL_TABELA) {
      int32_t q;
      if (S->perf_n < 2 || !S->perf_r || !S->perf_rho) {
        js_erro(&j, "profile table exige r_km e rho_table com 2+ pontos");
        dae_spec_free(S); return DAE_ERR_JSON;
      }
      for (q = 1; q < S->perf_n; ++q)
        if (!(S->perf_r[q] > S->perf_r[q - 1])) {
          js_erro(&j, "r_km precisa ser estritamente crescente");
          dae_spec_free(S); return DAE_ERR_JSON;
        }
      S->nu_perfil.r_tab = S->perf_r;
      S->nu_perfil.rho_tab = S->perf_rho;
      S->nu_perfil.ye_tab = S->perf_ye;
      S->nu_perfil.n_tab = S->perf_n;
      if (S->nu_perfil.r_start <= 0.0) S->nu_perfil.r_start = S->perf_r[0];
    }
    if (S->nu_perfil.kind == DAE_NU_PERFIL_POTENCIA &&
        !(S->nu_perfil.r_ref > 0.0 && S->nu_perfil.r_start > 0.0)) {
      js_erro(&j, "power_law exige r_ref_km e r_start_km positivos");
      dae_spec_free(S); return DAE_ERR_JSON;
    }
  } else if (S->nu_ativo) {
    /* Sem bloco `matter`, o perfil e o constante que `neutrino.rho_g_cm3` diz.
       Um caminho so daqui para baixo: dae_run nunca le nu_rho direto. */
    S->nu_perfil.kind = DAE_NU_PERFIL_CONSTANTE;
    S->nu_perfil.rho = S->nu_rho;
    S->nu_perfil.ye = S->nu_ye;
  }
  if (S->numag_ativo || S->numag_tem_traj) {
    int32_t q;
    const char *msg = NULL;
    if (S->format_version < 4)
      msg = "os blocos magnetic e trajectory exigem format_version 4";
    else if (!S->numag_ativo || !S->numag_tem_traj)
      msg = "magnetic e trajectory vem juntos: um sem o outro nao define o problema";
    else if (S->nu_ativo || S->nu_perfil_ativo || S->scan_pontos > 0)
      msg = "magnetic e exclusivo com neutrino, matter e scan: e outra fisica";
    else if (S->n_traj > 0 || S->saida_traj != DAE_SAIDA_INDEFINIDA)
      msg = "magnetic nao aceita o bloco trajectories (defasagem)";
    else if (S->numag.natureza != DAE_NUMAG_DIRAC || S->numag.sabor.antineutrino)
      msg = "magnetic: so neutrino de Dirac nesta versao; Majorana e antineutrino ainda nao sao suportados";
    else if (!(S->numag.mu_re[0] >= 0.0))
      msg = "magnetic.mu_nu_bohr precisa ser >= 0";
    else if (!(S->numag.sabor.energy_eV > 0.0))
      msg = "magnetic.energy_MeV precisa ser positiva";
    else if (!(S->numag_tol > 0.0))
      msg = "magnetic.tolerance precisa ser positiva";
    else if (S->numag_integrador == 1 && S->numag_coerencia < 0)
      msg = "magnetic.integrator windows exige coherence: \"coherent\" ou \"incoherent\" (escolha fisica, sem padrao)";
    else if (S->numag_integrador == 0 && S->numag_coerencia == 1)
      msg = "magnetic.coherence incoherent so existe com integrator windows";
    else if (!(S->numag_gamma_eV >= 0.0))
      msg = "magnetic.collision_rate_eV precisa ser >= 0";
    else if ((S->numag_gamma_eV > 0.0 || S->traj_gamma) && S->numag_integrador != 0)
      msg = "colisoes so com integrator brute_force: a evolucao amortecida DURANTE a janela nao e o modo incoerente";
    else if ((S->numag_gamma_eV > 0.0 || S->traj_gamma || S->numag_n_traj > 0) && S->numag_n_traj < 2)
      msg = "colisoes exigem magnetic.n_trajectories >= 2 (o erro estatistico precisa de variancia)";
    else if (S->numag_integrador == 1 &&
             (!(S->numag_jan.K > 0.0) || !(S->numag_jan.eps_max > 0.0) ||
              !(S->numag_jan.eps_tol > 0.0) || !(S->numag_jan.tol_fase > 0.0) ||
              S->numag_jan.sub < 1))
      msg = "magnetic: window_K, window_eps_max, window_eps_tol, phase_tolerance_rad precisam ser positivos e samples_per_interval >= 1";
    else if (!S->numag_base_declarada)
      msg = "trajectory.transverse_basis e obrigatorio: \"parallel_transported\"";
    else if (S->traj_n < 2 || !S->traj_s || !S->traj_rho || !S->traj_ye ||
             !S->traj_b1 || !S->traj_b2)
      msg = "trajectory exige s_km, rho_table, ye_table, b1_gauss e b2_gauss com 2+ pontos";
    if (!msg) {
      for (q = 1; q < S->traj_n && !msg; ++q)
        if (!(S->traj_s[q] > S->traj_s[q - 1])) msg = "trajectory.s_km precisa ser estritamente crescente";
      for (q = 0; q < S->traj_n && !msg; ++q) {
        if (!(S->traj_rho[q] >= 0.0) || !(S->traj_ye[q] >= 0.0))
          msg = "trajectory: rho_table e ye_table precisam ser >= 0";
        else if (S->traj_b1[q] != S->traj_b1[q] || S->traj_b2[q] != S->traj_b2[q])
          msg = "trajectory: b1_gauss e b2_gauss nao podem ser NaN";
        else if (S->traj_ynu && !(fabs(S->traj_ynu[q]) <= 1.0))
          msg = "trajectory.ynu_table precisa ter |Y_nu| <= 1";
        else if (S->traj_e && !(S->traj_e[q] > 0.0))
          msg = "trajectory.energy_ratio precisa ser positiva";
        else if (S->traj_gamma && !(S->traj_gamma[q] >= 0.0))
          msg = "trajectory.collision_rate_table_eV precisa ser >= 0";
      }
    }
    if (!msg && !(S->t1 > 0.0 && S->traj_s[0] + S->t1 <= S->traj_s[S->traj_n - 1]))
      msg = "time.t1 (km) precisa ser positivo e caber na trajetoria: s_km[0] + t1 <= s_km[fim]";
    if (msg) { js_erro(&j, msg); dae_spec_free(S); return DAE_ERR_JSON; }
    S->numag.n_ativos = 1;
    S->numag_jan.tol = S->numag_tol;
    S->numag_traj.s_km = S->traj_s;
    S->numag_traj.rho = S->traj_rho;
    S->numag_traj.ye = S->traj_ye;
    S->numag_traj.b1 = S->traj_b1;
    S->numag_traj.b2 = S->traj_b2;
    S->numag_traj.razao_energia = S->traj_e;
    S->numag_traj.ynu = S->traj_ynu;
    S->numag_traj.n = S->traj_n;
  }
  return DAE_OK;
}

/* ----------------------------------------------------------- canonico --- */

typedef struct { char *b; int32_t cap, n; } dae_buf;

/* `m` e o tamanho PRETENDIDO por snprintf: sem o clamp, um truncamento vira
   leitura fora do buffer. */
static void buf_emite(dae_buf *B, const char *tmp, int m, int32_t lim)
{
  int32_t i;
  if (m < 0) return;
  if (m > lim) m = lim;
  for (i = 0; i < m; ++i) {
    if (B->n + 1 < B->cap) B->b[B->n] = tmp[i];
    ++B->n;
  }
}

/* Tres helpers tipados em vez de um so adivinhando pelo formato. Escolher
   entre int e double olhando o ULTIMO caractere de fmt erraria em
   ",\"cols\":%d}" — termina em '}', nao em 'd'. */
static void bufs(dae_buf *B, const char *fmt, const char *v)
{ char t[512]; buf_emite(B, t, snprintf(t, sizeof(t), fmt, v), (int32_t)sizeof(t) - 1); }
static void bufi(dae_buf *B, const char *fmt, int32_t v)
{ char t[64]; buf_emite(B, t, snprintf(t, sizeof(t), fmt, (int)v), (int32_t)sizeof(t) - 1); }
static void bufd(dae_buf *B, const char *fmt, double v)
{ char t[64]; buf_emite(B, t, snprintf(t, sizeof(t), fmt, v), (int32_t)sizeof(t) - 1); }

int32_t dae_spec_canonical(const dae_spec *S, char *buf, int32_t cap)
{
  dae_buf B;
  if (!S) return 0;
  B.b = buf; B.cap = cap; B.n = 0;
  if (cap > 0) buf[0] = '\0';

  bufi(&B, "{\"format_version\":%d", S->format_version);
  bufs(&B, ",\"core_hash\":\"%s\"", DAE_CORE_HASH);
  bufd(&B, ",\"seed\":%.0f", (double)S->seed);
  bufs(&B, ",\"graph\":{\"generator\":\"%s\"", nome_de(MAPA_GER, (int)S->gen.kind));
  if (S->n_edges > 0) {
    bufi(&B, ",\"n\":%d", S->gen.n);
    bufs(&B, ",\"edges\":[%s", "");
    { int32_t e;
      for (e = 0; e < S->n_edges; ++e) {
        bufs(&B, "%s", e ? "," : "");
        bufi(&B, "[%d", S->ei[e]);
        bufi(&B, ",%d", S->ej[e]);
        bufd(&B, ",%.17g]", S->w[e]);
      } }
    bufs(&B, "]%s", "");
  } else {
    bufi(&B, ",\"params\":{\"n_par\":%d", S->gen.n_par);
    bufi(&B, ",\"n_perp\":%d", S->gen.n_perp);
    bufi(&B, ",\"seam_shift\":%d", S->gen.seam_shift);
    bufs(&B, ",\"longitudinal_closed\":%s", S->gen.longitudinal_closed ? "true" : "false");
    bufd(&B, ",\"j_par\":%.17g", S->gen.j_par);
    bufd(&B, ",\"j_perp\":%.17g", S->gen.j_perp);
    bufi(&B, ",\"n_modules\":%d", S->gen.n_modules);
    bufd(&B, ",\"ws_p\":%.17g", S->gen.ws_p);
    bufs(&B, ",\"conn_mode\":\"%s\"", nome_de(MAPA_CONN, (int)S->gen.conn_mode));
    bufd(&B, ",\"p_in\":%.17g", S->gen.p_in);
    bufd(&B, ",\"p_out\":%.17g", S->gen.p_out);
    bufi(&B, ",\"n\":%d", S->gen.n);
    bufi(&B, ",\"dim\":%d", S->gen.dim);
    bufi(&B, ",\"rows\":%d", S->gen.rows);
    bufi(&B, ",\"cols\":%d}", S->gen.cols);
  }
  bufs(&B, "}%s", "");
  bufs(&B, ",\"hamiltonian\":{\"kind\":\"%s\"", nome_de(MAPA_HAM, (int)S->ham));
  bufd(&B, ",\"gamma\":%.17g", S->gamma);
  bufs(&B, ",\"normalization\":\"%s\"", nome_de(MAPA_NORM, (int)S->norm));
  bufi(&B, ",\"lanczos_steps\":%d}", S->lanczos_steps);
  if (S->init_site >= 0) bufi(&B, ",\"initial\":{\"site\":%d}", S->init_site);
  else {
    bufs(&B, ",\"initial\":{\"vector\":[%s", "");
    { int32_t i;
      for (i = 0; i < S->init_n; ++i) {
        bufs(&B, "%s", i ? ",[" : "[");
        bufd(&B, "%.17g", S->init_re[i]);
        bufd(&B, ",%.17g]", S->init_im[i]);
      } }
    bufs(&B, "]}%s", "");
  }
  bufd(&B, ",\"time\":{\"t1\":%.17g", S->t1);
  bufi(&B, ",\"nt\":%d}", S->nt);
  bufi(&B, ",\"observables\":{\"target\":%d", S->target);
  bufs(&B, ",\"population\":%s", S->want_pop ? "true" : "false");
  bufs(&B, ",\"module_concurrence\":%s", S->want_conc_mod ? "true" : "false");
  bufs(&B, ",\"full_concurrence\":%s", S->want_conc_full ? "true" : "false");
  bufi(&B, ",\"pop_stride\":%d}", S->pop_stride);
  /* O bloco neutrino no canonico e o ponto do GHOST: e ele que faz o CSV
     carregar o CENARIO FISICO, e nao seis nos com pesos numericos. Sem isto,
     de um resultado ninguem recupera theta13, Dm2, E ou rho. */
  if (S->nu_ativo) {
    const double GRAU = 180.0 / 3.14159265358979323846;
    bufi(&B, ",\"neutrino\":{\"n_flavors\":%d", S->nu.n_flavors);
    bufd(&B, ",\"theta12_deg\":%.17g", S->nu.theta12 * GRAU);
    bufd(&B, ",\"theta13_deg\":%.17g", S->nu.theta13 * GRAU);
    bufd(&B, ",\"theta23_deg\":%.17g", S->nu.theta23 * GRAU);
    bufd(&B, ",\"delta_cp_deg\":%.17g", S->nu.delta_cp * GRAU);
    bufd(&B, ",\"dm21_sq_eV2\":%.17g", S->nu.dm21_sq);
    bufd(&B, ",\"dm31_sq_eV2\":%.17g", S->nu.dm31_sq);
    bufd(&B, ",\"energy_MeV\":%.17g", S->nu.energy_eV / 1.0e6);
    bufs(&B, ",\"antineutrino\":%s", S->nu.antineutrino ? "true" : "false");
    bufs(&B, ",\"flavor\":\"%s\"", nome_de(MAPA_SABOR, (int)S->nu_sabor));
    bufd(&B, ",\"rho_g_cm3\":%.17g", S->nu_rho);
    bufd(&B, ",\"Ye\":%.17g}", S->nu_ye);
  }
  if (S->nu_perfil_ativo) {
    bufs(&B, ",\"matter\":{\"profile\":\"%s\"", nome_de(MAPA_PERFIL, (int)S->nu_perfil.kind));
    bufd(&B, ",\"Ye\":%.17g", S->nu_perfil.ye);
    if (S->nu_perfil.kind == DAE_NU_PERFIL_CONSTANTE) {
      bufd(&B, ",\"rho_g_cm3\":%.17g", S->nu_perfil.rho);
    } else {
      bufd(&B, ",\"rho_ref_g_cm3\":%.17g", S->nu_perfil.rho_ref);
      bufd(&B, ",\"r_ref_km\":%.17g", S->nu_perfil.r_ref);
      bufd(&B, ",\"exponent\":%.17g", S->nu_perfil.expoente);
      bufd(&B, ",\"r_start_km\":%.17g", S->nu_perfil.r_start);
    }
    if (S->nu_perfil.kind == DAE_NU_PERFIL_TABELA) {
      /* A tabela vai INTEIRA para o canonico, e portanto para o cabecalho do
         CSV. E cara — centenas de pontos —, e e o preco da procedencia: um
         resultado cujo perfil de materia ficou num arquivo de fora nao e
         reproduzivel por quem so tem o CSV. Reamostre ANTES de rodar. */
      int32_t q;
      bufd(&B, ",\"r_start_km\":%.17g", S->nu_perfil.r_start);
      bufs(&B, "%s", ",\"r_km\":[");
      for (q = 0; q < S->perf_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->perf_r[q]);
      bufs(&B, "%s", "],\"rho_table\":[");
      for (q = 0; q < S->perf_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->perf_rho[q]);
      bufs(&B, "%s", "]");
      if (S->perf_ye) {
        bufs(&B, "%s", ",\"ye_table\":[");
        for (q = 0; q < S->perf_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->perf_ye[q]);
        bufs(&B, "%s", "]");
      }
    }
    bufd(&B, ",\"tolerance\":%.17g", S->nu_tol);
    if (S->nu_perfil.ruido_rms > 0.0) {
      /* O ruido vai INTEIRO para o canonico — rms, expoente, faixa de modos e
         SEMENTE. E a semente que torna cada realizacao regeneravel, e sem ela
         um resultado sobre um campo aleatorio nao e reproduzivel por ninguem. */
      bufd(&B, ",\"noise_rms\":%.17g", S->nu_perfil.ruido_rms);
      bufd(&B, ",\"noise_exponent\":%.17g", S->nu_perfil.ruido_p);
      bufi(&B, ",\"noise_n_min\":%d", S->nu_perfil.ruido_n_min);
      bufi(&B, ",\"noise_n_max\":%d", S->nu_perfil.ruido_n_max);
      bufd(&B, ",\"noise_seed\":%.17g", (double)S->nu_perfil.ruido_semente);
    }
    bufs(&B, "%s", "}");
  }
  if (S->scan_pontos > 0) {
    bufd(&B, ",\"scan\":{\"energy_min_MeV\":%.17g", S->scan_min);
    bufd(&B, ",\"energy_max_MeV\":%.17g", S->scan_max);
    bufi(&B, ",\"points\":%d", S->scan_pontos);
    bufs(&B, ",\"log\":%s", S->scan_log ? "true" : "false");
    bufs(&B, ",\"both_signs\":%s}", S->scan_ambos ? "true" : "false");
  }
  if (S->numag_ativo) {
    /* As tabelas vão INTEIRAS, como as de matter: o resultado tem de ser
       reproduzível por quem só tem o CSV — inclusive a declaração da base. */
    int32_t q;
    bufs(&B, ",\"magnetic\":{\"nature\":\"%s\"", nome_de(MAPA_NATUREZA, (int)S->numag.natureza));
    bufd(&B, ",\"mu_nu_bohr\":%.17g", S->numag.mu_re[0]);
    bufd(&B, ",\"energy_MeV\":%.17g", S->numag.sabor.energy_eV / 1.0e6);
    bufs(&B, ",\"antineutrino\":%s", S->numag.sabor.antineutrino ? "true" : "false");
    bufs(&B, ",\"initial_helicity\":\"%s\"", nome_de(MAPA_HELICIDADE, (int)S->numag_helicidade));
    bufd(&B, ",\"tolerance\":%.17g", S->numag_tol);
    bufs(&B, ",\"integrator\":\"%s\"", nome_de(MAPA_INTEGRADOR, S->numag_integrador));
    if (S->numag_integrador == 1) {
      bufs(&B, ",\"coherence\":\"%s\"", nome_de(MAPA_COERENCIA, S->numag_coerencia));
      bufd(&B, ",\"window_K\":%.17g", S->numag_jan.K);
      bufd(&B, ",\"window_eps_max\":%.17g", S->numag_jan.eps_max);
      bufd(&B, ",\"window_eps_tol\":%.17g", S->numag_jan.eps_tol);
      bufd(&B, ",\"phase_tolerance_rad\":%.17g", S->numag_jan.tol_fase);
      bufi(&B, ",\"samples_per_interval\":%d", S->numag_jan.sub);
    }
    if (S->numag_n_traj > 0) {
      bufd(&B, ",\"collision_rate_eV\":%.17g", S->numag_gamma_eV);
      bufi(&B, ",\"n_trajectories\":%d", S->numag_n_traj);
    }
    bufs(&B, "%s", "}");
    bufs(&B, ",\"trajectory\":{\"transverse_basis\":\"%s\"", nome_de(MAPA_BASE, S->numag_base_declarada));
    bufs(&B, "%s", ",\"s_km\":[");
    for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_s[q]);
    bufs(&B, "%s", "],\"rho_table\":[");
    for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_rho[q]);
    bufs(&B, "%s", "],\"ye_table\":[");
    for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_ye[q]);
    bufs(&B, "%s", "],\"b1_gauss\":[");
    for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_b1[q]);
    bufs(&B, "%s", "],\"b2_gauss\":[");
    for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_b2[q]);
    bufs(&B, "%s", "]");
    /* so quando existe: um spec sem fundo de neutrinos nao carrega o campo
       no canonico */
    if (S->traj_ynu) {
      bufs(&B, "%s", ",\"ynu_table\":[");
      for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_ynu[q]);
      bufs(&B, "%s", "]");
    }
    if (S->traj_e) {
      bufs(&B, "%s", ",\"energy_ratio\":[");
      for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_e[q]);
      bufs(&B, "%s", "]");
    }
    if (S->traj_gamma) {
      bufs(&B, "%s", ",\"collision_rate_table_eV\":[");
      for (q = 0; q < S->traj_n; ++q) bufd(&B, q ? ",%.17g" : "%.17g", S->traj_gamma[q]);
      bufs(&B, "%s", "]");
    }
    bufs(&B, "%s", "}");
  }
  bufi(&B, ",\"realizations\":%d", S->realizations);
  /* O bloco so entra no canonico quando existe de fato: um spec unitario nao
     deve carregar campos de trajetoria zerados, que pareceriam escolha. */
  if (S->n_traj > 0 || S->saida_traj != DAE_SAIDA_INDEFINIDA) {
    bufi(&B, ",\"trajectories\":{\"n_traj\":%d", S->n_traj);
    bufd(&B, ",\"gamma_deph\":%.17g", S->gamma_deph);
    bufi(&B, ",\"rho_stride\":%d", S->rho_stride);
    bufs(&B, ",\"output_mode\":\"%s\"}", nome_de(MAPA_SAIDA, (int)S->saida_traj));
  }
  bufs(&B, "%s", "}");

  if (cap > 0) buf[B.n < cap ? B.n : cap - 1] = '\0';
  return B.n;
}
