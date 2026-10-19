#include "bridge.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define DAE_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define DAE_EXPORT
#endif

/* ARMADILHA DA FRONTEIRA DE MEMORIA, escrita aqui porque e aqui que ela nasce.
 *
 * Todo buffer devolvido por estas funcoes vive no heap do WASM. Qualquer malloc
 * deste lado que faca o heap crescer troca o ArrayBuffer inteiro, e QUALQUER
 * view que o JS tenha guardado passa a apontar para memoria detached — leitura
 * silenciosa de lixo, ou zeros. Isso funciona nos testes pequenos e falha
 * quando o usuario aumenta N, que e a pior hora de descobrir.
 *
 * O contrato e: o lado JS refaz a view a cada leitura, a partir de
 * wasmMemory.buffer. Ver wasm/daedalus.mjs e wasm/teste_memoria.mjs. */

typedef struct {
  dae_spec     S;
  dae_error    erro_json;
  dae_graph    G;
  dae_csr      H;
  dae_cheb     W;
  dae_obs      O;
  dae_obs_cfg  cfg;
  dae_metrics  M;
  dae_series   R;          /* mesmo tipo do nucleo: dae_csv escreve direto dele */
  double      *psire, *psiim;
  float       *popf;       /* n: f32, so para a textura WebGL */
  float       *xy_esp;     /* layout espectral, calculado sob demanda */
  char        *texto;      /* buffer de saida para spec canonico e CSV */
  int32_t      texto_cap;
  int32_t      cursor;
  int          tem_spec, tem_grafo, tem_prop;
  dae_status   ultimo;
} dae_sessao;

DAE_EXPORT int32_t dae_ws_ncol(void)       { return (int32_t)DAE_S_NCOL; }
DAE_EXPORT const char *dae_ws_hash(void)   { return DAE_CORE_HASH; }
DAE_EXPORT const char *dae_ws_versao(void) { return DAE_VERSION; }

DAE_EXPORT dae_sessao *dae_ws_nova(void)
{
  return (dae_sessao *)calloc(1, sizeof(dae_sessao));
}

static void solta_prop(dae_sessao *s)
{
  if (!s->tem_prop) return;
  dae_obs_free(&s->O);
  dae_cheb_free(&s->W);
  dae_csr_free(&s->H);
  dae_series_free(&s->R);
  free(s->psire); free(s->psiim); free(s->popf); free(s->xy_esp);
  s->psire = NULL; s->psiim = NULL; s->popf = NULL; s->xy_esp = NULL;
  s->tem_prop = 0;
}

DAE_EXPORT void dae_ws_libera(dae_sessao *s)
{
  if (!s) return;
  solta_prop(s);
  if (s->tem_grafo) dae_graph_free(&s->G);
  if (s->tem_spec) dae_spec_free(&s->S);
  free(s->texto);
  free(s);
}

/* UM caminho de entrada: o texto do spec.json, interpretado por dae_spec.c —
 * o mesmo parser que o .cpp exportado usa sobre o mesmo texto. A ponte nao tem
 * tabela de parametros propria, e por isso nao pode divergir dos padroes do
 * nucleo. */
DAE_EXPORT int32_t dae_ws_spec(dae_sessao *s, const char *json)
{
  dae_metrics_cfg mc;
  dae_status st;
  int32_t i, npop;

  if (!s || !json) return DAE_ERR_PARAM;
  solta_prop(s);
  if (s->tem_grafo) { dae_graph_free(&s->G); s->tem_grafo = 0; }
  if (s->tem_spec)  { dae_spec_free(&s->S); s->tem_spec = 0; }

  st = dae_spec_parse(&s->S, json, &s->erro_json);
  s->ultimo = st;
  if (st != DAE_OK) return st;
  s->tem_spec = 1;

  if (s->S.n_edges > 0)
    st = dae_graph_from_edges(&s->G, s->S.gen.n, s->S.ei, s->S.ej, s->S.w, s->S.n_edges);
  else
    st = dae_graph_build(&s->G, &s->S.gen);
  s->ultimo = st;
  if (st != DAE_OK) return st;
  s->tem_grafo = 1;

  dae_metrics_cfg_default(&mc);
  dae_metrics_compute(&s->G, &mc, &s->M);

  st = dae_hamiltonian(&s->H, &s->G.A, s->S.ham, s->S.gamma, s->S.norm,
                       s->S.lanczos_steps, &s->R.info.escala);
  if (st != DAE_OK) { s->ultimo = st; return st; }
  st = dae_cheb_init(&s->W, &s->H, s->S.lanczos_steps);
  if (st != DAE_OK) { dae_csr_free(&s->H); s->ultimo = st; return st; }

  s->R.info.lo = s->W.lo; s->R.info.hi = s->W.hi;
  s->R.info.a = s->W.a;   s->R.info.b = s->W.b;
  s->R.info.dt = s->S.t1 / (double)s->S.nt;
  s->R.info.alpha = s->W.a * s->R.info.dt;
  s->R.info.lanczos_used = s->W.lanczos_used;
  s->R.info.fingerprint = dae_graph_fingerprint(&s->G);

  s->cfg.module_of = s->G.module_of;
  s->cfg.nmod = s->G.nmod;
  s->cfg.target = (s->S.target < s->G.n) ? s->S.target : -1;
  s->cfg.want_pop = 1;                      /* a interface sempre desenha */
  s->cfg.want_conc_mod = s->S.want_conc_mod;
  s->cfg.want_conc_full = 0;
  st = dae_obs_alloc(&s->O, &s->cfg, s->G.n);
  if (st != DAE_OK) { dae_cheb_free(&s->W); dae_csr_free(&s->H); s->ultimo = st; return st; }

  npop = 0;
  s->cursor = 0;
  s->R.nt = s->S.nt; s->R.n = s->G.n; s->R.nmod = s->G.nmod; s->R.npop = npop;
  s->R.t    = (double *)calloc((size_t)s->S.nt, sizeof(double));
  s->R.scal = (double *)calloc((size_t)s->S.nt * DAE_S_NCOL, sizeof(double));
  s->R.pmod = (double *)calloc((size_t)s->S.nt * (size_t)s->G.nmod, sizeof(double));
  s->R.psi_re = (double *)calloc((size_t)s->G.n, sizeof(double));
  s->R.psi_im = (double *)calloc((size_t)s->G.n, sizeof(double));
  s->psire = (double *)calloc((size_t)s->G.n, sizeof(double));
  s->psiim = (double *)calloc((size_t)s->G.n, sizeof(double));
  s->popf  = (float  *)calloc((size_t)s->G.n, sizeof(float));
  if (!s->R.t || !s->R.scal || !s->R.pmod || !s->R.psi_re || !s->R.psi_im ||
      !s->psire || !s->psiim || !s->popf) {
    solta_prop(s); s->ultimo = DAE_ERR_ALLOC; return DAE_ERR_ALLOC;
  }
  for (i = 0; i < s->S.nt * DAE_S_NCOL; ++i) s->R.scal[i] = NAN;

  if (s->S.init_site >= 0) {
    if (s->S.init_site >= s->G.n) { solta_prop(s); return DAE_ERR_PARAM; }
    s->psire[s->S.init_site] = 1.0;
  } else {
    double nrm = 0.0;
    if (s->S.init_n != s->G.n) { solta_prop(s); return DAE_ERR_PARAM; }
    for (i = 0; i < s->G.n; ++i) {
      s->psire[i] = s->S.init_re[i]; s->psiim[i] = s->S.init_im[i];
      nrm += s->psire[i] * s->psire[i] + s->psiim[i] * s->psiim[i];
    }
    nrm = sqrt(nrm);
    if (!(nrm > 0.0)) { solta_prop(s); return DAE_ERR_PARAM; }
    for (i = 0; i < s->G.n; ++i) { s->psire[i] /= nrm; s->psiim[i] /= nrm; }
  }
  s->tem_prop = 1;
  s->R.info.k_used = 0;

  dae_obs_eval(s->psire, s->psiim, s->G.n, &s->cfg, &s->O);
  for (i = 0; i < s->G.n; ++i) s->popf[i] = (float)s->O.pop[i];
  s->ultimo = DAE_OK;
  return DAE_OK;
}

/* Valida um spec.json SEM tocar na sessao. Existe para a reimportacao: o CSV
 * que volta do cluster carrega o spec canonico no cabecalho, e ele passa pelo
 * MESMO parser estrito que tudo o mais. Sem isso, um CSV de uma versao
 * anterior, ou com colunas em outra ordem, plotaria e plotaria errado. */
DAE_EXPORT int32_t dae_ws_valida(dae_sessao *s, const char *json)
{
  dae_spec S;
  dae_status st;
  if (!s || !json) return DAE_ERR_PARAM;
  st = dae_spec_parse(&S, json, &s->erro_json);
  if (st == DAE_OK) dae_spec_free(&S);
  return st;
}

DAE_EXPORT int32_t dae_ws_json_linha(const dae_sessao *s)  { return s ? s->erro_json.line : 0; }
DAE_EXPORT int32_t dae_ws_json_coluna(const dae_sessao *s) { return s ? s->erro_json.col : 0; }
DAE_EXPORT const char *dae_ws_json_msg(const dae_sessao *s) { return s ? s->erro_json.msg : ""; }

static char *garante_texto(dae_sessao *s, int32_t precisa)
{
  if (s->texto_cap < precisa) {
    char *n = (char *)realloc(s->texto, (size_t)precisa);
    if (!n) return NULL;
    s->texto = n; s->texto_cap = precisa;
  }
  return s->texto;
}

/* Devolve ponteiro para o JSON canonico (o mesmo que vai no CSV e no .cpp). */
DAE_EXPORT const char *dae_ws_spec_canonico(dae_sessao *s)
{
  int32_t precisa;
  if (!s || !s->tem_spec) return "";
  precisa = dae_spec_canonical(&s->S, NULL, 0) + 1;
  if (!garante_texto(s, precisa)) return "";
  dae_spec_canonical(&s->S, s->texto, precisa);
  return s->texto;
}

/* O CSV do navegador sai do MESMO escritor que o nativo e o .cpp exportado. */
DAE_EXPORT const char *dae_ws_csv(dae_sessao *s, int32_t incluir_estado)
{
  int32_t precisa;
  if (!s || !s->tem_prop) return "";
  precisa = dae_csv(&s->S, &s->G, &s->M, &s->R, incluir_estado, NULL, 0) + 1;
  if (!garante_texto(s, precisa)) return "";
  dae_csv(&s->S, &s->G, &s->M, &s->R, incluir_estado, s->texto, precisa);
  return s->texto;
}

/* Volta ao estado inicial SEM reconstruir o grafo.
 *
 * Existe porque "propagar" de novo tem de recomecar, e nao continuar de onde
 * parou. Sem isto, a segunda propagacao encontraria o cursor no fim, nao
 * avancaria nada, e ainda assim emitiria quadros — todos do mesmo estado final.
 * Reconstruir o grafo resolveria tambem, mas pagaria Lanczos e metricas de
 * novo por nada. */
DAE_EXPORT int32_t dae_ws_reiniciar(dae_sessao *s)
{
  int32_t i;
  if (!s || !s->tem_prop) return DAE_ERR_PARAM;
  for (i = 0; i < s->G.n; ++i) { s->psire[i] = 0.0; s->psiim[i] = 0.0; }
  if (s->S.init_site >= 0) {
    s->psire[s->S.init_site] = 1.0;
  } else {
    double nrm = 0.0;
    for (i = 0; i < s->G.n; ++i) {
      s->psire[i] = s->S.init_re[i]; s->psiim[i] = s->S.init_im[i];
      nrm += s->psire[i] * s->psire[i] + s->psiim[i] * s->psiim[i];
    }
    nrm = sqrt(nrm);
    if (!(nrm > 0.0)) return DAE_ERR_PARAM;
    for (i = 0; i < s->G.n; ++i) { s->psire[i] /= nrm; s->psiim[i] /= nrm; }
  }
  for (i = 0; i < s->S.nt * DAE_S_NCOL; ++i) s->R.scal[i] = NAN;
  for (i = 0; i < s->S.nt * s->G.nmod; ++i) s->R.pmod[i] = 0.0;
  s->cursor = 0;
  s->R.info.k_used = 0;
  dae_obs_eval(s->psire, s->psiim, s->G.n, &s->cfg, &s->O);
  for (i = 0; i < s->G.n; ++i) s->popf[i] = (float)s->O.pop[i];
  return DAE_OK;
}

/* Avanca ate `quantos` passos e devolve quantos de fato avancou. O cancelamento
 * mora do lado JS: quem chama decide se pede o proximo bloco. */
DAE_EXPORT int32_t dae_ws_avanca(dae_sessao *s, int32_t quantos)
{
  int32_t feitos = 0, i;
  dae_cheb_info ci;
  if (!s || !s->tem_prop) return 0;
  ci.k_used = s->R.info.k_used;
  while (feitos < quantos && s->cursor < s->R.nt) {
    const int32_t cursor = s->cursor;
    s->ultimo = dae_cheb_step(&s->W, &s->H, s->R.info.dt, s->psire, s->psiim, &ci);
    if (s->ultimo != DAE_OK) return -1;
    dae_obs_eval(s->psire, s->psiim, s->G.n, &s->cfg, &s->O);
    s->R.t[cursor] = (double)(cursor + 1) * s->R.info.dt;
    { double *r = s->R.scal + (size_t)cursor * DAE_S_NCOL;
      r[DAE_S_NORMA]  = s->O.norm;
      r[DAE_S_IPR]    = s->O.ipr;
      r[DAE_S_COH_L1] = s->O.coh_l1;
      r[DAE_S_P_ALVO] = s->O.p_target; }
    for (i = 0; i < s->G.nmod; ++i)
      s->R.pmod[(size_t)cursor * (size_t)s->G.nmod + (size_t)i] = s->O.pmod[i];
    s->cursor = cursor + 1;
    ++feitos;
  }
  s->R.info.k_used = ci.k_used;
  for (i = 0; i < s->G.n; ++i) s->popf[i] = (float)s->O.pop[i];
  memcpy(s->R.psi_re, s->psire, (size_t)s->G.n * sizeof(double));
  memcpy(s->R.psi_im, s->psiim, (size_t)s->G.n * sizeof(double));
  return feitos;
}

DAE_EXPORT int32_t dae_ws_cursor(const dae_sessao *s) { return s ? s->cursor : 0; }
DAE_EXPORT int32_t dae_ws_nt(const dae_sessao *s)     { return s ? s->R.nt : 0; }
DAE_EXPORT int32_t dae_ws_erro(const dae_sessao *s)   { return s ? (int32_t)s->ultimo : DAE_ERR_PARAM; }
DAE_EXPORT const char *dae_ws_erro_texto(int32_t c)   { return dae_strerror((dae_status)c); }

DAE_EXPORT int32_t dae_ws_n(const dae_sessao *s)      { return (s && s->tem_grafo) ? s->G.n : 0; }
DAE_EXPORT int32_t dae_ws_nnz(const dae_sessao *s)    { return (s && s->tem_grafo) ? s->G.A.nnz : 0; }
DAE_EXPORT int32_t dae_ws_nmod(const dae_sessao *s)   { return (s && s->tem_grafo) ? s->G.nmod : 0; }
DAE_EXPORT int32_t dae_ws_geom(const dae_sessao *s)
{ return (s && s->tem_grafo) ? s->G.geom : 0; }

/* Layout espectral, calculado SOB DEMANDA: e um Lanczos deflacionado a mais, e
 * so faz sentido pagar por ele quando a vista espectral for de fato pedida.
 * Devolve 0 se nao deu (grafo grande demais, ou falha). */
DAE_EXPORT const float *dae_ws_espectral(dae_sessao *s)
{
  if (!s || !s->tem_grafo) return 0;
  if (s->xy_esp) return s->xy_esp;
  if (s->G.n < 2 || s->G.n > 20000) return 0;
  s->xy_esp = (float *)calloc(2u * (size_t)s->G.n, sizeof(float));
  if (!s->xy_esp) return 0;
  if (dae_layout_espectral(&s->G, s->xy_esp) != DAE_OK) {
    free(s->xy_esp); s->xy_esp = 0; return 0;
  }
  return s->xy_esp;
}

DAE_EXPORT int32_t dae_ws_npar(const dae_sessao *s)   { return (s && s->tem_grafo) ? s->G.n_par : 0; }
DAE_EXPORT int32_t dae_ws_nperp(const dae_sessao *s)  { return (s && s->tem_grafo) ? s->G.n_perp : 0; }
DAE_EXPORT int32_t dae_ws_descartadas(const dae_sessao *s) { return (s && s->tem_grafo) ? s->G.n_dropped : 0; }
DAE_EXPORT int32_t dae_ws_religa_falhas(const dae_sessao *s) { return (s && s->tem_grafo) ? s->G.n_rewire_failed : 0; }

DAE_EXPORT const int32_t *dae_ws_rowptr(const dae_sessao *s)  { return s->G.A.rowptr; }
DAE_EXPORT const int32_t *dae_ws_colind(const dae_sessao *s)  { return s->G.A.colind; }
/* Os pesos, para o emissor montar a lista de arestas explicita dos oraculos
   .wl e .py. Eles RECEBEM o grafo em vez de regenera-lo: reimplementar o
   gerador em Wolfram seria uma segunda implementacao do que a amalgamacao
   existe para manter unica, e um oraculo que reproduz o mesmo engano do gerador
   nao e oraculo. Em compensacao validam so o PROPAGADOR — quem valida o emissor
   e o .cpp, que regenera tudo a partir do spec.json. */
DAE_EXPORT const double  *dae_ws_valores(const dae_sessao *s) { return s->G.A.val; }
DAE_EXPORT const int32_t *dae_ws_modulos(const dae_sessao *s) { return s->G.module_of; }
DAE_EXPORT const float   *dae_ws_xy(const dae_sessao *s)      { return s->G.xy; }
DAE_EXPORT const float   *dae_ws_pop(const dae_sessao *s)     { return s->popf; }
DAE_EXPORT const double  *dae_ws_scal(const dae_sessao *s)    { return s->R.scal; }
DAE_EXPORT const double  *dae_ws_pmod(const dae_sessao *s)    { return s->R.pmod; }
DAE_EXPORT const double  *dae_ws_conc_mod(const dae_sessao *s) { return s->O.conc_mod; }

DAE_EXPORT double  dae_ws_escala(const dae_sessao *s) { return s ? s->R.info.escala : 1.0; }
DAE_EXPORT double  dae_ws_dt(const dae_sessao *s)     { return s ? s->R.info.dt : 0.0; }
DAE_EXPORT double  dae_ws_alpha(const dae_sessao *s)  { return s ? s->R.info.alpha : 0.0; }
DAE_EXPORT double  dae_ws_cheb_a(const dae_sessao *s) { return s ? s->R.info.a : 0.0; }
DAE_EXPORT double  dae_ws_fingerprint_lo(const dae_sessao *s)
{ return s ? (double)(uint32_t)(s->R.info.fingerprint & 0xFFFFFFFFULL) : 0.0; }
DAE_EXPORT double  dae_ws_fingerprint_hi(const dae_sessao *s)
{ return s ? (double)(uint32_t)(s->R.info.fingerprint >> 32) : 0.0; }

DAE_EXPORT double  dae_ws_lambda2(const dae_sessao *s)     { return s->M.lambda2; }
DAE_EXPORT double  dae_ws_lambda2_res(const dae_sessao *s) { return s->M.lambda2_residual; }
DAE_EXPORT int32_t dae_ws_lambda2_ok(const dae_sessao *s)  { return s->M.lambda2_converged; }
DAE_EXPORT double  dae_ws_Q(const dae_sessao *s)           { return s->M.modularity_Q; }
DAE_EXPORT double  dae_ws_grau(const dae_sessao *s)        { return s->M.mean_degree; }
DAE_EXPORT double  dae_ws_caminho(const dae_sessao *s)     { return s->M.mean_path_len; }
DAE_EXPORT int32_t dae_ws_caminho_exato(const dae_sessao *s) { return s->M.path_len_exact; }
DAE_EXPORT int32_t dae_ws_componentes(const dae_sessao *s) { return s->M.n_components; }
DAE_EXPORT int32_t dae_ws_arestas(const dae_sessao *s)     { return s->M.n_edges; }
