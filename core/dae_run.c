#include "dae_run.h"

#include "dae_cheb.h"
#include "dae_csr.h"
#include "dae_ham.h"
#include "dae_obs.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void dae_series_free(dae_series *R)
{
  if (!R) return;
  free(R->t); free(R->scal); free(R->pmod); free(R->pop); free(R->conc_mod);
  free(R->psi_re); free(R->psi_im);
  free(R->nu_sabor); free(R->nu_massa); free(R->nu_inc);
  memset(R, 0, sizeof(*R));
}

dae_status dae_run(const dae_spec *S, uint64_t seed,
                   dae_graph *G, dae_metrics *M, dae_series *R,
                   dae_progress_fn cb, void *user)
{
  dae_csr H;
  dae_cheb W;
  dae_cheb_info ci;
  dae_obs_cfg cfg;
  dae_obs O;
  double *re = NULL, *im = NULL, dt;
  int32_t i, s, npop = 0;
  dae_status st;
  int nu_merg;
  double nu_ure[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  double nu_uim[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
  int perfil_variavel;
  double passo_km;
  int32_t n_obs;
  double *obs_re = NULL, *obs_im = NULL;

  if (!S || !G || !R) return DAE_ERR_PARAM;
  /* O setor magnético roda por dae_numag_run: outra física, outro CSV.
     Deixar passar aqui rodaria o grafo padrão em silêncio. */
  if (S->numag_ativo) return DAE_ERR_PARAM;
  memset(G, 0, sizeof(*G));
  memset(R, 0, sizeof(*R));
  memset(&H, 0, sizeof(H));
  memset(&W, 0, sizeof(W));
  memset(&O, 0, sizeof(O));
  memset(&ci, 0, sizeof(ci));
  nu_merg = 0;
  perfil_variavel = 0;
  passo_km = 0.0;

  /* --- grafo ---
     Com o bloco `neutrino`, quem monta a adjacencia e a FISICA: PMNS, massas,
     energia e potencial de materia. `graph` e ignorado de proposito — nao ha
     grafo nenhum aqui, e fingir que ha esconderia a fisica em pesos
     numericos. O mergulho 2n x 2n, quando delta_cp exige, acontece dentro
     de dae_nu_adjacencia e o resto do pipeline nao sabe. */
  if (S->nu_ativo) {
    double hre[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
    double him[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
    st = dae_nu_hamiltoniano(&S->nu, S->nu_perfil.rho, S->nu_perfil.ye, hre, him);
    if (st != DAE_OK) return st;
    st = dae_nu_adjacencia(hre, him, S->nu.n_flavors, &G->A, &nu_merg);
    if (st != DAE_OK) return st;
    G->n = G->A.n;
    G->nmod = 1;
    G->module_of = (int32_t *)calloc((size_t)G->n, sizeof(int32_t));
    if (!G->module_of) { dae_csr_free(&G->A); return DAE_ERR_ALLOC; }
    perfil_variavel = (S->nu_perfil.kind != DAE_NU_PERFIL_CONSTANTE);
    passo_km = S->t1 / (double)S->nt;
  } else if (S->n_edges > 0) {
    if (S->gen.n <= 0) return DAE_ERR_PARAM;
    st = dae_graph_from_edges(G, S->gen.n, S->ei, S->ej, S->w, S->n_edges);
  } else {
    dae_gen_params p = S->gen;
    p.seed = seed;
    st = dae_graph_build(G, &p);
  }
  if (st != DAE_OK) return st;
  R->info.fingerprint = dae_graph_fingerprint(G);

  if (M) {
    dae_metrics_cfg mc;
    dae_metrics_cfg_default(&mc);
    st = dae_metrics_compute(G, &mc, M);
    if (st != DAE_OK) goto falhou;
  }

  /* --- hamiltoniano e propagador --- */
  /* Com o bloco neutrino, a normalizacao espectral tem de ficar DESLIGADA.
     Ela e a disciplina de CONVENTIONS 3.2 — comparar topologias a escala de
     energia fixa — e por isso vem ligada por padrao. Mas aqui a escala de
     energia E a fisica: H esta em eV, e reescalonar apagaria a relacao entre
     Dm2/2E e a linha de base em km. O sintoma, se alguem religar, e imediato e
     nao e silencioso: alpha = a*dt explode e o propagador devolve
     DAE_ERR_TOOBIG. */
  if (S->nu_ativo) {
    st = dae_hamiltonian(&H, &G->A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE,
                         0, &R->info.escala);
  } else {
    st = dae_hamiltonian(&H, &G->A, S->ham, S->gamma, S->norm,
                         S->lanczos_steps, &R->info.escala);
  }
  if (st != DAE_OK) goto falhou;
  st = dae_cheb_init(&W, &H, S->lanczos_steps);
  if (st != DAE_OK) goto falhou;

  /* Com o bloco neutrino, t1 esta em QUILOMETROS e H em eV: converte para a
     unidade do propagador (1/eV). E o unico campo que muda de unidade, e o
     canonico registra qual dos dois regimes valeu. */
  dt = (S->nu_ativo ? S->t1 * DAE_NU_T_POR_KM : S->t1) / (double)S->nt;
  R->info.lo = W.lo; R->info.hi = W.hi; R->info.a = W.a; R->info.b = W.b;
  R->info.dt = dt;
  R->info.alpha = W.a * dt;
  R->info.lanczos_used = W.lanczos_used;

  /* --- estado inicial --- */
  re = (double *)calloc((size_t)G->n, sizeof(double));
  im = (double *)calloc((size_t)G->n, sizeof(double));
  if (!re || !im) { st = DAE_ERR_ALLOC; goto falhou; }
  if (S->nu_ativo) {
    dae_nu_inicial(S->nu_sabor, S->nu.n_flavors, nu_merg, re, im);
  } else if (S->init_site >= 0) {
    if (S->init_site >= G->n) { st = DAE_ERR_PARAM; goto falhou; }
    re[S->init_site] = 1.0;
  } else {
    double nrm = 0.0;
    if (S->init_n != G->n) { st = DAE_ERR_PARAM; goto falhou; }
    for (i = 0; i < G->n; ++i) {
      re[i] = S->init_re[i]; im[i] = S->init_im[i];
      nrm += re[i] * re[i] + im[i] * im[i];
    }
    nrm = sqrt(nrm);
    if (!(nrm > 0.0)) { st = DAE_ERR_PARAM; goto falhou; }
    for (i = 0; i < G->n; ++i) { re[i] /= nrm; im[i] /= nrm; }
  }

  /* --- observáveis e buffers ---
     No mergulho o propagador anda em 2n sitios, mas o observavel e do espaco
     FISICO de n sabores: |phi_a|^2 vale |psi_a|^2 / 2, e reportar isso seria
     errado por um fator 2 exato — o tipo de erro que passa despercebido porque
     as curvas tem a forma certa. Projeta antes de medir. */
  n_obs = (S->nu_ativo && nu_merg) ? S->nu.n_flavors : G->n;
  if (S->nu_ativo && nu_merg) {
    obs_re = (double *)calloc((size_t)n_obs, sizeof(double));
    obs_im = (double *)calloc((size_t)n_obs, sizeof(double));
    if (!obs_re || !obs_im) { st = DAE_ERR_ALLOC; goto falhou; }
  }
  cfg.module_of = G->module_of;
  cfg.nmod = G->nmod;
  cfg.target = (S->target < n_obs) ? S->target : -1;
  cfg.want_pop = S->want_pop;
  cfg.want_conc_mod = S->want_conc_mod;
  cfg.want_conc_full = S->want_conc_full;
  st = dae_obs_alloc(&O, &cfg, n_obs);
  if (st != DAE_OK) goto falhou;

  if (S->nu_ativo) {
    const size_t nb = (size_t)S->nt * (size_t)S->nu.n_flavors * sizeof(double);
    st = dae_nu_pmns(&S->nu, nu_ure, nu_uim);
    if (st != DAE_OK) goto falhou;
    R->nu_nf = S->nu.n_flavors;
    R->nu_sabor = (double *)malloc(nb);
    R->nu_massa = (double *)malloc(nb);
    R->nu_inc   = (double *)malloc(nb);
    if (!R->nu_sabor || !R->nu_massa || !R->nu_inc) { st = DAE_ERR_ALLOC; goto falhou; }
  }
  if (S->want_pop) npop = (S->nt + S->pop_stride - 1) / S->pop_stride;
  R->nt = S->nt; R->n = n_obs; R->nmod = G->nmod; R->npop = npop;
  R->psi_re = (double *)malloc((size_t)n_obs * sizeof(double));
  R->psi_im = (double *)malloc((size_t)n_obs * sizeof(double));
  if (!R->psi_re || !R->psi_im) { st = DAE_ERR_ALLOC; goto falhou; }
  R->t    = (double *)malloc((size_t)S->nt * sizeof(double));
  R->scal = (double *)malloc((size_t)S->nt * DAE_S_NCOL * sizeof(double));
  R->pmod = (double *)malloc((size_t)S->nt * (size_t)G->nmod * sizeof(double));
  if (!R->t || !R->scal || !R->pmod) { st = DAE_ERR_ALLOC; goto falhou; }
  if (npop > 0) {
    R->pop = (double *)malloc((size_t)npop * (size_t)n_obs * sizeof(double));
    if (!R->pop) { st = DAE_ERR_ALLOC; goto falhou; }
  }
  if (S->want_conc_mod) {
    R->conc_mod = (double *)malloc((size_t)S->nt * (size_t)G->nmod
                                   * (size_t)G->nmod * sizeof(double));
    if (!R->conc_mod) { st = DAE_ERR_ALLOC; goto falhou; }
  }

  /* --- propagação --- */
  for (s = 0; s < S->nt; ++s) {
    if (perfil_variavel) {
      /* Perfil variavel: H muda dentro do passo, entao quem anda e
         dae_nu_propaga, com o guarda de adiabaticidade decidindo quantos
         sub-passos a fatia precisa. O `dt` fixo do caminho generico nao vale
         aqui — nem o cache de coeficientes de Chebyshev, que pressupoe H fixo
         (CONVENTIONS.md 6.4). */
      const double ra = S->nu_perfil.r_start + (double)s * passo_km;
      const double rb = S->nu_perfil.r_start + (double)(s + 1) * passo_km;
      int32_t sub = 0;
      st = dae_nu_propaga(&S->nu, &S->nu_perfil, ra, rb, S->nu_tol,
                          nu_merg, re, im, &sub);
      if (st != DAE_OK) goto falhou;
      if (sub > R->info.k_used) R->info.k_used = sub;   /* diagnostico */
    } else {
      st = dae_cheb_step(&W, &H, dt, re, im, &ci);
      if (st != DAE_OK) goto falhou;
    }
    if (obs_re) {
      dae_nu_projeta(re, im, S->nu.n_flavors, 1, obs_re, obs_im);
      st = dae_obs_eval(obs_re, obs_im, n_obs, &cfg, &O);
    } else {
      st = dae_obs_eval(re, im, G->n, &cfg, &O);
    }
    if (st != DAE_OK) goto falhou;

    /* O eixo volta para a unidade da ENTRADA. Com o bloco neutrino, t1 veio
       em km; reportar o CSV em 1/eV vazaria a unidade interna do propagador
       para o resultado, e ninguem leria "L = 4.9e12 km" como um numero certo
       em outra unidade — leria como um bug. */
    R->t[s] = perfil_variavel
      ? S->nu_perfil.r_start + (double)(s + 1) * passo_km
      : (double)(s + 1) * (S->nu_ativo ? S->t1 / (double)S->nt : dt);
    { double *r = R->scal + (size_t)s * DAE_S_NCOL;
      r[DAE_S_NORMA]  = O.norm;
      r[DAE_S_IPR]    = O.ipr;
      r[DAE_S_COH_L1] = O.coh_l1;
      r[DAE_S_P_ALVO] = O.p_target; }
    for (i = 0; i < G->nmod; ++i)
      R->pmod[(size_t)s * (size_t)G->nmod + (size_t)i] = O.pmod[i];
    if (R->conc_mod) {
      const size_t m2 = (size_t)G->nmod * (size_t)G->nmod;
      memcpy(R->conc_mod + (size_t)s * m2, O.conc_mod, m2 * sizeof(double));
    }
    if (R->nu_sabor) {
      const int32_t nf = S->nu.n_flavors;
      const double *pr = obs_re ? obs_re : re;
      const double *pi = obs_im ? obs_im : im;
      int32_t aa, mm;
      for (aa = 0; aa < nf; ++aa)
        R->nu_sabor[(size_t)s * (size_t)nf + (size_t)aa] =
          pr[aa] * pr[aa] + pi[aa] * pi[aa];
      /* |<nu_m|psi>|^2 = |sum_a conj(U_am) psi_a|^2 */
      for (mm = 0; mm < nf; ++mm) {
        double sr = 0.0, si = 0.0;
        for (aa = 0; aa < nf; ++aa) {
          const double cr = nu_ure[aa * nf + mm], cim = -nu_uim[aa * nf + mm];
          sr += cr * pr[aa] - cim * pi[aa];
          si += cr * pi[aa] + cim * pr[aa];
        }
        R->nu_massa[(size_t)s * (size_t)nf + (size_t)mm] = sr * sr + si * si;
      }
      for (aa = 0; aa < nf; ++aa) {
        double acc = 0.0;
        for (mm = 0; mm < nf; ++mm) {
          const double ur = nu_ure[aa * nf + mm], ui = nu_uim[aa * nf + mm];
          acc += (ur * ur + ui * ui) * R->nu_massa[(size_t)s * (size_t)nf + (size_t)mm];
        }
        R->nu_inc[(size_t)s * (size_t)nf + (size_t)aa] = acc;
      }
    }
    if (R->pop && (s % S->pop_stride) == 0) {
      const int32_t k = s / S->pop_stride;
      memcpy(R->pop + (size_t)k * (size_t)n_obs, O.pop,
             (size_t)n_obs * sizeof(double));
    }
    if (cb && ((s & 31) == 0 || s == S->nt - 1)) {
      if (!cb(user, s + 1, S->nt)) { st = DAE_ERR_CANCELLED; goto falhou; }
    }
  }
  /* Com perfil variavel, k_used ja e o maximo de subpassos por fatia,
     acumulado no laco; `ci` so e preenchido no caminho de H fixo, e
     sobrescrever o maximo com ele leria `ci` nao inicializado. Ver t108. */
  if (!perfil_variavel) R->info.k_used = ci.k_used;
  if (obs_re) {
    dae_nu_projeta(re, im, S->nu.n_flavors, 1, obs_re, obs_im);
    memcpy(R->psi_re, obs_re, (size_t)n_obs * sizeof(double));
    memcpy(R->psi_im, obs_im, (size_t)n_obs * sizeof(double));
  } else {
    memcpy(R->psi_re, re, (size_t)n_obs * sizeof(double));
    memcpy(R->psi_im, im, (size_t)n_obs * sizeof(double));
  }

  free(re); free(im); free(obs_re); free(obs_im);
  dae_obs_free(&O); dae_cheb_free(&W); dae_csr_free(&H);
  return DAE_OK;

falhou:
  free(re); free(im); free(obs_re); free(obs_im);
  dae_obs_free(&O); dae_cheb_free(&W); dae_csr_free(&H);
  dae_series_free(R);
  dae_graph_free(G);
  return st;
}

/* ------------------------------------------------------- setor magnético -- */

void dae_numag_serie_free(dae_numag_serie *R)
{
  if (!R) return;
  free(R->s_km); free(R->p_l); free(R->p_r); free(R->norma); free(R->subdiv);
  free(R->p_r_erro); free(R->coh_re); free(R->coh_im);
  memset(R, 0, sizeof(*R));
}

dae_status dae_numag_run(const dae_spec *S, dae_numag_serie *R,
                         dae_progress_fn cb, void *user)
{
  double re[4], im[4], pre[2], pim[2], pop[2];
  dae_numag_janela_aberta ab;
  dae_status st = DAE_OK;
  int32_t s, sub;
  double passo, s0;
  /* Cópia local: o callback recebe `void *` não-const, e o spec é const. */
  dae_numag_trajetoria traj;

  if (!S || !R) return DAE_ERR_PARAM;
  memset(R, 0, sizeof(*R));
  if (!S->numag_ativo || S->traj_n < 2) return DAE_ERR_PARAM;
  traj = S->numag_traj;

  R->nt = S->nt;
  R->s_km = (double *)malloc((size_t)S->nt * sizeof(double));
  R->p_l = (double *)malloc((size_t)S->nt * sizeof(double));
  R->p_r = (double *)malloc((size_t)S->nt * sizeof(double));
  R->norma = (double *)malloc((size_t)S->nt * sizeof(double));
  R->subdiv = (int32_t *)malloc((size_t)S->nt * sizeof(int32_t));
  if (!R->s_km || !R->p_l || !R->p_r || !R->norma || !R->subdiv) {
    dae_numag_serie_free(R);
    return DAE_ERR_ALLOC;
  }

  s0 = S->traj_s[0];
  passo = S->t1 / (double)S->nt;

  if (S->numag_n_traj > 0) {
    /* ENSEMBLE com colisões: taxa de salto do desdobramento simétrico
       gamma = Gamma T / 2, por nó se houver tabela. */
    dae_numag_colisao col;
    double *gtab = NULL, *var = NULL, pr0[2], pi0[2];
    int32_t q;
    col.gamma0_km = dae_numag_salto_por_km(S->numag_gamma_eV);
    col.gamma_km = NULL;
    if (S->traj_gamma) {
      gtab = (double *)malloc((size_t)S->traj_n * sizeof(double));
      if (!gtab) { dae_numag_serie_free(R); return DAE_ERR_ALLOC; }
      for (q = 0; q < S->traj_n; ++q) gtab[q] = dae_numag_salto_por_km(S->traj_gamma[q]);
      col.gamma_km = gtab;
    }
    R->p_r_erro = (double *)malloc((size_t)S->nt * sizeof(double));
    R->coh_re = (double *)malloc((size_t)S->nt * sizeof(double));
    R->coh_im = (double *)malloc((size_t)S->nt * sizeof(double));
    var = (double *)malloc((size_t)S->nt * sizeof(double));
    if (!R->p_r_erro || !R->coh_re || !R->coh_im || !var) {
      free(gtab); free(var); dae_numag_serie_free(R); return DAE_ERR_ALLOC;
    }
    pr0[0] = S->numag_helicidade ? 0.0 : 1.0; pr0[1] = 1.0 - pr0[0];
    pi0[0] = pi0[1] = 0.0;
    st = dae_numag_ensemble(&S->numag, &traj, s0, s0 + S->t1, &col, S->numag_tol,
                            S->numag_n_traj, S->seed, pr0, pi0, S->nt,
                            R->p_r, var, R->coh_re, R->coh_im, &R->saltos_medio);
    free(gtab);
    if (st != DAE_OK) { free(var); dae_numag_serie_free(R); return st; }
    for (s = 0; s < S->nt; ++s) {
      R->s_km[s] = s0 + (double)(s + 1) * passo;
      R->p_l[s] = 1.0 - R->p_r[s];      /* a defasagem conserva a norma */
      R->norma[s] = 1.0;
      R->subdiv[s] = 0;
      R->p_r_erro[s] = sqrt(var[s] / (double)S->numag_n_traj);
    }
    free(var);
    R->n_traj = S->numag_n_traj;
    R->salto_por_km = col.gamma0_km;
    R->tem_estado = 0;
    return DAE_OK;
  }

  dae_nu_inicial(S->numag_helicidade, 2, 1, re, im);
  {
    /* Modo incoerente: populações dos autoestados de H~ em s0, a partir da
       helicidade inicial. |<+~|R^dag e_h>|^2 = cos^2 th (h = L) ou sin^2 th. */
    double th, fi;
    dae_numag_angulos(&S->numag, &traj, s0, 0, &th, &fi);
    pop[0] = S->numag_helicidade ? sin(th) * sin(th) : cos(th) * cos(th);
    pop[1] = 1.0 - pop[0];
    memset(&ab, 0, sizeof ab);
  }
  for (s = 0; s < S->nt; ++s) {
    const double a = s0 + (double)s * passo, b = s0 + (double)(s + 1) * passo;
    sub = 0;
    if (S->numag_integrador == 0) {
      st = dae_numag_propaga_tabela(&S->numag, &traj, a, b, S->numag_tol, re, im, &sub);
    } else if (S->numag_coerencia == 0) {
      const int32_t antes = R->janelas.subdiv_max;
      R->janelas.subdiv_max = 0;
      st = dae_numag_propaga_janelas(&S->numag, &traj, a, b, &S->numag_jan,
                                     re, im, &R->janelas);
      sub = R->janelas.subdiv_max;
      if (antes > R->janelas.subdiv_max) R->janelas.subdiv_max = antes;
    } else {
      /* Incoerente: uma janela que atravessa a fronteira entre fatias fica
         ABERTA (`ab`), com o operador parcial, em vez de ser fechada por
         probabilidade no meio da ressonância (ver dae_numag.h). */
      const int32_t antes = R->janelas.subdiv_max;
      R->janelas.subdiv_max = 0;
      st = dae_numag_propaga_janelas_inc(&S->numag, &traj, a, b, &S->numag_jan,
                                         pop, &ab, &R->janelas);
      sub = R->janelas.subdiv_max;
      if (antes > R->janelas.subdiv_max) R->janelas.subdiv_max = antes;
    }
    if (st != DAE_OK) break;
    R->s_km[s] = b;
    if (S->numag_integrador == 1 && S->numag_coerencia == 1) {
      /* populações adiabáticas -> helicidade, INCOERENTE: sem termo cruzado */
      double th, fi;
      dae_numag_angulos(&S->numag, &traj, b, 1, &th, &fi);
      R->p_l[s] = pop[0] * cos(th) * cos(th) + pop[1] * sin(th) * sin(th);
      R->p_r[s] = pop[0] * sin(th) * sin(th) + pop[1] * cos(th) * cos(th);
    } else {
      dae_nu_projeta(re, im, 2, 1, pre, pim);
      R->p_l[s] = pre[0] * pre[0] + pim[0] * pim[0];
      R->p_r[s] = pre[1] * pre[1] + pim[1] * pim[1];
    }
    R->norma[s] = R->p_l[s] + R->p_r[s];
    R->subdiv[s] = sub;
    if (sub > R->subdiv_max) R->subdiv_max = sub;
    if (cb && ((s & 31) == 0 || s == S->nt - 1)) {
      if (!cb(user, s + 1, S->nt)) { st = DAE_ERR_CANCELLED; break; }
    }
  }
  if (st != DAE_OK) { dae_numag_serie_free(R); return st; }
  R->tem_estado = !(S->numag_integrador == 1 && S->numag_coerencia == 1);
  if (R->tem_estado) dae_nu_projeta(re, im, 2, 1, R->psi_re, R->psi_im);
  return DAE_OK;
}
