/* cli.c — `ghost run spec.json`.
 *
 * O mesmo núcleo, o mesmo parser, o mesmo CSV que o navegador e o .cpp
 * exportado. Este binário existe para o teste 7 ter um dos lados da comparação
 * e para rodar um spec.json à mão sem abrir o navegador.
 */
#include "dae.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static char *le_arquivo(const char *caminho)
{
  FILE *f = fopen(caminho, "rb");
  long n;
  char *b;
  if (!f) return NULL;
  fseek(f, 0, SEEK_END);
  n = ftell(f);
  fseek(f, 0, SEEK_SET);
  if (n < 0) { fclose(f); return NULL; }
  b = (char *)malloc((size_t)n + 1);
  if (!b) { fclose(f); return NULL; }
  if (fread(b, 1, (size_t)n, f) != (size_t)n) { free(b); fclose(f); return NULL; }
  b[n] = '\0';
  fclose(f);
  return b;
}

static int escreve_csv(const dae_spec *S, const dae_graph *G, const dae_metrics *M,
                       const dae_series *R, int estado, const char *saida)
{
  const int32_t precisa = dae_csv(S, G, M, R, estado, NULL, 0) + 1;
  char *buf = (char *)malloc((size_t)precisa);
  FILE *f;
  if (!buf) return 0;
  dae_csv(S, G, M, R, estado, buf, precisa);
  if (saida) {
    f = fopen(saida, "wb");
    if (!f) { free(buf); return 0; }
    fwrite(buf, 1, strlen(buf), f);
    fclose(f);
  } else {
    fwrite(buf, 1, strlen(buf), stdout);
  }
  free(buf);
  return 1;
}

/* Setor magnetico: outra fisica, outro CSV (dae_numag_run, dae_numag_csv). */
static int rodar_magnetico(const dae_spec *S, const char *saida)
{
  dae_numag_serie R;
  dae_status st = dae_numag_run(S, &R, NULL, NULL);
  int32_t precisa;
  char *buf;
  FILE *f;
  if (st != DAE_OK) {
    fprintf(stderr, "setor magnetico falhou: %s\n", dae_strerror(st));
    if (st == DAE_ERR_TOOBIG)
      fprintf(stderr,
        "  o guarda de passo nao convergiu numa fatia de %g km.\n"
        "  aumente time.nt (fatias mais curtas) ou, se souber justificar,\n"
        "  afrouxe magnetic.tolerance.\n", S->t1 / (double)S->nt);
    return 1;
  }
  precisa = dae_numag_csv(S, &R, NULL, 0) + 1;
  buf = (char *)malloc((size_t)precisa);
  if (!buf) { dae_numag_serie_free(&R); return 1; }
  dae_numag_csv(S, &R, buf, precisa);
  if (saida) {
    f = fopen(saida, "wb");
    if (!f) { free(buf); dae_numag_serie_free(&R); return 1; }
    fwrite(buf, 1, strlen(buf), f);
    fclose(f);
  } else {
    fwrite(buf, 1, strlen(buf), stdout);
  }
  free(buf);
  dae_numag_serie_free(&R);
  return 0;
}

/* `ghost traj` — roda o ensemble de trajetorias e escreve rho.
 *
 * Existe para a validacao das trajetorias: a media delas tem de reproduzir a
 * solucao exata de Lindblad, e quem calcula a solucao exata e o pacote Wolfram,
 * por um metodo que nao e trajetoria nenhuma. Sem este arquivo os dois lados
 * nao se encontram, e "as trajetorias convergiram" fica indistinguivel de "as
 * trajetorias convergiram para a coisa errada" — com barra de erro respeitavel
 * nos dois casos. */
static int rodar_trajetorias(const dae_spec *S, const char *saida)
{
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  dae_traj_cfg cfg;
  dae_rho_acc A;
  dae_status st;
  FILE *f = NULL;
  int32_t na, s2, i, j, precisa;
  int nu_merg = 0;
  double nu_re[2 * DAE_NU_MAX_SABORES], nu_im[2 * DAE_NU_MAX_SABORES];
  char *cab;
  double escala = 1.0;

  if (S->n_traj <= 0) {
    fprintf(stderr, "este spec nao pede trajetorias (trajectories.n_traj = 0)\n");
    return 1;
  }
  if (S->saida_traj != DAE_SAIDA_ACUMULAR_RHO) {
    fprintf(stderr, "`traj` escreve rho, entao exige output_mode accumulate_rho\n");
    return 1;
  }
  /* MESMO CRITERIO DE dae_run.c: uma lista de arestas explicita vence o
     gerador. Sem isto `traj` rejeita todo spec de edgelist — dae_graph_build
     devolve DAE_ERR_PARAM para DAE_G_EDGELIST (dae_graph.c, "use
     dae_graph_from_edges") — enquanto `run`, que passa por dae_run, aceita o
     mesmo arquivo. Os dois caminhos precisam concordar sobre o que e um grafo
     valido, senao o spec deixa de ser a moeda unica. */
  if (S->nu_ativo) {
    /* Mesmo caminho de dae_run: com o bloco `neutrino`, quem monta a
       adjacencia e a fisica. O estado inicial vira VETOR (no mergulho ele e
       (psi, -i psi)/sqrt(2), que nao e um sitio) e o colapso passa a ser no
       PAR — ver dae_traj.h, campo nf_par. */
    double hre[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
    double him[DAE_NU_MAX_SABORES * DAE_NU_MAX_SABORES];
    /* memset ANTES de preencher: dae_graph_build e dae_graph_from_edges zeram
       o grafo na entrada, e este caminho nao passa por nenhum dos dois. Sem
       isto os campos que a fisica de neutrinos nao usa (coordenadas, rotulos)
       ficam com lixo da pilha, e dae_graph_free tenta liberar ponteiro
       invalido — segfault DEPOIS do resultado sair certo. */
    memset(&G, 0, sizeof(G));
    st = dae_nu_hamiltoniano(&S->nu, S->nu_perfil.rho, S->nu_perfil.ye, hre, him);
    if (st != DAE_OK) { fprintf(stderr, "H de neutrino: %s\n", dae_strerror(st)); return 1; }
    st = dae_nu_adjacencia(hre, him, S->nu.n_flavors, &G.A, &nu_merg);
    if (st == DAE_OK) {
      G.n = G.A.n;
      G.nmod = 1;
      G.module_of = (int32_t *)calloc((size_t)G.n, sizeof(int32_t));
      if (!G.module_of) st = DAE_ERR_ALLOC;
    }
  } else if (S->n_edges > 0) {
    if (S->gen.n <= 0) {
      fprintf(stderr, "grafo: `edges` exige `n` em graph\n");
      return 1;
    }
    st = dae_graph_from_edges(&G, S->gen.n, S->ei, S->ej, S->w, S->n_edges);
  } else {
    dae_gen_params gp = S->gen; gp.seed = S->seed;
    st = dae_graph_build(&G, &gp);
  }
  if (st != DAE_OK) { fprintf(stderr, "grafo: %s\n", dae_strerror(st)); return 1; }
  if (S->nu_ativo)
    st = dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, 1.0, DAE_NORM_NONE, 0, &escala);
  else
    st = dae_hamiltonian(&H, &G.A, S->ham, S->gamma, S->norm, S->lanczos_steps, &escala);
  if (st != DAE_OK) { fprintf(stderr, "H: %s\n", dae_strerror(st)); return 1; }
  st = dae_cheb_init(&W, &H, S->lanczos_steps);
  if (st != DAE_OK) { fprintf(stderr, "cheb: %s\n", dae_strerror(st)); return 1; }

  dae_traj_cfg_default(&cfg);
  if (S->nu_ativo) {
    dae_nu_inicial(S->nu_sabor, S->nu.n_flavors, nu_merg, nu_re, nu_im);
    cfg.init_re = nu_re;
    cfg.init_im = nu_im;
    cfg.nf_par = nu_merg ? S->nu.n_flavors : 0;
  }
  cfg.gamma_deph = S->gamma_deph;
  cfg.n_traj = S->n_traj;
  cfg.rho_stride = S->rho_stride;
  cfg.nt = S->nt;
  /* Com o bloco neutrino, t1 esta em km e o propagador anda em 1/eV — a mesma
     conversao que dae_run faz. Sem isto o dt sairia mil vezes menor que o
     pedido e a defasagem pareceria fraca sem motivo. */
  cfg.dt = (S->nu_ativo ? S->t1 * DAE_NU_T_POR_KM : S->t1) / (double)S->nt;
  cfg.sitio_inicial = S->init_site >= 0 ? S->init_site : 0;
  cfg.saida = DAE_SAIDA_ACUMULAR_RHO;
  na = dae_traj_amostras(&cfg);

  st = dae_rho_acc_init(&A, G.n, na);
  if (st != DAE_OK) { fprintf(stderr, "acumulador: %s\n", dae_strerror(st)); return 1; }
  st = dae_traj_ensemble(&cfg, S->seed, &W, &H, &A, NULL, NULL);
  if (st != DAE_OK) { fprintf(stderr, "ensemble: %s\n", dae_strerror(st)); return 1; }
  dae_rho_acc_finalizar(&A);

  f = saida ? fopen(saida, "wb") : stdout;
  if (!f) { fprintf(stderr, "nao consegui abrir %s\n", saida); return 1; }
  precisa = dae_csv_cabecalho_psi(S, S->seed, G.n, na, S->n_traj, NULL, 0) + 1;
  cab = (char *)malloc((size_t)precisa);
  dae_csv_cabecalho_psi(S, S->seed, G.n, na, S->n_traj, cab, precisa);
  fwrite(cab, 1, strlen(cab), f);
  free(cab);
  fprintf(f, "#! conteudo rho_medio\n");
  fprintf(f, "#! dt %.17g\n", cfg.dt);
  /* A escala do nucleo vai JUNTO: ela e estimativa de Lanczos, e o Wolfram
     calcularia o raio espectral exato. Sem impo-la, os dois compararao
     dinamicas de hamiltonianos diferentes e a discordancia apareceria como
     desdobramento errado. */
  fprintf(f, "#! scale %.17g\n", escala);
  fprintf(f, "#! gamma_deph %.17g\n", cfg.gamma_deph);
  fprintf(f, "s,t,i,j,re,im\n");
  for (s2 = 0; s2 < na; ++s2) {
    const double tt = (double)(s2 * cfg.rho_stride) * cfg.dt;
    const size_t off = (size_t)s2 * (size_t)G.n * (size_t)G.n;
    for (i = 0; i < G.n; ++i)
      for (j = 0; j < G.n; ++j)
        fprintf(f, "%d,%.17g,%d,%d,%.17g,%.17g\n", s2, tt, i, j,
                A.re[off + (size_t)i * (size_t)G.n + (size_t)j],
                A.im[off + (size_t)i * (size_t)G.n + (size_t)j]);
  }
  if (saida) fclose(f);
  fprintf(stderr, "  traj: n=%d  amostras=%d  trajetorias=%d  gamma=%g -> %s\n",
          G.n, na, S->n_traj, S->gamma_deph, saida ? saida : "(stdout)");

  dae_rho_acc_free(&A);
  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
  return 0;
}


/* `ghost run` com bloco `scan`: um espectro, nao um ponto.
 *
 * Experimento mede espectro, e a razao nao e conveniencia: com delta_cp = 0 e
 * materia a assimetria nu/antinu ja vale +0,42, entao um ponto so NAO separa
 * fase CP de efeito de materia. Quem separa e a dependencia em E — as duas
 * causas entram com formas diferentes. Por isso `both_signs` roda os dois
 * sinais na mesma varredura e emite A_CP junto: comparar duas corridas feitas
 * em momentos diferentes convidaria a comparar parametros diferentes.
 */
static int rodar_varredura(dae_spec *S, const char *saida)
{
  dae_graph G;
  dae_metrics M;
  dae_series R;
  FILE *f;
  char *cab;
  int32_t p, a, precisa, nf = S->nu.n_flavors;
  double e_ini = S->scan_min, e_fim = S->scan_max;
  dae_status st;
  static const char *NOME[3] = { "e", "mu", "tau" };

  f = saida ? fopen(saida, "wb") : stdout;
  if (!f) { fprintf(stderr, "nao consegui abrir %s\n", saida); return 1; }

  precisa = dae_csv(S, NULL, NULL, NULL, 0, NULL, 0);
  (void)precisa;
  fprintf(f, "#! ghost %s\n#! core_hash %s\n#! implementacao c\n",
          DAE_VERSION, DAE_CORE_HASH);
  fprintf(f, "#! conteudo espectro\n");
  { int32_t need = dae_spec_canonical(S, NULL, 0);
    cab = (char *)malloc((size_t)need + 1);
    if (cab) { dae_spec_canonical(S, cab, need + 1);
               fprintf(f, "#! spec %s\n", cab); free(cab); } }
  fprintf(f, "E_MeV");
  for (a = 0; a < nf; ++a) fprintf(f, ",p_%s", NOME[a]);
  for (a = 0; a < nf; ++a) fprintf(f, ",m_%d", a + 1);
  for (a = 0; a < nf; ++a) fprintf(f, ",inc_%s", NOME[a]);
  if (S->scan_ambos) {
    for (a = 0; a < nf; ++a) fprintf(f, ",bar_p_%s", NOME[a]);
    for (a = 0; a < nf; ++a) fprintf(f, ",acp_%s", NOME[a]);
  }
  fprintf(f, "\n");

  for (p = 0; p < S->scan_pontos; ++p) {
    const double u = (double)p / (double)(S->scan_pontos - 1);
    const double E = S->scan_log ? e_ini * pow(e_fim / e_ini, u)
                                 : e_ini + u * (e_fim - e_ini);
    double nu[3 * 3];
    int sinal;

    for (sinal = 0; sinal <= (S->scan_ambos ? 1 : 0); ++sinal) {
      S->nu.energy_eV = E * 1.0e6;
      S->nu.antineutrino = sinal;
      st = dae_run(S, S->seed, &G, &M, &R, NULL, NULL);
      if (st != DAE_OK) {
        fprintf(stderr, "E = %g MeV falhou: %s\n", E, dae_strerror(st));
        if (saida) fclose(f);
        return 1;
      }
      { const size_t o = (size_t)(R.nt - 1) * (size_t)nf;
        if (sinal == 0) {
          fprintf(f, "%.17g", E);
          for (a = 0; a < nf; ++a) { nu[a] = R.nu_sabor[o + (size_t)a];
                                     fprintf(f, ",%.17g", nu[a]); }
          for (a = 0; a < nf; ++a) fprintf(f, ",%.17g", R.nu_massa[o + (size_t)a]);
          for (a = 0; a < nf; ++a) fprintf(f, ",%.17g", R.nu_inc[o + (size_t)a]);
        } else {
          for (a = 0; a < nf; ++a) fprintf(f, ",%.17g", R.nu_sabor[o + (size_t)a]);
          for (a = 0; a < nf; ++a) {
            const double b = R.nu_sabor[o + (size_t)a], soma = nu[a] + b;
            /* A_CP = (P - Pbar)/(P + Pbar). Soma nula sai como nan, nunca 0:
               0 e valor fisicamente valido e significaria "sem assimetria".
               Mesma regra de p_alvo, CONVENTIONS.md 7.2. */
            if (soma > 0.0) fprintf(f, ",%.17g", (nu[a] - b) / soma);
            else            fprintf(f, ",nan");
          }
        }
      }
      dae_series_free(&R); dae_graph_free(&G);
    }
    fprintf(f, "\n");
  }
  if (saida) fclose(f);
  fprintf(stderr, "  espectro: %d energias de %g a %g MeV%s -> %s\n",
          S->scan_pontos, e_ini, e_fim,
          S->scan_ambos ? " (nu e antinu)" : "", saida ? saida : "(stdout)");
  return 0;
}

int main(int argc, char **argv)
{
  dae_spec S;
  dae_error err;
  dae_graph G;
  dae_metrics M;
  dae_series R;
  char *texto;
  const char *saida = NULL;
  int estado = 0, i, realizacoes = -1;
  dae_status st;

  if (argc < 3 || (strcmp(argv[1], "run") != 0 && strcmp(argv[1], "traj") != 0)) {
    fprintf(stderr,
      "ghost %s (nucleo %s)\n"
      "uso: ghost run  <spec.json> [--saida arq.csv] [--estado] [--realizacoes N]\n"
      "     ghost traj <spec.json> [--saida rho.csv]\n",
      DAE_VERSION, DAE_CORE_HASH);
    return 2;
  }
  for (i = 3; i < argc; ++i) {
    if (strcmp(argv[i], "--saida") == 0 && i + 1 < argc) saida = argv[++i];
    else if (strcmp(argv[i], "--estado") == 0) estado = 1;
    else if (strcmp(argv[i], "--realizacoes") == 0 && i + 1 < argc) realizacoes = atoi(argv[++i]);
    else { fprintf(stderr, "argumento desconhecido: %s\n", argv[i]); return 2; }
  }

  texto = le_arquivo(argv[2]);
  if (!texto) { fprintf(stderr, "nao consegui ler %s\n", argv[2]); return 1; }

  st = dae_spec_parse(&S, texto, &err);
  free(texto);
  if (st != DAE_OK) {
    fprintf(stderr, "%s:%d:%d: %s\n", argv[2], (int)err.line, (int)err.col, err.msg);
    return 1;
  }
  if (S.core_hash[0] && !dae_spec_hash_confere(&S)) {
    fprintf(stderr, "aviso: spec gerado pelo nucleo %s, executando com %s — "
                    "o resultado nao e rastreavel ate aquele codigo\n",
            S.core_hash, DAE_CORE_HASH);
  }
  if (realizacoes > 0) S.realizations = realizacoes;

  if (S.numag_ativo) {
    int r;
    if (strcmp(argv[1], "run") != 0) {
      fprintf(stderr, "o setor magnetico roda com `ghost run`, nao com `%s`\n", argv[1]);
      dae_spec_free(&S);
      return 2;
    }
    r = rodar_magnetico(&S, saida);
    dae_spec_free(&S);
    return r;
  }

  if (S.scan_pontos > 0 && strcmp(argv[1], "run") == 0) {
    const int r = rodar_varredura(&S, saida);
    dae_spec_free(&S);
    return r;
  }

  if (strcmp(argv[1], "traj") == 0) {
    const int r = rodar_trajetorias(&S, saida);
    dae_spec_free(&S);
    return r;
  }

  for (i = 0; i < S.realizations; ++i) {
    const uint64_t semente = S.seed + (uint64_t)i;
    char nome[512];
    st = dae_run(&S, semente, &G, &M, &R, NULL, NULL);
    if (st != DAE_OK) {
      fprintf(stderr, "realizacao %d falhou: %s\n", i, dae_strerror(st));
      if (st == DAE_ERR_TOOBIG && S.nu_perfil_ativo) {
        /* O guarda de adiabaticidade recusou. Recusar e o comportamento CERTO —
           a alternativa e um numero errado em silencio —, mas a mensagem tem de
           dizer o que fazer, senao parece defeito do programa. */
        fprintf(stderr,
          "  o guarda de adiabaticidade nao convergiu numa fatia.\n"
          "  cada fatia mede t1/nt = %g km, e perto da fonte a densidade exige\n"
          "  passos bem menores. tente, nesta ordem:\n"
          "    - aumentar time.nt: fatias mais curtas;\n"
          "    - aumentar matter.r_start_km: comecar onde a densidade ja caiu.\n"
          "      com V dominando Dm2/2E o nu_e JA E autoestado de materia, e\n"
          "      aquele trecho nao converte sabor nenhum — so custa fase;\n"
          "    - afrouxar matter.tolerance, se souber justificar.\n",
          S.t1 / (double)S.nt);
      }
      dae_spec_free(&S);
      return 1;
    }
    if (S.realizations > 1 && saida) {
      snprintf(nome, sizeof(nome), "%s.r%04d.csv", saida, i);
      escreve_csv(&S, &G, &M, &R, estado, nome);
    } else {
      escreve_csv(&S, &G, &M, &R, estado, saida);
    }
    dae_series_free(&R);
    dae_graph_free(&G);
  }
  dae_spec_free(&S);
  return 0;
}
