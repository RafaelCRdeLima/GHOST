#include "dae_csv.h"

#include <stdio.h>      /* snprintf; ver a emenda em dae_spec.c */
#include <stdlib.h>
#include <string.h>

typedef struct { char *b; int32_t cap, n; } dae_cbuf;

/* `m` e o tamanho PRETENDIDO por snprintf, que pode passar do buffer quando
   houve truncamento — copiar `m` bytes de um buffer de 256 e leitura fora de
   limite. O clamp para `lim` fecha isso; strings longas nao passam por aqui,
   passam por cput. */
static void ce(dae_cbuf *B, const char *t, int m, int32_t lim)
{
  int32_t i;
  if (m < 0) return;
  if (m > lim) m = lim;
  for (i = 0; i < m; ++i) {
    if (B->n + 1 < B->cap) B->b[B->n] = t[i];
    ++B->n;
  }
}

/* Anexa uma string de qualquer tamanho, sem buffer intermediario. O spec
   canonico tem centenas de bytes e nao cabe num temporario de pilha, que o
   truncaria. */
static void cput(dae_cbuf *B, const char *s)
{
  for (; *s; ++s) {
    if (B->n + 1 < B->cap) B->b[B->n] = *s;
    ++B->n;
  }
}
static void cs(dae_cbuf *B, const char *f, const char *v)
{ char t[256]; ce(B, t, snprintf(t, sizeof(t), f, v), (int32_t)sizeof(t) - 1); }
static void ci(dae_cbuf *B, const char *f, int32_t v)
{ char t[64]; ce(B, t, snprintf(t, sizeof(t), f, (int)v), (int32_t)sizeof(t) - 1); }
static void cd(dae_cbuf *B, const char *f, double v)
{ char t[64]; ce(B, t, snprintf(t, sizeof(t), f, v), (int32_t)sizeof(t) - 1); }
static void cu(dae_cbuf *B, const char *f, uint64_t v)
{ char t[64]; ce(B, t, snprintf(t, sizeof(t), f, (unsigned long long)v), (int32_t)sizeof(t) - 1); }

int32_t dae_csv_cabecalho_psi(const dae_spec *S, uint64_t semente_base,
                              int32_t n, int32_t n_amostras, int32_t n_traj,
                              char *buf, int32_t cap)
{
  dae_cbuf B;
  int32_t need;
  char *tmp;
  if (!S) return 0;
  B.b = buf; B.cap = cap; B.n = 0;
  if (cap > 0) buf[0] = '\0';

  cs(&B, "#! ghost %s\n", DAE_VERSION);
  cs(&B, "#! core_hash %s\n", DAE_CORE_HASH);
  cs(&B, "#! implementacao %s\n", "c");
  cs(&B, "#! conteudo %s\n", "psi_trajetorias");
  /* A semente-base fica EXPLICITA aqui e nao so dentro do spec: e ela, com o
     indice da trajetoria, que gera cada fluxo (dae_traj_semente). Sem ela o
     acervo nao e regeneravel, e um cache que nao se regenera e um acervo. */
  {
    char sb[32];
    int32_t k = 0;
    uint64_t v = semente_base;
    char rev[32];
    if (v == 0) { rev[k++] = '0'; }
    while (v > 0 && k < 20) { rev[k++] = (char)('0' + (int)(v % 10ULL)); v /= 10ULL; }
    { int32_t q; for (q = 0; q < k; ++q) sb[q] = rev[k - 1 - q]; sb[k] = '\0'; }
    cs(&B, "#! semente_base %s\n", sb);
  }
  ci(&B, "#! n %d\n", n);
  ci(&B, "#! n_amostras %d\n", n_amostras);
  ci(&B, "#! n_traj %d\n", n_traj);
  cput(&B, "#! layout f64 re[n_amostras*n] depois im[n_amostras*n], por trajetoria\n");

  need = dae_spec_canonical(S, NULL, 0);
  tmp = (char *)malloc((size_t)need + 1);
  if (tmp) {
    dae_spec_canonical(S, tmp, need + 1);
    cput(&B, "#! spec ");
    cput(&B, tmp);
    cput(&B, "\n");
    free(tmp);
  }
  /* FECHA A STRING. Nem `ce` nem `cput` escrevem o terminador — ele e do
     emissor. Sem ele, `strlen` le fora do buffer e o lixo entra no arquivo
     como bytes binarios no meio do cabecalho. `t97` exige strlen == devolvido. */
  if (cap > 0) buf[B.n < cap ? B.n : cap - 1] = '\0';
  return B.n;
}

int32_t dae_csv(const dae_spec *S, const dae_graph *G, const dae_metrics *M,
                const dae_series *R, int incluir_estado, char *buf, int32_t cap)
{
  dae_cbuf B;
  int32_t s, i, m;
  if (!S || !G || !R) return 0;
  B.b = buf; B.cap = cap; B.n = 0;
  if (cap > 0) buf[0] = '\0';

  cs(&B, "#! ghost %s\n", DAE_VERSION);
  cs(&B, "#! core_hash %s\n", DAE_CORE_HASH);
  /* QUEM CALCULOU, que e outra pergunta que `core_hash` nao responde. O
     core_hash diz de qual spec e de qual versao do projeto o resultado veio; a
     implementacao diz qual programa produziu os numeros. O pacote Wolfram
     escreve "wolfram" aqui. Sem esta linha, um CSV do Wolfram reimportado se
     apresentaria como saida do nucleo em C — procedencia falsa, que e
     exatamente o que o cabecalho existe para impedir. */
  cs(&B, "#! implementacao %s\n", "c");
  cu(&B, "#! graph_fingerprint %llu\n", R->info.fingerprint);
  ci(&B, "#! n %d\n", G->n);
  ci(&B, "#! nnz %d\n", G->A.nnz);
  ci(&B, "#! nmod %d\n", G->nmod);
  cd(&B, "#! scale %.17g\n", R->info.escala);
  cd(&B, "#! spectral_lo %.17g\n", R->info.lo);
  cd(&B, "#! spectral_hi %.17g\n", R->info.hi);
  cd(&B, "#! cheb_a %.17g\n", R->info.a);
  cd(&B, "#! cheb_b %.17g\n", R->info.b);
  cd(&B, "#! dt %.17g\n", R->info.dt);
  /* alpha explícito: é o argumento que a fórmula da ordem usa, e a grade
     exportada o varre numa faixa que os testes do núcleo nunca viram. */
  cd(&B, "#! alpha %.17g\n", R->info.alpha);
  ci(&B, "#! cheb_k %d\n", R->info.k_used);
  ci(&B, "#! lanczos_used %d\n", R->info.lanczos_used);
  if (M) {
    cd(&B, "#! lambda2 %.17g\n", M->lambda2);
    cd(&B, "#! lambda2_residual %.17g\n", M->lambda2_residual);
    ci(&B, "#! lambda2_converged %d\n", M->lambda2_converged);
    cd(&B, "#! modularity_Q %.17g\n", M->modularity_Q);
    cd(&B, "#! mean_degree %.17g\n", M->mean_degree);
    cd(&B, "#! mean_path_len %.17g\n", M->mean_path_len);
    ci(&B, "#! n_edges %d\n", M->n_edges);
    ci(&B, "#! n_components %d\n", M->n_components);
  }
  /* Religacoes que nao acharam destino livre. Sai SEMPRE, inclusive zero: uma
     tentativa que colide consome sorteios extras, entao este contador e a
     unica testemunha externa de quantos sorteios o gerador gastou. Duas
     implementacoes podem chegar a mesma digital por caminhos diferentes de
     sorteio so se este numero tambem bater. */
  ci(&B, "#! rewire_failed %d\n", G->n_rewire_failed);
  { /* o spec canônico entra em linha única: um resultado solto continua
       rastreável até o JSON que o gerou */
    int32_t precisa = dae_spec_canonical(S, NULL, 0) + 1;
    char *tmp = (char *)malloc((size_t)precisa);
    if (tmp) {
      dae_spec_canonical(S, tmp, precisa);
      cput(&B, "#! spec ");
      cput(&B, tmp);
      cput(&B, "\n");
      free(tmp);
    }
  }

  cput(&B, "t,norm,ipr,coh_l1,p_target");
  if (R->nu_sabor) {
    /* Colunas NOMEADAS. `pmod0` obrigaria quem le a saber, de fora do arquivo,
       que o modulo 0 e o nu_e — e o arquivo existe justamente para nao
       precisar disso. As tres familias sao diferentes e o nome diz qual:
         p_*    sabor coerente, no ponto;
         m_*    base de massa, a metrica de convergencia de um perfil;
         inc_*  INCOERENTE, o que um detector astrofisico de fato recebe. */
    static const char *NOME[3] = { "e", "mu", "tau" };
    for (m = 0; m < R->nu_nf; ++m) cs(&B, ",p_%s", NOME[m]);
    for (m = 0; m < R->nu_nf; ++m) ci(&B, ",m_%d", m + 1);
    for (m = 0; m < R->nu_nf; ++m) cs(&B, ",inc_%s", NOME[m]);
  } else {
    for (m = 0; m < G->nmod; ++m) ci(&B, ",pmod%d", m);
  }
  cput(&B, "\n");
  for (s = 0; s < R->nt; ++s) {
    const double *r = R->scal + (size_t)s * DAE_S_NCOL;
    cd(&B, "%.17g", R->t[s]);
    cd(&B, ",%.17g", r[DAE_S_NORMA]);
    cd(&B, ",%.17g", r[DAE_S_IPR]);
    cd(&B, ",%.17g", r[DAE_S_COH_L1]);
    /* p_alvo sem alvo é NaN, e sai como "nan" — nunca como 0, que é valor
       fisicamente válido. Ver CONVENTIONS.md, parte 7.2. */
    if (r[DAE_S_P_ALVO] == r[DAE_S_P_ALVO]) cd(&B, ",%.17g", r[DAE_S_P_ALVO]);
    else cput(&B, ",nan");
    if (R->nu_sabor) {
      const size_t off = (size_t)s * (size_t)R->nu_nf;
      for (m = 0; m < R->nu_nf; ++m) cd(&B, ",%.17g", R->nu_sabor[off + (size_t)m]);
      for (m = 0; m < R->nu_nf; ++m) cd(&B, ",%.17g", R->nu_massa[off + (size_t)m]);
      for (m = 0; m < R->nu_nf; ++m) cd(&B, ",%.17g", R->nu_inc[off + (size_t)m]);
    } else {
      for (m = 0; m < G->nmod; ++m)
        cd(&B, ",%.17g", R->pmod[(size_t)s * (size_t)G->nmod + (size_t)m]);
    }
    cput(&B, "\n");
  }

  if (incluir_estado && R->psi_re) {
    cput(&B, "\n# estado final\nj,re,im\n");
    for (i = 0; i < G->n; ++i) {
      ci(&B, "%d", i);
      cd(&B, ",%.17g", R->psi_re[i]);
      cd(&B, ",%.17g\n", R->psi_im[i]);
    }
  }
  if (cap > 0) buf[B.n < cap ? B.n : cap - 1] = '\0';
  return B.n;
}

int32_t dae_numag_csv(const dae_spec *S, const dae_numag_serie *R,
                      char *buf, int32_t cap)
{
  dae_cbuf B;
  int32_t s;
  if (!S || !R) return 0;
  B.b = buf; B.cap = cap; B.n = 0;
  if (cap > 0) buf[0] = '\0';

  cs(&B, "#! ghost %s\n", DAE_VERSION);
  cs(&B, "#! core_hash %s\n", DAE_CORE_HASH);
  cs(&B, "#! implementacao %s\n", "c");
  cs(&B, "#! setor %s\n", "magnetico");
  ci(&B, "#! subdiv_max %d\n", R->subdiv_max);
  cs(&B, "#! integrador %s\n", S->numag_integrador ? "janelas" : "forca_bruta");
  if (R->n_traj > 0) {
    ci(&B, "#! trajetorias %d\n", R->n_traj);
    cd(&B, "#! taxa_colisao_eV %.17g\n", S->numag_gamma_eV);
    cd(&B, "#! salto_por_km %.17g\n", R->salto_por_km);
    cd(&B, "#! saltos_medio %.17g\n", R->saltos_medio);
  }
  if (S->numag_integrador) {
    const dae_numag_janelas_info *J = &R->janelas;
    cs(&B, "#! coerencia %s\n", S->numag_coerencia ? "incoerente" : "coerente");
    ci(&B, "#! janelas %d\n", J->n_janelas);
    ci(&B, "#! trechos_adiabaticos %d\n", J->n_trechos);
    cd(&B, "#! comprimento_janela_km %.17g\n", J->comprimento_janela_km);
    /* Apesar do nome, "cota_adiabatica" é um ESTIMADOR de 1a ordem com a
       variação total AMOSTRADA, não uma cota certificada (CONVENTIONS 13.5b) */
    cd(&B, "#! cota_adiabatica %.17g\n", J->cota_adiabatica);
    cd(&B, "#! eps_fora_max %.17g\n", J->eps_fora_max);
    cd(&B, "#! saltos_tabela %.17g\n", J->saltos_tabela);
    ci(&B, "#! cruzamentos %d\n", J->n_cruzamentos);
    cd(&B, "#! cruzamento_s_km %.17g\n", J->cruzamento_s_km);
    ci(&B, "#! contatos %d\n", J->n_contatos);
    ci(&B, "#! intervalos_nulos %d\n", J->n_intervalos_nulos);
    cd(&B, "#! comprimento_nulo_km %.17g\n", J->comprimento_nulo_km);
  }
  {
    int32_t precisa = dae_spec_canonical(S, NULL, 0) + 1;
    char *tmp = (char *)malloc((size_t)precisa);
    if (tmp) {
      dae_spec_canonical(S, tmp, precisa);
      cput(&B, "#! spec ");
      cput(&B, tmp);
      cput(&B, "\n");
      free(tmp);
    }
  }
  if (R->n_traj > 0) {
    cput(&B, "s_km,norm,p_L,p_R,p_R_err,rho_LR_re,rho_LR_im\n");
    for (s = 0; s < R->nt; ++s) {
      cd(&B, "%.17g", R->s_km[s]);
      cd(&B, ",%.17g", R->norma[s]);
      cd(&B, ",%.17g", R->p_l[s]);
      cd(&B, ",%.17g", R->p_r[s]);
      cd(&B, ",%.17g", R->p_r_erro[s]);
      cd(&B, ",%.17g", R->coh_re[s]);
      cd(&B, ",%.17g\n", R->coh_im[s]);
    }
  } else {
    cput(&B, "s_km,norm,p_L,p_R,subdiv\n");
    for (s = 0; s < R->nt; ++s) {
      cd(&B, "%.17g", R->s_km[s]);
      cd(&B, ",%.17g", R->norma[s]);
      cd(&B, ",%.17g", R->p_l[s]);
      cd(&B, ",%.17g", R->p_r[s]);
      ci(&B, ",%d\n", R->subdiv[s]);
    }
  }
  if (R->tem_estado) {
    /* O estado final COMPLEXO, para o operador de transferência: as
       populações não veem fase relativa. */
    int32_t a;
    cput(&B, "\n# estado final\nj,re,im\n");
    for (a = 0; a < 2; ++a) {
      ci(&B, "%d", a);
      cd(&B, ",%.17g", R->psi_re[a]);
      cd(&B, ",%.17g\n", R->psi_im[a]);
    }
  }
  if (cap > 0) buf[B.n < cap ? B.n : cap - 1] = '\0';
  return B.n;
}
