/* t97_modos_saida.c — os dois modos de saída das trajetórias dão o MESMO
 * resultado, bit a bit, e o teste que afirma isso pode falhar.
 *
 * A igualdade entre `accumulate_rho` e `archive_psi` só é possível porque a
 * DINÂMICA é a mesma e a ORDEM DA SOMA é contrato: `dae_rho_acc_somar` recusa
 * índice fora de sequência em vez de devolver outro número. Soma de ponto
 * flutuante não é associativa, e sem esse contrato a igualdade seria uma
 * coincidência que dura até alguém trocar o escalonador.
 *
 * ANTI-VACUIDADE, DUAS COMPANHEIRAS. A concordância passaria trivialmente se um
 * modo caísse silenciosamente no outro, e passaria trivialmente se os dois
 * produzissem observáveis nulos. Então:
 *   (a) os observáveis são comprovadamente não triviais, e com valor de
 *       referência: em gamma = 0 o ensemble tem de reproduzir EXATAMENTE o
 *       resultado unitário que o núcleo já calcula por outro caminho;
 *   (b) sabotagens deliberadas em UM dos caminhos, exigindo que a comparação
 *       morda. Se não morder, ela não está testando os dois caminhos.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---- coletor do modo ARQUIVAR_PSI ---------------------------------------
   Guarda em memória o que na produção iria para disco. O formato é o mesmo:
   n_amostras * n amplitudes por trajetória, re e im separados. */
typedef struct {
  double  *re, *im;
  int32_t  na, n, cap, guardadas;
} acervo;

static dae_status coletar(void *user, int32_t idx, const double *re,
                          const double *im, int32_t na, int32_t n)
{
  acervo *S = (acervo *)user;
  const size_t bloco = (size_t)na * (size_t)n;
  if (idx >= S->cap) return DAE_ERR_PARAM;
  S->na = na; S->n = n;
  memcpy(S->re + (size_t)idx * bloco, re, bloco * sizeof(double));
  memcpy(S->im + (size_t)idx * bloco, im, bloco * sizeof(double));
  S->guardadas += 1;
  return DAE_OK;
}

/* Coerência l1 de rho: sum_{i!=j} |rho_ij|. É o observável central das
   trajetórias, e é dele que a igualdade entre modos precisa ser bit a bit. */
static double coh_l1_de_rho(const dae_rho_acc *A, int32_t amostra)
{
  const size_t off = (size_t)amostra * (size_t)A->n * (size_t)A->n;
  double soma = 0.0;
  int32_t i, j;
  for (i = 0; i < A->n; ++i)
    for (j = 0; j < A->n; ++j) {
      if (i == j) continue;
      soma += sqrt(A->re[off + (size_t)i * (size_t)A->n + (size_t)j] *
                   A->re[off + (size_t)i * (size_t)A->n + (size_t)j] +
                   A->im[off + (size_t)i * (size_t)A->n + (size_t)j] *
                   A->im[off + (size_t)i * (size_t)A->n + (size_t)j]);
    }
  return soma;
}

static double traco(const dae_rho_acc *A, int32_t amostra)
{
  const size_t off = (size_t)amostra * (size_t)A->n * (size_t)A->n;
  double s = 0.0;
  int32_t i;
  for (i = 0; i < A->n; ++i) s += A->re[off + (size_t)i * (size_t)A->n + (size_t)i];
  return s;
}

int main(void)
{
  dae_test T;
  dae_graph G;
  dae_csr H;
  dae_cheb W;
  dae_gen_params gp;
  dae_traj_cfg cfg;
  dae_rho_acc A1, A2;
  acervo S;
  dae_status st;
  int32_t na, idx, k;
  size_t bloco, tot;
  int iguais;

  dae_test_begin(&T, "t97 modos de saida das trajetorias");

  /* Rede pequena e com costura: topologia de verdade, tamanho de teste. */
  dae_gen_params_default(&gp);
  gp.kind = DAE_G_MICROTUBULE;
  gp.n_par = 8; gp.n_perp = 5; gp.seam_shift = 2; gp.n_modules = 2;
  gp.ws_p = 0.0;
  gp.seed = 20260828ULL;
  st = dae_graph_build(&G, &gp);
  dae_test_ok(&T, st == DAE_OK, "gerou a rede (%s)", dae_strerror(st));

  st = dae_hamiltonian(&H, &G.A, DAE_H_ADJACENCY, 1.0, DAE_NORM_SPECTRAL, 20, NULL);
  dae_test_ok(&T, st == DAE_OK, "montou H (%s)", dae_strerror(st));
  st = dae_cheb_init(&W, &H, 20);
  dae_test_ok(&T, st == DAE_OK, "inicializou Chebyshev (%s)", dae_strerror(st));

  /* A config tem ponteiros (init_re, gamma_tab): preencher campo a
     campo deixaria os demais com lixo da pilha. */
  dae_traj_cfg_default(&cfg);
  cfg.gamma_deph = 0.35;
  cfg.n_traj = 24;
  cfg.rho_stride = 5;
  cfg.nt = 20;
  cfg.dt = 0.25;
  cfg.sitio_inicial = 0;
  na = dae_traj_amostras(&cfg);
  dae_test_ok(&T, na == 5, "5 amostras de tempo (obtive %d)", na);

  bloco = (size_t)na * (size_t)G.n;
  S.cap = cfg.n_traj; S.guardadas = 0; S.na = 0; S.n = 0;
  S.re = (double *)calloc(bloco * (size_t)cfg.n_traj, sizeof(double));
  S.im = (double *)calloc(bloco * (size_t)cfg.n_traj, sizeof(double));

  /* ---------- modo ACUMULAR_RHO ---------- */
  cfg.saida = DAE_SAIDA_ACUMULAR_RHO;
  st = dae_rho_acc_init(&A1, G.n, na);
  dae_test_ok(&T, st == DAE_OK, "acumulador 1 (%s)", dae_strerror(st));
  st = dae_traj_ensemble(&cfg, 777ULL, &W, &H, &A1, NULL, NULL);
  dae_test_ok(&T, st == DAE_OK, "ensemble acumulando rho (%s)", dae_strerror(st));

  /* ---------- modo ARQUIVAR_PSI ---------- */
  cfg.saida = DAE_SAIDA_ARQUIVAR_PSI;
  st = dae_traj_ensemble(&cfg, 777ULL, &W, &H, NULL, coletar, &S);
  dae_test_ok(&T, st == DAE_OK, "ensemble arquivando psi (%s)", dae_strerror(st));
  dae_test_ok(&T, S.guardadas == cfg.n_traj,
              "arquivou as %d trajetorias (obtive %d)", cfg.n_traj, S.guardadas);

  /* rho montado na ANALISE, em ordem de indice — o mesmo contrato. */
  st = dae_rho_acc_init(&A2, G.n, na);
  dae_test_ok(&T, st == DAE_OK, "acumulador 2 (%s)", dae_strerror(st));
  for (idx = 0; idx < cfg.n_traj; ++idx) {
    st = dae_rho_acc_somar(&A2, idx, S.re + (size_t)idx * bloco,
                                     S.im + (size_t)idx * bloco);
    if (st != DAE_OK) break;
  }
  dae_test_ok(&T, st == DAE_OK, "montou rho a partir do acervo (%s)", dae_strerror(st));

  dae_rho_acc_finalizar(&A1);
  dae_rho_acc_finalizar(&A2);

  /* ---------- A AFIRMACAO: bit a bit ---------- */
  tot = (size_t)na * (size_t)G.n * (size_t)G.n;
  iguais = (memcmp(A1.re, A2.re, tot * sizeof(double)) == 0) &&
           (memcmp(A1.im, A2.im, tot * sizeof(double)) == 0);
  dae_test_ok(&T, iguais, "os dois modos dao rho IDENTICO bit a bit");
  for (k = 0; k < na; ++k) {
    dae_test_ok(&T, coh_l1_de_rho(&A1, k) == coh_l1_de_rho(&A2, k),
                "coerencia l1 identica na amostra %d", k);
  }

  /* ---------- COMPANHEIRA (a): os observaveis nao sao triviais ---------- */
  {
    const double tr = traco(&A1, na - 1);
    const double coh = coh_l1_de_rho(&A1, na - 1);
    dae_test_near(&T, tr, 1.0, 1e-12, "traco de rho vale 1 (obtive %.17g)", tr);
    dae_test_ok(&T, coh > 0.5, "coerencia l1 nao e nula: %.6f", coh);
    dae_test_ok(&T, coh < (double)G.n - 1.0,
                "e nao e a maxima tampouco: %.6f < %d", coh, G.n - 1);
    dae_test_note("  coerencia l1 em t final: %.6f  (maximo possivel %d)", coh, G.n - 1);
  }

  /* Valor de REFERENCIA, e vem de outro caminho do nucleo: com gamma = 0 nao ha
     salto nenhum, toda trajetoria e a mesma, rho e puro, e a coerencia l1 tem de
     bater com a que dae_obs calcula sobre psi propagado unitariamente. Duas
     implementacoes, mesmo numero — e se o desdobramento estiver errado, ele
     erra AQUI, onde ha resposta conhecida. */
  {
    dae_rho_acc A0;
    dae_obs O;
    dae_obs_cfg oc;
    double *pr, *pi;
    int32_t p;
    double coh_traj, coh_unit;
    dae_cheb_info info;

    cfg.gamma_deph = 0.0;
    cfg.saida = DAE_SAIDA_ACUMULAR_RHO;
    dae_rho_acc_init(&A0, G.n, na);
    st = dae_traj_ensemble(&cfg, 777ULL, &W, &H, &A0, NULL, NULL);
    dae_test_ok(&T, st == DAE_OK, "ensemble com gamma = 0 (%s)", dae_strerror(st));
    dae_rho_acc_finalizar(&A0);
    coh_traj = coh_l1_de_rho(&A0, na - 1);

    pr = (double *)calloc((size_t)G.n, sizeof(double));
    pi = (double *)calloc((size_t)G.n, sizeof(double));
    pr[cfg.sitio_inicial] = 1.0;
    for (p = 0; p < cfg.nt; ++p) dae_cheb_step(&W, &H, cfg.dt, pr, pi, &info);
    memset(&oc, 0, sizeof(oc));
    oc.target = -1; oc.module_of = G.module_of; oc.nmod = G.nmod;
    dae_obs_alloc(&O, &oc, G.n);
    dae_obs_eval(pr, pi, G.n, &oc, &O);
    coh_unit = O.coh_l1;

    dae_test_near(&T, coh_traj, coh_unit, 1e-9,
                  "gamma = 0 reproduz o unitario: %.17g contra %.17g",
                  coh_traj, coh_unit);
    /* E a defasagem de fato mudou alguma coisa — senao o par acima estaria
       comparando o mesmo numero consigo mesmo por outro caminho. */
    dae_test_ok(&T, fabs(coh_l1_de_rho(&A1, na - 1) - coh_unit) > 1e-3,
                "e com gamma = 0.35 o resultado e OUTRO (%.6f contra %.6f)",
                coh_l1_de_rho(&A1, na - 1), coh_unit);

    dae_obs_free(&O);
    free(pr); free(pi);
    dae_rho_acc_free(&A0);
  }

  /* ---------- COMPANHEIRA (b): o teste morde? ---------- */
  /* (b.1) Somar fora de ordem NAO pode ser aceito em silencio: se fosse, a
     igualdade bit a bit acima seria coincidencia de um laco, e mudaria com o
     escalonador. */
  {
    dae_rho_acc Ax;
    dae_rho_acc_init(&Ax, G.n, na);
    st = dae_rho_acc_somar(&Ax, 1, S.re + bloco, S.im + bloco);
    dae_test_ok(&T, st == DAE_ERR_PARAM,
                "somar fora de ordem e RECUSADO (obtive %s)", dae_strerror(st));
    st = dae_rho_acc_somar(&Ax, 0, S.re, S.im);
    dae_test_ok(&T, st == DAE_OK, "e em ordem e aceito (%s)", dae_strerror(st));
    dae_rho_acc_free(&Ax);
  }

  /* (b.2) Sabotagens que MODELAM o que pode dar errado num arquivo de psi, e
     a comparacao tem de morder nas duas. Se nao morder, ela nao esta testando
     o caminho do acervo e a assercao principal e decorativa.

     Duas sabotagens que parecem naturais NAO mordem, e vale saber por que:

       1. Sabotar a amostra 0, onde psi = delta_{sitio inicial} e quase toda
          amplitude vale exatamente zero. nextafter(0) da um denormal de
          5e-324 que, multiplicado pelas outras amplitudes, desaparece antes
          de chegar a rho. "Alterei um valor" nao e "alterei um valor que
          importa".
       2. Sabotar 1 ULP da maior amplitude. Tambem nao morde, e o motivo e
          quantitativo: a media e sobre 24 trajetorias, entao uma perturbacao
          relativa de ~1e-16 numa contribuicao vira ~4e-18 no resultado, abaixo
          do ULP do proprio acumulado. Isto NAO e defeito da comparacao — e a
          resolucao dela, medida. Vale saber qual e.

     As duas sabotagens abaixo estao acima dessa resolucao, e correspondem a
     falhas reais: precisao perdida na escrita, e arquivo truncado. */
  {
    dae_rho_acc Ay;
    size_t alvo = bloco * 3;
    double guardado, melhor = -1.0;
    int32_t q;
    for (q = 0; q < G.n; ++q) {
      const size_t c = bloco * 3 + (size_t)(na - 1) * (size_t)G.n + (size_t)q;
      const double m = fabs(S.re[c]);
      if (m > melhor) { melhor = m; alvo = c; }
    }
    guardado = S.re[alvo];
    dae_test_ok(&T, melhor > 1e-3,
                "a sabotagem cai numa amplitude que importa: |re| = %.6g", melhor);

    /* (i) PRECISAO PERDIDA NA ESCRITA: o acervo gravado em f32 em vez de f64.
       E o erro mais provavel de um formato de arquivo, e o mais silencioso. */
    S.re[alvo] = (double)(float)guardado;
    dae_test_ok(&T, S.re[alvo] != guardado, "f32 de fato perde bits aqui");
    dae_rho_acc_init(&Ay, G.n, na);
    for (idx = 0; idx < cfg.n_traj; ++idx)
      dae_rho_acc_somar(&Ay, idx, S.re + (size_t)idx * bloco,
                                  S.im + (size_t)idx * bloco);
    dae_rho_acc_finalizar(&Ay);
    dae_test_ok(&T, memcmp(A1.re, Ay.re, tot * sizeof(double)) != 0,
                "precisao de f32 no acervo QUEBRA a igualdade — o teste morde");
    dae_rho_acc_free(&Ay);
    S.re[alvo] = guardado;

    /* (ii) ARQUIVO TRUNCADO: uma trajetoria a menos. O acumulador tem de
       terminar com outro numero E com outra contagem — as duas coisas, porque
       so a contagem seria conferivel sem olhar os numeros. */
    dae_rho_acc_init(&Ay, G.n, na);
    for (idx = 0; idx < cfg.n_traj - 1; ++idx)
      dae_rho_acc_somar(&Ay, idx, S.re + (size_t)idx * bloco,
                                  S.im + (size_t)idx * bloco);
    dae_test_ok(&T, Ay.somadas == cfg.n_traj - 1,
                "o acumulador conta quantas somou (%d)", Ay.somadas);
    dae_rho_acc_finalizar(&Ay);
    dae_test_ok(&T, memcmp(A1.re, Ay.re, tot * sizeof(double)) != 0,
                "acervo truncado QUEBRA a igualdade — o teste morde");
    dae_rho_acc_free(&Ay);
  }

  /* ---------- o acervo carrega procedencia ---------- */
  {
    dae_spec Sp;
    char *cab;
    int32_t need;
    dae_spec_default(&Sp);
    Sp.n_traj = cfg.n_traj;
    Sp.gamma_deph = cfg.gamma_deph;
    Sp.rho_stride = cfg.rho_stride;
    Sp.saida_traj = DAE_SAIDA_ARQUIVAR_PSI;
    need = dae_csv_cabecalho_psi(&Sp, 777ULL, G.n, na, cfg.n_traj, NULL, 0);
    cab = (char *)malloc((size_t)need + 1);
    dae_csv_cabecalho_psi(&Sp, 777ULL, G.n, na, cfg.n_traj, cab, need + 1);
    /* O tamanho devolvido e o tamanho da string. Sem isto, os strstr abaixo
       passam por sorte quando o malloc devolve memoria zerada, e a falta do
       terminador passaria despercebida. */
    dae_test_ok(&T, (int32_t)strlen(cab) == need,
                "o cabecalho esta fechado: strlen %d contra devolvido %d",
                (int)strlen(cab), (int)need);
    dae_test_ok(&T, strstr(cab, "#! core_hash ") != NULL, "acervo traz core_hash");
    dae_test_ok(&T, strstr(cab, "#! implementacao c") != NULL, "acervo traz implementacao");
    dae_test_ok(&T, strstr(cab, "#! semente_base 777") != NULL, "acervo traz a semente-base");
    dae_test_ok(&T, strstr(cab, "#! spec {") != NULL, "acervo traz o spec canonico");
    dae_test_ok(&T, strstr(cab, "\"output_mode\":\"archive_psi\"") != NULL,
                "e o spec canonico registra o modo escolhido");
    /* ANTI-VACUIDADE: o mesmo emissor com o outro modo tem de dizer o outro
       modo — senao a linha acima passaria com o campo fixo no texto. */
    Sp.saida_traj = DAE_SAIDA_ACUMULAR_RHO;
    dae_csv_cabecalho_psi(&Sp, 777ULL, G.n, na, cfg.n_traj, cab, need + 1);
    dae_test_ok(&T, strstr(cab, "\"output_mode\":\"accumulate_rho\"") != NULL,
                "e com o outro modo diz o outro modo");
    free(cab);
    dae_spec_free(&Sp);
  }

  free(S.re); free(S.im);
  dae_rho_acc_free(&A1); dae_rho_acc_free(&A2);
  dae_cheb_free(&W); dae_csr_free(&H); dae_graph_free(&G);
  return dae_test_end(&T);
}
