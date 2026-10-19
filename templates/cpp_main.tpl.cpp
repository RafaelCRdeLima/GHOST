/*__DAE_BANNER__*/

/*__DAE_CORE__*/

/* ------------------------------------------------------------------------
 * main gerado por templates/cpp_main.tpl.cpp
 *
 * Compile e rode sem instalar nada:
 *     g++ -O3 -fopenmp __ARQUIVO__ -o run && ./run
 *
 * O nucleo acima e o MESMO que roda no navegador — concatenado, nao
 * reescrito. Por isso comparar este programa com a pagina compara o mesmo
 * codigo sob dois compiladores, e nao duas implementacoes que concordam hoje.
 * ---------------------------------------------------------------------- */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

static const char *DAE_SPEC_JSON = R"DAEDALUSJSON(
/*__DAE_SPEC__*/
)DAEDALUSJSON";

static int escreve(const dae_spec *S, const dae_graph *G, const dae_metrics *M,
                   const dae_series *R, int estado, const char *saida)
{
  const int32_t precisa = dae_csv(S, G, M, R, estado, NULL, 0) + 1;
  char *buf = (char *)malloc((size_t)precisa);
  if (!buf) return 0;
  dae_csv(S, G, M, R, estado, buf, precisa);
  if (saida) {
    FILE *f = fopen(saida, "wb");
    if (!f) { free(buf); return 0; }
    fwrite(buf, 1, strlen(buf), f);
    fclose(f);
  } else {
    fwrite(buf, 1, strlen(buf), stdout);
  }
  free(buf);
  return 1;
}

int main(int argc, char **argv)
{
  dae_spec S0;
  dae_error err;
  const char *prefixo = 0;
  int estado = 0, realizacoes = -1, i, falhas = 0, nreal = 1;

  for (i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--saida") && i + 1 < argc) prefixo = argv[++i];
    else if (!strcmp(argv[i], "--estado")) estado = 1;
    else if (!strcmp(argv[i], "--realizacoes") && i + 1 < argc) realizacoes = atoi(argv[++i]);
    else {
      fprintf(stderr, "uso: %s [--saida prefixo] [--estado] [--realizacoes N]\n", argv[0]);
      return 2;
    }
  }

  if (dae_spec_parse(&S0, DAE_SPEC_JSON, &err) != DAE_OK) {
    fprintf(stderr, "spec embutido invalido em %d:%d: %s\n",
            (int)err.line, (int)err.col, err.msg);
    return 1;
  }
  if (realizacoes > 0) S0.realizations = realizacoes;
  if (S0.core_hash[0] && !dae_spec_hash_confere(&S0)) {
    fprintf(stderr, "aviso: spec gerado pelo nucleo %s, embutido aqui esta o %s\n",
            S0.core_hash, DAE_CORE_HASH);
  }
  fprintf(stderr, "daedalus %s  nucleo %s  %d realizacao(oes)\n",
          DAE_VERSION, DAE_CORE_HASH, (int)S0.realizations);
  nreal = (int)S0.realizations;      /* copiado ANTES de liberar S0 */
  dae_spec_free(&S0);

  /* Cada realizacao analisa o SEU proprio spec. dae_spec tem ponteiros de heap,
     entao compartilhar uma instancia entre threads convidaria a uma liberacao
     dupla; e o parser custa microssegundos contra a propagacao. O nucleo nao
     tem estado global mutavel — e essa regra que torna este laco possivel. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:falhas)
#endif
  for (i = 0; i < nreal; ++i) {
    dae_spec S;
    dae_error e2;
    dae_graph G;
    dae_metrics M;
    dae_series R;
    char nome[512];
    if (dae_spec_parse(&S, DAE_SPEC_JSON, &e2) != DAE_OK) { ++falhas; continue; }
    if (dae_run(&S, S.seed + (uint64_t)i, &G, &M, &R, 0, 0) != DAE_OK) {
      ++falhas;
      dae_spec_free(&S);
      continue;
    }
    if (prefixo && nreal > 1) {
      snprintf(nome, sizeof(nome), "%s.r%04d.csv", prefixo, i);
      escreve(&S, &G, &M, &R, estado, nome);
    } else if (prefixo) {
      escreve(&S, &G, &M, &R, estado, prefixo);
    } else {
#ifdef _OPENMP
#pragma omp critical
#endif
      escreve(&S, &G, &M, &R, estado, 0);
    }
    dae_series_free(&R);
    dae_graph_free(&G);
    dae_spec_free(&S);
  }
  if (falhas) fprintf(stderr, "%d realizacao(oes) falharam\n", falhas);
  return falhas ? 1 : 0;
}
