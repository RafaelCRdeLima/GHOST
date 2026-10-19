/* t95_layout.c — o layout espectral separa? E, quando NAO deveria, nao separa?
 *
 * A asserção "os modulos ficam separados no plano" passa trivialmente numa
 * implementacao que empilhasse os vertices por indice, porque o SBM numera os
 * modulos em blocos contiguos. Por isso ela vem sempre em par com o caso
 * homogeneo — p_in = p_out, sem modulo nenhum — onde a MESMA medida tem de
 * dizer que nao ha separacao. Sem o par, a fixture nao testa nada.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

/* Razao entre a distancia media ENTRE modulos e a distancia media DENTRO de um
   modulo, no plano do embedding. Grande = separou. */
static double razao_separacao(const dae_graph *G, const float *xy)
{
  double dentro = 0.0, entre = 0.0;
  int64_t nd = 0, ne = 0;
  int32_t i, j;
  for (i = 0; i < G->n; ++i) {
    for (j = i + 1; j < G->n; ++j) {
      const double dx = (double)xy[2 * i] - (double)xy[2 * j];
      const double dy = (double)xy[2 * i + 1] - (double)xy[2 * j + 1];
      const double d = sqrt(dx * dx + dy * dy);
      if (G->module_of[i] == G->module_of[j]) { dentro += d; ++nd; }
      else                                    { entre += d; ++ne; }
    }
  }
  if (nd == 0 || ne == 0) return 1.0;
  dentro /= (double)nd;
  entre  /= (double)ne;
  return (dentro > 1e-14) ? entre / dentro : 1e9;
}

static double caso(dae_test *T, const char *nome, double p_in, double p_out)
{
  dae_gen_params p;
  dae_graph G;
  float *xy;
  double razao;

  dae_gen_params_default(&p);
  p.kind = DAE_G_SBM; p.n = 240; p.n_modules = 4;
  p.p_in = p_in; p.p_out = p_out; p.seed = 2026ULL;
  if (dae_graph_build(&G, &p) != DAE_OK) { dae_test_ok(T, 0, "%s: montagem", nome); return 1.0; }

  xy = (float *)calloc(2u * (size_t)G.n, sizeof(float));
  dae_test_ok(T, dae_layout_espectral(&G, xy) == DAE_OK, "%s: layout", nome);
  razao = razao_separacao(&G, xy);

  { /* anti-vacuidade basica: as coordenadas nao sao todas iguais */
    float lo = xy[0], hi = xy[0];
    int32_t i;
    for (i = 0; i < 2 * G.n; ++i) { if (xy[i] < lo) lo = xy[i]; if (xy[i] > hi) hi = xy[i]; }
    dae_test_ok(T, hi - lo > 1e-6, "%s: o embedding nao colapsou num ponto (%g)",
                nome, (double)(hi - lo));
  }

  dae_test_note("%-26s p_in=%.2f p_out=%.3f  entre/dentro = %.3f", nome, p_in, p_out, razao);
  free(xy);
  dae_graph_free(&G);
  return razao;
}

int main(void)
{
  dae_test T;
  double modular, homogeneo;
  dae_test_begin(&T, "layout espectral: separa quando ha modulo, e so quando ha");

  modular   = caso(&T, "SBM modular", 0.35, 0.004);
  homogeneo = caso(&T, "SBM homogeneo", 0.10, 0.10);

  dae_test_ok(&T, modular > 2.0,
              "grafo modular: os modulos se separam no plano (razao %.3f, esperado > 2)",
              modular);

  /* A COMPANHEIRA. Com p_in = p_out nao existe modulo: a particao continua
     rotulando os vertices, mas nao corresponde a estrutura nenhuma. Se a
     medida acusasse separacao aqui, ela estaria medindo a numeracao dos
     vertices, nao a topologia — e a assercao acima seria vazia. */
  dae_test_ok(&T, homogeneo < 1.3,
              "grafo homogeneo: NAO ha separacao (razao %.3f, esperado < 1.3)",
              homogeneo);

  dae_test_ok(&T, modular > 3.0 * homogeneo,
              "a separacao distingue os dois casos por mais de 3x (%.3f contra %.3f)",
              modular, homogeneo);

  /* Grafo com geometria propria continua com a dela: o layout espectral e o
     recurso de QUEM NAO TEM, nao um substituto universal. */
  {
    dae_gen_params p;
    dae_graph G;
    dae_gen_params_default(&p);
    p.kind = DAE_G_MICROTUBULE; p.n_par = 30; p.n_perp = 13; p.n_modules = 3;
    dae_graph_build(&G, &p);
    dae_test_ok(&T, G.n_par == 30 && G.n_perp == 13,
                "microtubulo mantem a geometria desenrolada (%d x %d)",
                (int)G.n_par, (int)G.n_perp);
    dae_test_ok(&T, G.xy[2 * (5 * 13 + 7)] == 5.0f && G.xy[2 * (5 * 13 + 7) + 1] == 7.0f,
                "e as coordenadas dele sao (m, q), nao o embedding");
    dae_graph_free(&G);
  }

  return dae_test_end(&T);
}
