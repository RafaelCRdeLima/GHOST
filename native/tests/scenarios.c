#include "scenarios.h"

#include <string.h>

static dae_scenario TAB[8];
static int PRONTO = 0;

static void monta(void)
{
  int i;
  for (i = 0; i < 8; ++i) {
    memset(&TAB[i], 0, sizeof(TAB[i]));
    dae_gen_params_default(&TAB[i].gen);
    TAB[i].ham = DAE_H_ADJACENCY;
    TAB[i].norm = DAE_NORM_NONE;
    TAB[i].gamma = 1.0;
    TAB[i].t1 = 10.0;
    TAB[i].nt = 40;
    TAB[i].lanczos = 0;
    TAB[i].init_site = 0;
    TAB[i].target = -1;
  }

  TAB[0].nome = "microtubulo 40x13 seam=3";
  TAB[0].gen.kind = DAE_G_MICROTUBULE;
  TAB[0].gen.n_par = 40; TAB[0].gen.n_perp = 13; TAB[0].gen.seam_shift = 3;
  TAB[0].gen.n_modules = 4; TAB[0].gen.seed = 2026ULL;
  TAB[0].init_site = 6 * 13 + 6; TAB[0].target = 33 * 13 + 2;

  TAB[1].nome = "microtubulo religado p=0.15";
  TAB[1].gen.kind = DAE_G_MICROTUBULE;
  TAB[1].gen.n_par = 60; TAB[1].gen.n_perp = 13; TAB[1].gen.seam_shift = 3;
  TAB[1].gen.n_modules = 6; TAB[1].gen.ws_p = 0.15;
  TAB[1].gen.conn_mode = DAE_REWIRE; TAB[1].gen.seed = 4242ULL;
  TAB[1].init_site = 0; TAB[1].target = 59 * 13 + 12;
  TAB[1].norm = DAE_NORM_SPECTRAL; TAB[1].lanczos = 40;

  TAB[2].nome = "microtubulo laplaciana, grau medio";
  TAB[2].gen.kind = DAE_G_MICROTUBULE;
  TAB[2].gen.n_par = 30; TAB[2].gen.n_perp = 13; TAB[2].gen.seam_shift = 0;
  TAB[2].gen.n_modules = 3; TAB[2].gen.seed = 7ULL;
  TAB[2].ham = DAE_H_LAPLACIAN; TAB[2].norm = DAE_NORM_MEAN_DEGREE;
  TAB[2].init_site = 15 * 13; TAB[2].target = 4;

  TAB[3].nome = "SBM modular";
  TAB[3].gen.kind = DAE_G_SBM;
  TAB[3].gen.n = 300; TAB[3].gen.n_modules = 4;
  TAB[3].gen.p_in = 0.3; TAB[3].gen.p_out = 0.01; TAB[3].gen.seed = 99ULL;
  TAB[3].target = 299; TAB[3].t1 = 5.0;

  TAB[4].nome = "linha longa, dt grande";
  TAB[4].gen.kind = DAE_G_PATH; TAB[4].gen.n = 801;
  TAB[4].init_site = 400; TAB[4].t1 = 60.0; TAB[4].nt = 12;

  TAB[5].nome = "ciclo";
  TAB[5].gen.kind = DAE_G_CYCLE; TAB[5].gen.n = 256;
  TAB[5].target = 128; TAB[5].t1 = 20.0;

  TAB[6].nome = "hipercubo Q_9";
  TAB[6].gen.kind = DAE_G_HYPERCUBE; TAB[6].gen.dim = 9;
  TAB[6].target = 511; TAB[6].t1 = 8.0;
  TAB[6].norm = DAE_NORM_SPECTRAL; TAB[6].lanczos = 30;

  TAB[7].nome = "grade 2D acrescentada p=0.3";
  TAB[7].gen.kind = DAE_G_GRID2D; TAB[7].gen.rows = 20; TAB[7].gen.cols = 25;
  TAB[7].gen.ws_p = 0.3; TAB[7].gen.conn_mode = DAE_ADD; TAB[7].gen.seed = 31337ULL;
  TAB[7].target = 499; TAB[7].t1 = 15.0;

  PRONTO = 1;
}

int dae_scenario_count(void) { return 8; }

const dae_scenario *dae_scenario_get(int i)
{
  if (!PRONTO) monta();
  if (i < 0 || i >= 8) return 0;
  return &TAB[i];
}
