#include "dae_obs.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

dae_status dae_obs_alloc(dae_obs *O, const dae_obs_cfg *c, int32_t n)
{
  int32_t nmod;
  if (!O || !c || n <= 0) return DAE_ERR_PARAM;
  memset(O, 0, sizeof(*O));
  nmod = (c->module_of && c->nmod > 0) ? c->nmod : 1;
  O->n = n;
  O->nmod = nmod;

  if (c->want_pop) {
    O->pop = (double *)malloc((size_t)n * sizeof(double));
    if (!O->pop) { dae_obs_free(O); return DAE_ERR_ALLOC; }
  }
  O->pmod  = (double *)malloc((size_t)nmod * sizeof(double));
  O->l1mod = (double *)malloc((size_t)nmod * sizeof(double));
  if (!O->pmod || !O->l1mod) { dae_obs_free(O); return DAE_ERR_ALLOC; }

  if (c->want_conc_mod) {
    O->conc_mod = (double *)malloc((size_t)nmod * (size_t)nmod * sizeof(double));
    if (!O->conc_mod) { dae_obs_free(O); return DAE_ERR_ALLOC; }
  }
  if (c->want_conc_full) {
    if (n > DAE_FULL_CONC_MAX) { dae_obs_free(O); return DAE_ERR_TOOBIG; }
    O->conc_full = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!O->conc_full) { dae_obs_free(O); return DAE_ERR_ALLOC; }
  }
  return DAE_OK;
}

void dae_obs_free(dae_obs *O)
{
  if (!O) return;
  free(O->pop); free(O->pmod); free(O->l1mod);
  free(O->conc_mod); free(O->conc_full);
  memset(O, 0, sizeof(*O));
}

dae_status dae_obs_eval(const double *psire, const double *psiim, int32_t n,
                        const dae_obs_cfg *c, dae_obs *O)
{
  const int32_t nmod = O ? O->nmod : 0;
  int32_t i, j, m;
  double s1 = 0.0, s2 = 0.0, s4 = 0.0;

  if (!psire || !psiim || !c || !O || n != O->n) return DAE_ERR_PARAM;
  if (c->target >= n) return DAE_ERR_PARAM;

  for (m = 0; m < nmod; ++m) { O->pmod[m] = 0.0; O->l1mod[m] = 0.0; }

  for (i = 0; i < n; ++i) {
    const double p = psire[i] * psire[i] + psiim[i] * psiim[i];
    const double a = sqrt(p);
    s1 += a; s2 += p; s4 += p * p;
    if (O->pop) O->pop[i] = p;
    if (c->module_of) {
      const int32_t mi = c->module_of[i];
      if (mi < 0 || mi >= nmod) return DAE_ERR_PARAM;
      O->pmod[mi]  += p;
      O->l1mod[mi] += a;
    } else {
      O->pmod[0]  += p;
      O->l1mod[0] += a;
    }
  }

  O->norm   = s2;
  O->ipr    = s4;
  O->coh_l1 = s1 * s1 - 1.0;
  if (c->target >= 0) {
    O->p_target = psire[c->target] * psire[c->target] +
                  psiim[c->target] * psiim[c->target];
  } else {
    O->p_target = NAN;      /* ausência de alvo não é p_alvo = 0 */
  }

  if (O->conc_mod) {
    for (m = 0; m < nmod; ++m) {
      for (j = 0; j < nmod; ++j) {
        O->conc_mod[(size_t)m * (size_t)nmod + (size_t)j] =
          (m == j) ? (O->l1mod[m] * O->l1mod[m] - O->pmod[m])
                   : (2.0 * O->l1mod[m] * O->l1mod[j]);
      }
    }
  }

  if (O->conc_full) {
    for (i = 0; i < n; ++i) {
      const double ai = sqrt(psire[i] * psire[i] + psiim[i] * psiim[i]);
      for (j = 0; j < n; ++j) {
        const double aj = sqrt(psire[j] * psire[j] + psiim[j] * psiim[j]);
        O->conc_full[(size_t)i * (size_t)n + (size_t)j] = (i == j) ? 0.0 : 2.0 * ai * aj;
      }
    }
  }
  return DAE_OK;
}
