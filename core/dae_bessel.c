#include "dae_bessel.h"

#include <math.h>

dae_status dae_bessel_j_array(double x, int32_t kmax, double *out)
{
  /* A escada de reescala é de 100 décadas, não 250: quanto menor o degrau,
     menos faixa dinâmica se perde nos termos já armazenados. */
  const double BIG = 1e100, SMALL = 1e-100;
  int32_t start, k, i;
  double jk, jkp1, sum;

  if (!out || kmax < 0 || !(x == x) || x < 0.0) return DAE_ERR_PARAM;
  for (i = 0; i <= kmax; ++i) out[i] = 0.0;
  if (x == 0.0) { out[0] = 1.0; return DAE_OK; }

  {
    const double top = (x > (double)kmax) ? x : (double)kmax;
    const double s = top + 40.0 + 10.0 * sqrt(top + 1.0);
    if (s > 1.0e8) return DAE_ERR_TOOBIG;
    start = (int32_t)s;
    if (start & 1) ++start;          /* par: a soma de normalização é sobre pares */
  }

  jkp1 = 0.0;                        /* j_{start+1} */
  jk   = SMALL;                      /* j_{start}   */
  sum  = 0.0;
  for (k = start; k >= 1; --k) {
    const double jm = (2.0 * (double)k / x) * jk - jkp1;   /* j_{k-1} */
    jkp1 = jk;
    jk   = jm;
    if (fabs(jk) > BIG) {
      jk *= SMALL; jkp1 *= SMALL; sum *= SMALL;
      for (i = k; i <= kmax; ++i) out[i] *= SMALL;
    }
    {
      const int32_t idx = k - 1;
      if (idx <= kmax)         out[idx] = jk;
      if (idx == 0)            sum += jk;
      else if ((idx & 1) == 0) sum += 2.0 * jk;
    }
  }
  if (!(sum == sum) || sum == 0.0) return DAE_ERR_PARAM;
  for (i = 0; i <= kmax; ++i) out[i] /= sum;
  return DAE_OK;
}
