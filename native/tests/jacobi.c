#include "jacobi.h"

#include <math.h>

int dae_jacobi(double *a, double *v, double *w, int32_t n)
{
  int32_t i, j, p, q, sweep;
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j) v[i * n + j] = (i == j) ? 1.0 : 0.0;

  for (sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (p = 0; p < n; ++p)
      for (q = p + 1; q < n; ++q) off += a[p * n + q] * a[p * n + q];
    if (off <= 1e-30) break;

    for (p = 0; p < n; ++p) {
      for (q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        double theta, t, c, s;
        if (fabs(apq) < 1e-300) continue;
        theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
        c = 1.0 / sqrt(t * t + 1.0);
        s = t * c;
        for (i = 0; i < n; ++i) {
          const double aip = a[i * n + p], aiq = a[i * n + q];
          a[i * n + p] = c * aip - s * aiq;
          a[i * n + q] = s * aip + c * aiq;
        }
        for (i = 0; i < n; ++i) {
          const double api = a[p * n + i], aqi = a[q * n + i];
          a[p * n + i] = c * api - s * aqi;
          a[q * n + i] = s * api + c * aqi;
        }
        for (i = 0; i < n; ++i) {
          const double vip = v[i * n + p], viq = v[i * n + q];
          v[i * n + p] = c * vip - s * viq;
          v[i * n + q] = s * vip + c * viq;
        }
      }
    }
  }
  for (i = 0; i < n; ++i) w[i] = a[i * n + i];
  return (sweep >= 100) ? 1 : 0;
}
