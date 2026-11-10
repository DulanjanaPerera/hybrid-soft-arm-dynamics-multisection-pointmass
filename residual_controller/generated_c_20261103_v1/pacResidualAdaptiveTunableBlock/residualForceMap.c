/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualForceMap.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "residualForceMap.h"
#include "pacResidualAdaptiveTunableBlock_data.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * RESIDUALFORCEMAP Activation re-expression keeps F_cfg+r_eff=B0*p+r0+deltaG.
 *
 * Arguments    : const double p[3]
 *                const double measured[6]
 *                double dt
 *                const double h[128]
 *                const double network_input[3072]
 *                const double network_recurrent[49152]
 *                const double network_bias[384]
 *                const double network_output[256]
 *                const double network_skip[16]
 *                const double network_outputBias[2]
 *                const double network_mean[8]
 *                const double network_std[8]
 *                double cfg_activation_bar
 *                double cfg_forcePerBar_N
 *                const double gravityShift[2]
 *                double total[2]
 *                double correction[2]
 * Return Type  : void
 */
void residualForceMap(
    const double p[3], const double measured[6], double dt, const double h[128],
    const double network_input[3072], const double network_recurrent[49152],
    const double network_bias[384], const double network_output[256],
    const double network_skip[16], const double network_outputBias[2],
    const double network_mean[8], const double network_std[8],
    double cfg_activation_bar, double cfg_forcePerBar_N,
    const double gravityShift[2], double total[2], double correction[2])
{
  __m128d r;
  double u[384];
  double w[384];
  double z[8];
  double ep[3];
  double a[2];
  double b_dv[2];
  double b_network_skip[2];
  double d;
  double d1;
  int i;
  int i1;
  int i2;
  int k;
  /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
   * evaluation. */
  /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history fixed. */
  _mm_storeu_pd(&z[0], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                             _mm_loadu_pd(&network_mean[0])),
                                  _mm_loadu_pd(&network_std[0])));
  _mm_storeu_pd(&z[2], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                             _mm_loadu_pd(&network_mean[2])),
                                  _mm_loadu_pd(&network_std[2])));
  _mm_storeu_pd(&z[4], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&p[0]),
                                             _mm_loadu_pd(&network_mean[4])),
                                  _mm_loadu_pd(&network_std[4])));
  z[6] = (p[2] - network_mean[6]) / network_std[6];
  z[7] = (dt - network_mean[7]) / network_std[7];
#pragma omp parallel for num_threads(omp_get_max_threads()) private(d, i2)

  for (i = 0; i < 384; i++) {
    d = 0.0;
    for (i2 = 0; i2 < 8; i2++) {
      d += network_input[i + 384 * i2] * z[i2];
    }
    u[i] = d + network_bias[i];
  }
  memset(&w[0], 0, 384U * sizeof(double));
  for (k = 0; k < 128; k++) {
    for (i1 = 0; i1 <= 382; i1 += 2) {
      r = _mm_loadu_pd(&w[i1]);
      _mm_storeu_pd(
          &w[i1],
          _mm_add_pd(r,
                     _mm_mul_pd(_mm_loadu_pd(&network_recurrent[i1 + 384 * k]),
                                _mm_set1_pd(h[k]))));
    }
  }
  memset(&a[0], 0, sizeof(double) << 1);
  for (k = 0; k < 128; k++) {
    d1 = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
    d1 = d1 * h[k] + (1.0 - d1) * tanh(u[k + 256] +
                                       0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0) *
                                           w[k + 256]);
    r = _mm_loadu_pd(&a[0]);
    _mm_storeu_pd(
        &a[0], _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network_output[k << 1]),
                                        _mm_set1_pd(d1))));
  }
  memset(&b_network_skip[0], 0, sizeof(double) << 1);
  for (k = 0; k < 8; k++) {
    r = _mm_loadu_pd(&b_network_skip[0]);
    _mm_storeu_pd(&b_network_skip[0],
                  _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network_skip[k << 1]),
                                           _mm_set1_pd(z[k]))));
  }
  __m128d r1;
  double d2;
  double d3;
  a[0] =
      tanh(10.0 * ((a[0] + b_network_skip[0]) + network_outputBias[0]) / 80.0);
  a[1] =
      tanh(10.0 * ((a[1] + b_network_skip[1]) + network_outputBias[1]) / 80.0);
  _mm_storeu_pd(&b_dv[0], _mm_sub_pd(_mm_loadu_pd(&p[0]),
                                     _mm_set1_pd(cfg_activation_bar)));
  b_network_skip[0] = fmax(b_dv[0], 0.0);
  b_network_skip[1] = fmax(b_dv[1], 0.0);
  r = _mm_loadu_pd(&b_network_skip[0]);
  _mm_storeu_pd(&ep[0], r);
  d1 = p[0];
  d2 = p[1];
  d3 = p[2];
  r = _mm_loadu_pd(&dv[2]);
  r = _mm_mul_pd(r, _mm_set1_pd(d2));
  r1 = _mm_add_pd(_mm_set1_pd(-13.273228961416876 * d1), r);
  r = _mm_loadu_pd(&dv[4]);
  r = _mm_mul_pd(r, _mm_set1_pd(d3));
  r1 = _mm_add_pd(r1, r);
  r = _mm_loadu_pd(&a[0]);
  r = _mm_mul_pd(_mm_set1_pd(80.0), r);
  r = _mm_add_pd(r1, r);
  r1 = _mm_loadu_pd(&gravityShift[0]);
  r = _mm_add_pd(r, r1);
  _mm_storeu_pd(&total[0], r);
  correction[0] = total[0] - cfg_forcePerBar_N * (ep[1] - ep[0]);
  correction[1] = total[1] - cfg_forcePerBar_N *
                                 (fmax(p[2] - cfg_activation_bar, 0.0) - ep[0]);
}

/*
 * File trailer for residualForceMap.c
 *
 * [EOF]
 */
