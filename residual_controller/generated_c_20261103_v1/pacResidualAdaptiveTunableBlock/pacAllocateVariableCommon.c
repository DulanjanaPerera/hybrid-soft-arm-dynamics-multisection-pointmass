/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacAllocateVariableCommon.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "pacAllocateVariableCommon.h"
#include "all.h"
#include "minOrMax.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
/*
 * PACALLOCATEVARIABLECOMMON Closest feasible nominal force, adjustable mean.
 *
 * Arguments    : const double requested[2]
 *                const double bias[2]
 *                double cfg_common_bar
 *                double cfg_sTrim_bar
 *                double cfg_aTrim_bar
 *                double cfg_minimum_bar
 *                double cfg_maximum_bar
 *                double cfg_maxDeltaS_bar
 *                double cfg_maxDeltaA_bar
 *                double cfg_activation_bar
 *                double cfg_forcePerBar_N
 *                double target[3]
 *                double equivalent[2]
 *                double residual[2]
 * Return Type  : boolean_T
 */
boolean_T pacAllocateVariableCommon(
    const double requested[2], const double bias[2], double cfg_common_bar,
    double cfg_sTrim_bar, double cfg_aTrim_bar, double cfg_minimum_bar,
    double cfg_maximum_bar, double cfg_maxDeltaS_bar, double cfg_maxDeltaA_bar,
    double cfg_activation_bar, double cfg_forcePerBar_N, double target[3],
    double equivalent[2], double residual[2])
{
  static const double b_a[20] = {1.0,  -1.0, 0.0, 0.0,  -1.0, 1.0, 0.5,
                                 -0.5, -0.5, 0.5, 0.0,  0.0,  1.0, -1.0,
                                 1.0,  -1.0, 0.5, -0.5, 0.5,  -0.5};
  __m128d r;
  __m128d r1;
  double b[10];
  double b_dv[10];
  double b_cfg_maximum_bar[3];
  double ep[3];
  double offset[3];
  double b_h[2];
  double best[2];
  double d[2];
  double dv1[2];
  double dv2[2];
  double y[2];
  double cost;
  double detA;
  double h;
  double lo;
  double span;
  int b_i;
  int i;
  int j;
  boolean_T a[10];
  /*  Remain on the all-active branch. The feasible pressure-difference polygon
   */
  /*  enforces rails and trim-relative s/a clamps. No optimizer toolbox
   * required. */
  lo = fmax(cfg_minimum_bar, cfg_activation_bar);
  span = cfg_maximum_bar - lo;
  b[0] = span;
  b[1] = span;
  b[2] = span;
  b[3] = span;
  b[4] = span;
  b[5] = span;
  b[6] = cfg_sTrim_bar + cfg_maxDeltaS_bar;
  b[7] = -cfg_sTrim_bar + cfg_maxDeltaS_bar;
  b[8] = cfg_aTrim_bar + cfg_maxDeltaA_bar;
  b[9] = -cfg_aTrim_bar + cfg_maxDeltaA_bar;
  _mm_storeu_pd(&d[0], _mm_div_pd(_mm_add_pd(_mm_loadu_pd(&requested[0]),
                                             _mm_loadu_pd(&bias[0])),
                                  _mm_set1_pd(cfg_forcePerBar_N)));
  _mm_storeu_pd(&best[0], _mm_set1_pd(0.0));
  cost = rtInf;
  span = d[0];
  h = d[1];
  for (i = 0; i < 10; i++) {
    detA = b[i] + 1.0E-10;
    b_dv[i] = detA;
    a[i] = (b_a[i] * span + b_a[i + 10] * h <= detA);
  }
  if (b_all(a)) {
    best[0] = d[0];
    best[1] = d[1];
    cost = 0.0;
  }
  /*  Euclidean projection onto every feasible edge and vertex of the polygon.
   */
  r = _mm_loadu_pd(&d[0]);
  for (b_i = 0; b_i < 10; b_i++) {
    __m128d r2;
    double b_d;
    int c_i;
    b_d = b_a[b_i];
    b_h[0] = b_d;
    span = b_d * d[0];
    h = b_d * b_d;
    b_d = b_a[b_i + 10];
    b_h[1] = b_d;
    span += b_d * d[1];
    h += b_d * b_d;
    r1 = _mm_loadu_pd(&b_h[0]);
    _mm_storeu_pd(
        &b_h[0],
        _mm_sub_pd(r, _mm_mul_pd(r1, _mm_set1_pd((span - b[b_i]) / h))));
    span = b_h[0];
    detA = b_h[1];
    for (i = 0; i < 10; i++) {
      a[i] = (b_a[i] * span + b_a[i + 10] * detA <= b_dv[i]);
    }
    if (b_all(a)) {
      r1 = _mm_loadu_pd(&b_h[0]);
      r2 = _mm_loadu_pd(&d[0]);
      r1 = _mm_sub_pd(r1, r2);
      _mm_storeu_pd(&y[0], _mm_mul_pd(r1, r1));
      span = y[0] + y[1];
      if (span < cost) {
        best[0] = b_h[0];
        best[1] = b_h[1];
        cost = span;
      }
    }
    c_i = 9 - b_i;
    for (j = 0; j < c_i; j++) {
      int b_j;
      b_j = (b_i + j) + 1;
      span = b_a[b_j + 10];
      detA = b_a[b_i] * span - b_d * b_a[b_j];
      if (fabs(detA) > 1.0E-12) {
        b_h[0] = (b[b_i] * span - b_d * b[b_j]) / detA;
        b_h[1] = (b_a[b_i] * b[b_j] - b[b_i] * b_a[b_j]) / detA;
        span = b_h[0];
        h = b_h[1];
        for (i = 0; i < 10; i++) {
          a[i] = (b_a[i] * span + b_a[i + 10] * h <= b_dv[i]);
        }
        if (b_all(a)) {
          r1 = _mm_loadu_pd(&b_h[0]);
          r2 = _mm_loadu_pd(&d[0]);
          r1 = _mm_sub_pd(r1, r2);
          _mm_storeu_pd(&y[0], _mm_mul_pd(r1, r1));
          span = y[0] + y[1];
          if (span < cost) {
            best[0] = b_h[0];
            best[1] = b_h[1];
            cost = span;
          }
        }
      }
    }
  }
  offset[0] = -(best[0] + best[1]) / 3.0;
  offset[1] = (2.0 * best[0] - best[1]) / 3.0;
  offset[2] = (2.0 * best[1] - best[0]) / 3.0;
  r1 = _mm_loadu_pd(&offset[0]);
  _mm_storeu_pd(&ep[0], _mm_sub_pd(_mm_set1_pd(lo), r1));
  _mm_storeu_pd(&b_cfg_maximum_bar[0],
                _mm_sub_pd(_mm_set1_pd(cfg_maximum_bar), r1));
  ep[2] = lo - offset[2];
  b_cfg_maximum_bar[2] = cfg_maximum_bar - offset[2];
  span = fmin(fmax(cfg_common_bar, maximum(ep)), minimum(b_cfg_maximum_bar));
  r1 = _mm_loadu_pd(&offset[0]);
  _mm_storeu_pd(&dv1[0], _mm_add_pd(_mm_set1_pd(span), r1));
  d[0] = fmin(fmax(dv1[0], lo), cfg_maximum_bar);
  d[1] = fmin(fmax(dv1[1], lo), cfg_maximum_bar);
  r1 = _mm_loadu_pd(&d[0]);
  _mm_storeu_pd(&target[0], r1);
  r1 = _mm_sub_pd(r1, _mm_set1_pd(cfg_activation_bar));
  _mm_storeu_pd(&offset[0], r1);
  _mm_storeu_pd(&dv2[0], r1);
  d[0] = fmax(dv2[0], 0.0);
  d[1] = fmax(dv2[1], 0.0);
  r1 = _mm_loadu_pd(&d[0]);
  _mm_storeu_pd(&ep[0], r1);
  span = fmin(fmax(span + offset[2], lo), cfg_maximum_bar);
  target[2] = span;
  equivalent[0] = cfg_forcePerBar_N * (ep[1] - ep[0]) - bias[0];
  equivalent[1] =
      cfg_forcePerBar_N * (fmax(span - cfg_activation_bar, 0.0) - ep[0]) -
      bias[1];
  r1 = _mm_loadu_pd(&equivalent[0]);
  _mm_storeu_pd(&residual[0], _mm_sub_pd(r1, _mm_loadu_pd(&requested[0])));
  return b_norm(residual) > 1.0E-8;
}

/*
 * File trailer for pacAllocateVariableCommon.c
 *
 * [EOF]
 */
