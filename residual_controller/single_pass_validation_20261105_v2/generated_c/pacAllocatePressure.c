/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacAllocatePressure.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 05-Nov-2026 20:37:30
 */

/* Include Files */
#include "pacAllocatePressure.h"
#include "all.h"
#include "any.h"
#include "minOrMax.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
/*
 * PACALLOCATEPRESSURE Fixed-mean local allocation with calibrated force offset.
 *  forceBias=estimated trim pressure force - nominal static anchor force.
 *  Force-equivalent delivered estimate is pressure force minus this bias.
 *  Uses a nominal max(p-activation,0) map; this is an estimate, not a sensor.
 *
 * Arguments    : const double requested[2]
 *                double cfg_common_bar
 *                double cfg_sTrim_bar
 *                double cfg_aTrim_bar
 *                double cfg_minimum_bar
 *                double cfg_maximum_bar
 *                double cfg_maxDeltaS_bar
 *                double cfg_maxDeltaA_bar
 *                double cfg_activation_bar
 *                double cfg_forcePerBar_N
 *                boolean_T cfg_enableCommonAdaptation
 *                double target[3]
 *                double equivalentForce[2]
 *                double residual[2]
 * Return Type  : boolean_T
 */
boolean_T
pacAllocatePressure(const double requested[2], double cfg_common_bar,
                    double cfg_sTrim_bar, double cfg_aTrim_bar,
                    double cfg_minimum_bar, double cfg_maximum_bar,
                    double cfg_maxDeltaS_bar, double cfg_maxDeltaA_bar,
                    double cfg_activation_bar, double cfg_forcePerBar_N,
                    boolean_T cfg_enableCommonAdaptation, double target[3],
                    double equivalentForce[2], double residual[2])
{
  static const double c_a[20] = {1.0,  -1.0, 0.0, 0.0,  -1.0, 1.0, 0.5,
                                 -0.5, -0.5, 0.5, 0.0,  0.0,  1.0, -1.0,
                                 1.0,  -1.0, 0.5, -0.5, 0.5,  -0.5};
  double b[10];
  double ep[3];
  double offset[3];
  double raw[3];
  double best[2];
  double delta[2];
  double dv1[2];
  double dv2[2];
  double h[2];
  double y[2];
  int b_i;
  int i;
  int j;
  boolean_T limited;
  if (cfg_enableCommonAdaptation) {
    __m128d r;
    __m128d r1;
    double dv[10];
    double a;
    double cost;
    double detA;
    double lo;
    double span;
    boolean_T b_a[10];
    /* PACALLOCATEVARIABLECOMMON Closest feasible nominal force, adjustable
     * mean. */
    /*  Remain on the all-active branch. The feasible pressure-difference
     * polygon */
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
    _mm_storeu_pd(&delta[0], _mm_div_pd(_mm_loadu_pd(&requested[0]),
                                        _mm_set1_pd(cfg_forcePerBar_N)));
    _mm_storeu_pd(&best[0], _mm_set1_pd(0.0));
    cost = rtInf;
    span = delta[0];
    detA = delta[1];
    for (i = 0; i < 10; i++) {
      a = b[i] + 1.0E-10;
      dv[i] = a;
      b_a[i] = (c_a[i] * span + c_a[i + 10] * detA <= a);
    }
    if (b_all(b_a)) {
      best[0] = delta[0];
      best[1] = delta[1];
      cost = 0.0;
    }
    /*  Euclidean projection onto every feasible edge and vertex of the polygon.
     */
    r = _mm_loadu_pd(&delta[0]);
    for (b_i = 0; b_i < 10; b_i++) {
      __m128d r2;
      double ds;
      int c_i;
      ds = c_a[b_i];
      h[0] = ds;
      span = ds * delta[0];
      detA = ds * ds;
      ds = c_a[b_i + 10];
      h[1] = ds;
      span += ds * delta[1];
      detA += ds * ds;
      r1 = _mm_loadu_pd(&h[0]);
      _mm_storeu_pd(
          &h[0],
          _mm_sub_pd(r, _mm_mul_pd(r1, _mm_set1_pd((span - b[b_i]) / detA))));
      span = h[0];
      detA = h[1];
      for (i = 0; i < 10; i++) {
        b_a[i] = (c_a[i] * span + c_a[i + 10] * detA <= dv[i]);
      }
      if (b_all(b_a)) {
        r1 = _mm_loadu_pd(&h[0]);
        r2 = _mm_loadu_pd(&delta[0]);
        r1 = _mm_sub_pd(r1, r2);
        _mm_storeu_pd(&y[0], _mm_mul_pd(r1, r1));
        span = y[0] + y[1];
        if (span < cost) {
          best[0] = h[0];
          best[1] = h[1];
          cost = span;
        }
      }
      c_i = 9 - b_i;
      for (j = 0; j < c_i; j++) {
        int nz;
        nz = (b_i + j) + 1;
        span = c_a[nz + 10];
        detA = c_a[b_i] * span - ds * c_a[nz];
        if (fabs(detA) > 1.0E-12) {
          h[0] = (b[b_i] * span - ds * b[nz]) / detA;
          h[1] = (c_a[b_i] * b[nz] - b[b_i] * c_a[nz]) / detA;
          span = h[0];
          a = h[1];
          for (i = 0; i < 10; i++) {
            b_a[i] = (c_a[i] * span + c_a[i + 10] * a <= dv[i]);
          }
          if (b_all(b_a)) {
            r1 = _mm_loadu_pd(&h[0]);
            r2 = _mm_loadu_pd(&delta[0]);
            r1 = _mm_sub_pd(r1, r2);
            _mm_storeu_pd(&y[0], _mm_mul_pd(r1, r1));
            span = y[0] + y[1];
            if (span < cost) {
              best[0] = h[0];
              best[1] = h[1];
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
    _mm_storeu_pd(&raw[0], _mm_sub_pd(_mm_set1_pd(cfg_maximum_bar), r1));
    ep[2] = lo - offset[2];
    raw[2] = cfg_maximum_bar - offset[2];
    span = fmin(fmax(cfg_common_bar, maximum(ep)), minimum(raw));
    r1 = _mm_loadu_pd(&offset[0]);
    r1 = _mm_add_pd(_mm_set1_pd(span), r1);
    _mm_storeu_pd(&offset[0], r1);
    _mm_storeu_pd(&dv1[0], r1);
    delta[0] = fmin(fmax(dv1[0], lo), cfg_maximum_bar);
    delta[1] = fmin(fmax(dv1[1], lo), cfg_maximum_bar);
    r1 = _mm_loadu_pd(&delta[0]);
    _mm_storeu_pd(&target[0], r1);
    _mm_storeu_pd(&dv2[0], _mm_sub_pd(r1, _mm_set1_pd(cfg_activation_bar)));
    delta[0] = fmax(dv2[0], 0.0);
    delta[1] = fmax(dv2[1], 0.0);
    r1 = _mm_loadu_pd(&delta[0]);
    _mm_storeu_pd(&ep[0], r1);
    span = fmin(fmax(span + offset[2], lo), cfg_maximum_bar);
    target[2] = span;
    equivalentForce[0] = cfg_forcePerBar_N * (ep[1] - ep[0]);
    equivalentForce[1] =
        cfg_forcePerBar_N * (fmax(span - cfg_activation_bar, 0.0) - ep[0]);
    r1 = _mm_loadu_pd(&equivalentForce[0]);
    _mm_storeu_pd(&residual[0], _mm_sub_pd(r1, _mm_loadu_pd(&requested[0])));
    limited = (b_norm(residual) > 1.0E-8);
  } else {
    __m128d r1;
    __m128d r2;
    double a;
    double cost;
    double detA;
    double ds;
    double lo;
    double span;
    int nz;
    boolean_T b_offset[3];
    boolean_T guard1;
    span = 2.0 * cfg_forcePerBar_N;
    lo = (requested[0] + requested[1]) / span;
    detA = (requested[1] - requested[0]) / span;
    ds = fmin(fmax(lo, -cfg_maxDeltaS_bar), cfg_maxDeltaS_bar);
    cost = fmin(fmax(detA, -cfg_maxDeltaA_bar), cfg_maxDeltaA_bar);
    span = cfg_sTrim_bar + ds;
    a = cfg_aTrim_bar + cost;
    raw[0] = cfg_common_bar - 2.0 * span / 3.0;
    span = cfg_common_bar + span / 3.0;
    raw[1] = span - a;
    raw[2] = span + a;
    span = fmin(fmax(raw[0], cfg_minimum_bar), cfg_maximum_bar);
    target[0] = span;
    b_offset[0] = (fabs(raw[0] - span) > 1.0E-10);
    span = fmin(fmax(raw[1], cfg_minimum_bar), cfg_maximum_bar);
    target[1] = span;
    b_offset[1] = (fabs(raw[1] - span) > 1.0E-10);
    span = fmin(fmax(raw[2], cfg_minimum_bar), cfg_maximum_bar);
    target[2] = span;
    b_offset[2] = (fabs(raw[2] - span) > 1.0E-10);
    if (any(b_offset) || (fabs(ds - lo) > 1.0E-10) ||
        (fabs(cost - detA) > 1.0E-10)) {
      limited = true;
    } else {
      limited = false;
    }
    /*  Invert the same activation map used below. Effective pressures must
     * satisfy */
    /*  e2-e1 and e3-e1 = (requested+forceBias)/forcePerBar. Retain fixed mean;
     */
    /*  inactive muscles can take up the remaining pressure below activation. */
    offset[0] = 0.0;
    span = requested[0] / cfg_forcePerBar_N;
    offset[1] = span;
    detA = span;
    span = requested[1] / cfg_forcePerBar_N;
    offset[2] = span;
    detA += span;
    span = fmax(cfg_common_bar - detA / 3.0, -minimum(offset));
    r1 = _mm_loadu_pd(&offset[0]);
    r1 = _mm_add_pd(_mm_set1_pd(span), r1);
    _mm_storeu_pd(&offset[0], r1);
    _mm_storeu_pd(&delta[0], r1);
    delta[0] = fmax(delta[0], 0.0);
    delta[1] = fmax(delta[1], 0.0);
    r1 = _mm_loadu_pd(&delta[0]);
    _mm_storeu_pd(&ep[0], r1);
    r2 = _mm_set1_pd(cfg_activation_bar);
    _mm_storeu_pd(&raw[0], _mm_add_pd(r2, r1));
    span = fmax(span + offset[2], 0.0);
    raw[2] = cfg_activation_bar + span;
    nz = ((ep[0] <= 1.0E-10) + (ep[1] <= 1.0E-10)) + (span <= 1.0E-10);
    detA = ((raw[0] + raw[1]) + raw[2]) - 3.0 * cfg_common_bar;
    guard1 = false;
    if (detA > 1.0E-9) {
      if (nz > 0) {
        if (ep[0] <= 1.0E-10) {
          raw[0] -= detA / (double)nz;
        }
        if (ep[1] <= 1.0E-10) {
          raw[1] -= detA / (double)nz;
        }
        if (span <= 1.0E-10) {
          raw[2] -= detA / (double)nz;
        }
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
    if (guard1) {
      b_offset[0] = (raw[0] >= cfg_minimum_bar - 1.0E-10);
      b_offset[1] = (raw[1] >= cfg_minimum_bar - 1.0E-10);
      b_offset[2] = (raw[2] >= cfg_minimum_bar - 1.0E-10);
      if (c_all(b_offset)) {
        b_offset[0] = (raw[0] <= cfg_maximum_bar + 1.0E-10);
        b_offset[1] = (raw[1] <= cfg_maximum_bar + 1.0E-10);
        b_offset[2] = (raw[2] <= cfg_maximum_bar + 1.0E-10);
        if (c_all(b_offset) &&
            (fabs(((raw[1] + raw[2]) / 2.0 - raw[0]) - cfg_sTrim_bar) <=
             cfg_maxDeltaS_bar + 1.0E-10) &&
            (fabs((raw[2] - raw[1]) / 2.0 - cfg_aTrim_bar) <=
             cfg_maxDeltaA_bar + 1.0E-10)) {
          target[0] = fmin(fmax(raw[0], cfg_minimum_bar), cfg_maximum_bar);
          target[1] = fmin(fmax(raw[1], cfg_minimum_bar), cfg_maximum_bar);
          target[2] = fmin(fmax(raw[2], cfg_minimum_bar), cfg_maximum_bar);
          limited = false;
        }
      }
    }
    r1 = _mm_loadu_pd(&target[0]);
    _mm_storeu_pd(&best[0], _mm_sub_pd(r1, r2));
    delta[0] = fmax(best[0], 0.0);
    delta[1] = fmax(best[1], 0.0);
    r1 = _mm_loadu_pd(&delta[0]);
    _mm_storeu_pd(&offset[0], r1);
    equivalentForce[0] = cfg_forcePerBar_N * (offset[1] - offset[0]);
    equivalentForce[1] =
        cfg_forcePerBar_N *
        (fmax(target[2] - cfg_activation_bar, 0.0) - offset[0]);
    r1 = _mm_loadu_pd(&equivalentForce[0]);
    _mm_storeu_pd(&residual[0], _mm_sub_pd(r1, _mm_loadu_pd(&requested[0])));
    /*  A rail/correction/activation-infeasible request is a constrained
     * allocation, */
    /*  even when the legacy bounded fallback happens not to touch a pressure
     * rail. */
    if (limited || (b_norm(residual) > 1.0E-8)) {
      limited = true;
    }
  }
  return limited;
}

/*
 * File trailer for pacAllocatePressure.c
 *
 * [EOF]
 */
