/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacAllocatePressure.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "pacAllocatePressure.h"
#include "all.h"
#include "allOrAny.h"
#include "minOrMax.h"
#include "norm.h"
#include "pacAllocateVariableCommon.h"
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
 *                const double anchorForce[2]
 *                const double forceBias[2]
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
pacAllocatePressure(const double requested[2], const double anchorForce[2],
                    const double forceBias[2], double cfg_common_bar,
                    double cfg_sTrim_bar, double cfg_aTrim_bar,
                    double cfg_minimum_bar, double cfg_maximum_bar,
                    double cfg_maxDeltaS_bar, double cfg_maxDeltaA_bar,
                    double cfg_activation_bar, double cfg_forcePerBar_N,
                    boolean_T cfg_enableCommonAdaptation, double target[3],
                    double equivalentForce[2], double residual[2])
{
  double b_candidate[3];
  double offset[3];
  double raw[3];
  double b_dv[2];
  double delta[2];
  double dv1[2];
  boolean_T limited;
  if (cfg_enableCommonAdaptation) {
    limited = pacAllocateVariableCommon(
        requested, forceBias, cfg_common_bar, cfg_sTrim_bar, cfg_aTrim_bar,
        cfg_minimum_bar, cfg_maximum_bar, cfg_maxDeltaS_bar, cfg_maxDeltaA_bar,
        cfg_activation_bar, cfg_forcePerBar_N, target, equivalentForce,
        residual);
  } else {
    __m128d r;
    __m128d r1;
    double a;
    double da;
    double daWanted;
    double ds;
    double dsWanted;
    double excess;
    double s;
    int nz;
    boolean_T offset_data[10];
    boolean_T guard1;
    _mm_storeu_pd(&delta[0], _mm_sub_pd(_mm_loadu_pd(&requested[0]),
                                        _mm_loadu_pd(&anchorForce[0])));
    s = 2.0 * cfg_forcePerBar_N;
    dsWanted = (delta[0] + delta[1]) / s;
    daWanted = (delta[1] - delta[0]) / s;
    ds = fmin(fmax(dsWanted, -cfg_maxDeltaS_bar), cfg_maxDeltaS_bar);
    da = fmin(fmax(daWanted, -cfg_maxDeltaA_bar), cfg_maxDeltaA_bar);
    s = cfg_sTrim_bar + ds;
    a = cfg_aTrim_bar + da;
    raw[0] = cfg_common_bar - 2.0 * s / 3.0;
    excess = cfg_common_bar + s / 3.0;
    raw[1] = excess - a;
    raw[2] = excess + a;
    s = fmin(fmax(raw[0], cfg_minimum_bar), cfg_maximum_bar);
    target[0] = s;
    offset_data[0] = (fabs(raw[0] - s) > 1.0E-10);
    s = fmin(fmax(raw[1], cfg_minimum_bar), cfg_maximum_bar);
    target[1] = s;
    offset_data[1] = (fabs(raw[1] - s) > 1.0E-10);
    s = fmin(fmax(raw[2], cfg_minimum_bar), cfg_maximum_bar);
    target[2] = s;
    offset_data[2] = (fabs(raw[2] - s) > 1.0E-10);
    if (vectorAny(offset_data, 3) || (fabs(ds - dsWanted) > 1.0E-10) ||
        (fabs(da - daWanted) > 1.0E-10)) {
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
    s = (requested[0] + forceBias[0]) / cfg_forcePerBar_N;
    offset[1] = s;
    excess = s;
    s = (requested[1] + forceBias[1]) / cfg_forcePerBar_N;
    offset[2] = s;
    excess += s;
    s = fmax((cfg_common_bar - cfg_activation_bar) - excess / 3.0,
             -minimum(offset));
    r = _mm_loadu_pd(&offset[0]);
    r = _mm_add_pd(_mm_set1_pd(s), r);
    _mm_storeu_pd(&offset[0], r);
    _mm_storeu_pd(&b_dv[0], r);
    delta[0] = fmax(b_dv[0], 0.0);
    delta[1] = fmax(b_dv[1], 0.0);
    r = _mm_loadu_pd(&delta[0]);
    _mm_storeu_pd(&raw[0], r);
    r1 = _mm_set1_pd(cfg_activation_bar);
    _mm_storeu_pd(&b_candidate[0], _mm_add_pd(r1, r));
    s = fmax(s + offset[2], 0.0);
    b_candidate[2] = cfg_activation_bar + s;
    nz = ((raw[0] <= 1.0E-10) + (raw[1] <= 1.0E-10)) + (s <= 1.0E-10);
    excess = ((b_candidate[0] + b_candidate[1]) + b_candidate[2]) -
             3.0 * cfg_common_bar;
    guard1 = false;
    if (excess > 1.0E-9) {
      if (nz > 0) {
        if (raw[0] <= 1.0E-10) {
          b_candidate[0] -= excess / (double)nz;
        }
        if (raw[1] <= 1.0E-10) {
          b_candidate[1] -= excess / (double)nz;
        }
        if (s <= 1.0E-10) {
          b_candidate[2] -= excess / (double)nz;
        }
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
    if (guard1) {
      boolean_T c_candidate[3];
      c_candidate[0] = (b_candidate[0] >= cfg_minimum_bar - 1.0E-10);
      c_candidate[1] = (b_candidate[1] >= cfg_minimum_bar - 1.0E-10);
      c_candidate[2] = (b_candidate[2] >= cfg_minimum_bar - 1.0E-10);
      if (all(c_candidate)) {
        c_candidate[0] = (b_candidate[0] <= cfg_maximum_bar + 1.0E-10);
        c_candidate[1] = (b_candidate[1] <= cfg_maximum_bar + 1.0E-10);
        c_candidate[2] = (b_candidate[2] <= cfg_maximum_bar + 1.0E-10);
        if (all(c_candidate) &&
            (fabs(((b_candidate[1] + b_candidate[2]) / 2.0 - b_candidate[0]) -
                  cfg_sTrim_bar) <= cfg_maxDeltaS_bar + 1.0E-10) &&
            (fabs((b_candidate[2] - b_candidate[1]) / 2.0 - cfg_aTrim_bar) <=
             cfg_maxDeltaA_bar + 1.0E-10)) {
          target[0] =
              fmin(fmax(b_candidate[0], cfg_minimum_bar), cfg_maximum_bar);
          target[1] =
              fmin(fmax(b_candidate[1], cfg_minimum_bar), cfg_maximum_bar);
          target[2] =
              fmin(fmax(b_candidate[2], cfg_minimum_bar), cfg_maximum_bar);
          limited = false;
        }
      }
    }
    r = _mm_loadu_pd(&target[0]);
    _mm_storeu_pd(&dv1[0], _mm_sub_pd(r, r1));
    delta[0] = fmax(dv1[0], 0.0);
    delta[1] = fmax(dv1[1], 0.0);
    r = _mm_loadu_pd(&delta[0]);
    _mm_storeu_pd(&raw[0], r);
    equivalentForce[0] = cfg_forcePerBar_N * (raw[1] - raw[0]) - forceBias[0];
    equivalentForce[1] =
        cfg_forcePerBar_N *
            (fmax(target[2] - cfg_activation_bar, 0.0) - raw[0]) -
        forceBias[1];
    r = _mm_loadu_pd(&equivalentForce[0]);
    _mm_storeu_pd(&residual[0], _mm_sub_pd(r, _mm_loadu_pd(&requested[0])));
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
