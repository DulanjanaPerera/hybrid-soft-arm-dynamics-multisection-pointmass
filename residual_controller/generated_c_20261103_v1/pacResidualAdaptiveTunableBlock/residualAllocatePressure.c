/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualAllocatePressure.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "residualAllocatePressure.h"
#include "all.h"
#include "allOrAny.h"
#include "mldivide.h"
#include "norm.h"
#include "pacAllocatePressure.h"
#include "pacResidualAdaptiveTunableBlock_data.h"
#include "pacResidualAdaptiveTunableBlock_internal_types.h"
#include "pacResidualAdaptiveTunableBlock_rtwutil.h"
#include "pacResidualAdaptiveTunableBlock_types.h"
#include "residualForceMap.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Declarations */
static void candidate(const double u[2], const double request[2],
                      const double anchor[2], const double bias[2],
                      const double previous[3], double step,
                      const double measured[6], double dt, const double h[128],
                      const f_struct_T *network, const e_struct_T *cfg,
                      const double shift[2], double p[3], double e[2],
                      double target[3]);

/* Function Definitions */
/*
 * Arguments    : const double u[2]
 *                const double request[2]
 *                const double anchor[2]
 *                const double bias[2]
 *                const double previous[3]
 *                double step
 *                const double measured[6]
 *                double dt
 *                const double h[128]
 *                const f_struct_T *network
 *                const e_struct_T *cfg
 *                const double shift[2]
 *                double p[3]
 *                double e[2]
 *                double target[3]
 * Return Type  : void
 */
static void candidate(const double u[2], const double request[2],
                      const double anchor[2], const double bias[2],
                      const double previous[3], double step,
                      const double measured[6], double dt, const double h[128],
                      const f_struct_T *network, const e_struct_T *cfg,
                      const double shift[2], double p[3], double e[2],
                      double target[3])
{
  __m128d r;
  __m128d r1;
  double b_u[384];
  double w[384];
  double z[8];
  double a[2];
  double a__7[2];
  double b_dv[2];
  double dv1[2];
  double dv2[2];
  double d;
  double d1;
  double d2;
  double d3;
  int i;
  int i1;
  int i2;
  int k;
  pacAllocatePressure(u, anchor, bias, cfg->common_bar, cfg->sTrim_bar,
                      cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar,
                      cfg->maxDeltaS_bar, cfg->maxDeltaA_bar,
                      cfg->activation_bar, cfg->forcePerBar_N,
                      cfg->enableCommonAdaptation, target, a, a__7);
  r = _mm_loadu_pd(&previous[0]);
  r1 = _mm_set1_pd(step);
  _mm_storeu_pd(&b_dv[0], _mm_sub_pd(r, r1));
  _mm_storeu_pd(&dv1[0], _mm_add_pd(r, r1));
  dv2[0] = fmin(fmax(target[0], b_dv[0]), dv1[0]);
  dv2[1] = fmin(fmax(target[1], b_dv[1]), dv1[1]);
  r = _mm_loadu_pd(&dv2[0]);
  _mm_storeu_pd(&p[0], r);
  p[2] = fmin(fmax(target[2], previous[2] - step), previous[2] + step);
  /* RESIDUALFORCEMAP Activation re-expression keeps F_cfg+r_eff=B0*p+r0+deltaG.
   */
  /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
   * evaluation. */
  /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history fixed. */
  _mm_storeu_pd(&z[0], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                             _mm_loadu_pd(&network->mean[0])),
                                  _mm_loadu_pd(&network->std[0])));
  _mm_storeu_pd(&z[2], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                             _mm_loadu_pd(&network->mean[2])),
                                  _mm_loadu_pd(&network->std[2])));
  r = _mm_loadu_pd(&p[0]);
  _mm_storeu_pd(&z[4],
                _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&network->mean[4])),
                           _mm_loadu_pd(&network->std[4])));
  z[6] = (p[2] - network->mean[6]) / network->std[6];
  z[7] = (dt - network->mean[7]) / network->std[7];
#pragma omp parallel for num_threads(omp_get_max_threads()) private(d, i2)

  for (i = 0; i < 384; i++) {
    d = 0.0;
    for (i2 = 0; i2 < 8; i2++) {
      d += network->input[i + 384 * i2] * z[i2];
    }
    b_u[i] = d + network->bias[i];
  }
  memset(&w[0], 0, 384U * sizeof(double));
  for (k = 0; k < 128; k++) {
    for (i1 = 0; i1 <= 382; i1 += 2) {
      r = _mm_loadu_pd(&w[i1]);
      _mm_storeu_pd(
          &w[i1],
          _mm_add_pd(r,
                     _mm_mul_pd(_mm_loadu_pd(&network->recurrent[i1 + 384 * k]),
                                _mm_set1_pd(h[k]))));
    }
  }
  memset(&a[0], 0, sizeof(double) << 1);
  for (k = 0; k < 128; k++) {
    d1 = 0.5 * (tanh(0.5 * (b_u[k + 128] + w[k + 128])) + 1.0);
    d1 = d1 * h[k] +
         (1.0 - d1) *
             tanh(b_u[k + 256] +
                  0.5 * (tanh(0.5 * (b_u[k] + w[k])) + 1.0) * w[k + 256]);
    r = _mm_loadu_pd(&a[0]);
    _mm_storeu_pd(
        &a[0], _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->output[k << 1]),
                                        _mm_set1_pd(d1))));
  }
  memset(&a__7[0], 0, sizeof(double) << 1);
  for (k = 0; k < 8; k++) {
    r = _mm_loadu_pd(&a__7[0]);
    _mm_storeu_pd(&a__7[0],
                  _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->skip[k << 1]),
                                           _mm_set1_pd(z[k]))));
  }
  d1 = p[0];
  d2 = p[1];
  d3 = p[2];
  for (k = 0; k < 2; k++) {
    double d4;
    d4 = tanh(10.0 * ((a[k] + a__7[k]) + network->outputBias[k]) / 80.0);
    a[k] = d4;
    e[k] = (((((-13.273228961416876 * d1 + dv[k + 2] * d2) + dv[k + 4] * d3) +
              80.0 * d4) +
             shift[k]) -
            bias[k]) -
           request[k];
  }
}

/*
 * RESIDUALALLOCATEPRESSURE Bounded nonlinear inverse around existing allocator.
 *  Existing fixed/variable-common allocation and its rails/s/a caps are reused.
 *  Slew is applied inside EACH candidate evaluation, not after solving another
 * p. History is immutable here. No Optimization Toolbox or hardware I/O.
 *
 * Arguments    : const double request[2]
 *                const double anchor[2]
 *                const double anchorResidual[2]
 *                const double bias[2]
 *                const double previous[3]
 *                const double lastRequest[2]
 *                const double lastEffective[2]
 *                double step
 *                const double measured[6]
 *                double featureDt
 *                const double h[128]
 *                const f_struct_T *network
 *                const e_struct_T *cfg
 *                const double shift[2]
 *                const double options[8]
 *                double pressure[3]
 *                double b_error[2]
 *                unsigned char *iterations
 *                double linearRequest[2]
 *                double target[3]
 * Return Type  : boolean_T
 */
boolean_T residualAllocatePressure(
    const double request[2], const double anchor[2],
    const double anchorResidual[2], const double bias[2],
    const double previous[3], const double lastRequest[2],
    const double lastEffective[2], double step, const double measured[6],
    double featureDt, const double h[128], const f_struct_T *network,
    const e_struct_T *cfg, const double shift[2], const double options[8],
    double pressure[3], double b_error[2], unsigned char *iterations,
    double linearRequest[2], double target[3])
{
  static const double c_y[4] = {1.0E-6, 0.0, 0.0, 1.0E-6};
  __m128d r;
  __m128d r1;
  __m128d r2;
  __m128d r3;
  double W[1152];
  double w[384];
  double c[128];
  double g[128];
  double z[8];
  double b_network[6];
  double b_y[6];
  double seeds[6];
  double J[4];
  double minus[3];
  double plus[3];
  double tar[3];
  double varargin_2[3];
  double a[2];
  double b_dv[2];
  double b_linearRequest[2];
  double dv1[2];
  double dv2[2];
  double dv3[2];
  double e[2];
  double move[2];
  double y[2];
  double best;
  double cost;
  int i;
  int it;
  int k;
  int y_tmp;
  boolean_T exitg1;
  boolean_T limited;
  residualForceMap(previous, measured, featureDt, h, network->input,
                   network->recurrent, network->bias, network->output,
                   network->skip, network->outputBias, network->mean,
                   network->std, cfg->activation_bar, cfg->forcePerBar_N, shift,
                   a, move);
  r = _mm_loadu_pd(&request[0]);
  r1 = _mm_loadu_pd(&anchorResidual[0]);
  r2 = _mm_sub_pd(r, r1);
  _mm_storeu_pd(&seeds[0], r2);
  r3 = _mm_loadu_pd(&move[0]);
  _mm_storeu_pd(&seeds[2], _mm_sub_pd(r, r3));
  _mm_storeu_pd(&seeds[4],
                _mm_sub_pd(_mm_add_pd(_mm_loadu_pd(&lastRequest[0]), r),
                           _mm_loadu_pd(&lastEffective[0])));
  _mm_storeu_pd(&linearRequest[0], r2);
  _mm_storeu_pd(&b_dv[0], _mm_sub_pd(_mm_loadu_pd(&anchor[0]), r1));
  candidate(&seeds[0], request, b_dv, bias, previous, step, measured, featureDt,
            h, network, cfg, shift, pressure, b_error, target);
  r = _mm_loadu_pd(&b_error[0]);
  _mm_storeu_pd(&y[0], _mm_mul_pd(r, r));
  best = y[0] + y[1];
  *iterations = 0U;
  for (k = 0; k < 2; k++) {
    y_tmp = (k + 1) << 1;
    candidate(&seeds[y_tmp], request, b_dv, bias, previous, step, measured,
              featureDt, h, network, cfg, shift, minus, e, tar);
    r = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&y[0], _mm_mul_pd(r, r));
    cost = y[0] + y[1];
    if (cost < best) {
      best = cost;
      pressure[0] = minus[0];
      pressure[1] = minus[1];
      pressure[2] = minus[2];
      b_error[0] = e[0];
      b_error[1] = e[1];
      target[0] = tar[0];
      target[1] = tar[1];
      target[2] = tar[2];
      linearRequest[0] = seeds[y_tmp];
      linearRequest[1] = seeds[y_tmp + 1];
    }
  }
  it = 0;
  exitg1 = false;
  while ((!exitg1) && ((it < 12) && (!((double)it + 1.0 > options[1])) &&
                       (!(sqrt(best) <= options[2])))) {
    double u[384];
    double b_g[128];
    double b_r[128];
    double c_tmp[128];
    double next_tmp[128];
    double b_y_tmp[4];
    double c_tmp_tmp;
    double d;
    double d1;
    double d2;
    double d3;
    boolean_T bv[2];
    *iterations = (unsigned char)(it + 1);
    /* RESIDUALFORCEMAP Activation re-expression keeps
     * F_cfg+r_eff=B0*p+r0+deltaG. */
    /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
     * evaluation. */
    /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history fixed.
     */
    _mm_storeu_pd(&z[0], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                               _mm_loadu_pd(&network->mean[0])),
                                    _mm_loadu_pd(&network->std[0])));
    _mm_storeu_pd(&z[2], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                               _mm_loadu_pd(&network->mean[2])),
                                    _mm_loadu_pd(&network->std[2])));
    r = _mm_loadu_pd(&pressure[0]);
    _mm_storeu_pd(&z[4],
                  _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&network->mean[4])),
                             _mm_loadu_pd(&network->std[4])));
    z[6] = (pressure[2] - network->mean[6]) / network->std[6];
    z[7] = (featureDt - network->mean[7]) / network->std[7];
    for (k = 0; k < 384; k++) {
      cost = 0.0;
      for (i = 0; i < 8; i++) {
        cost += network->input[k + 384 * i] * z[i];
      }
      u[k] = cost + network->bias[k];
    }
    memset(&w[0], 0, 384U * sizeof(double));
    for (k = 0; k < 128; k++) {
      for (i = 0; i <= 382; i += 2) {
        r = _mm_loadu_pd(&w[i]);
        _mm_storeu_pd(
            &w[i],
            _mm_add_pd(
                r, _mm_mul_pd(_mm_loadu_pd(&network->recurrent[i + 384 * k]),
                              _mm_set1_pd(h[k]))));
      }
    }
    memset(&a[0], 0, sizeof(double) << 1);
    for (k = 0; k < 128; k++) {
      cost = 0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0);
      b_r[k] = cost;
      d = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
      b_g[k] = d;
      c_tmp_tmp = w[k + 256];
      c_tmp[k] = c_tmp_tmp;
      cost = tanh(u[k + 256] + cost * c_tmp_tmp);
      c[k] = cost;
      next_tmp[k] = 1.0 - d;
      r = _mm_loadu_pd(&a[0]);
      _mm_storeu_pd(
          &a[0],
          _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->output[k << 1]),
                                   _mm_set1_pd(d * h[k] + (1.0 - d) * cost))));
    }
    memset(&y[0], 0, sizeof(double) << 1);
    for (k = 0; k < 8; k++) {
      r = _mm_loadu_pd(&y[0]);
      _mm_storeu_pd(
          &y[0], _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->skip[k << 1]),
                                          _mm_set1_pd(z[k]))));
    }
    a[0] = tanh(10.0 * ((a[0] + y[0]) + network->outputBias[0]) / 80.0);
    a[1] = tanh(10.0 * ((a[1] + y[1]) + network->outputBias[1]) / 80.0);
    for (k = 0; k < 3; k++) {
      for (i = 0; i <= 382; i += 2) {
        _mm_storeu_pd(
            &W[i + 384 * k],
            _mm_div_pd(_mm_loadu_pd(&network->input[i + 384 * (k + 4)]),
                       _mm_set1_pd(network->std[k + 4])));
      }
    }
    for (k = 0; k <= 126; k += 2) {
      r = _mm_loadu_pd(&c[k]);
      _mm_storeu_pd(&g[k], _mm_mul_pd(r, r));
    }
    r = _mm_loadu_pd(&a[0]);
    _mm_storeu_pd(&y[0],
                  _mm_mul_pd(_mm_set1_pd(10.0),
                             _mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(r, r))));
    cost = 2.0 * options[5];
    e[1] = 0.0;
    e[0] = options[5];
    r = _mm_loadu_pd(&linearRequest[0]);
    r1 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_add_pd(r, r1));
    pacAllocatePressure(
        b_linearRequest, b_dv, bias, cfg->common_bar, cfg->sTrim_bar,
        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar, cfg->maxDeltaS_bar,
        cfg->maxDeltaA_bar, cfg->activation_bar, cfg->forcePerBar_N,
        cfg->enableCommonAdaptation, tar, move, a);
    r = _mm_loadu_pd(&previous[0]);
    r2 = _mm_set1_pd(step);
    r1 = _mm_sub_pd(r, r2);
    _mm_storeu_pd(&minus[0], r1);
    _mm_storeu_pd(&dv1[0], r1);
    r = _mm_add_pd(r, r2);
    _mm_storeu_pd(&varargin_2[0], r);
    _mm_storeu_pd(&dv2[0], r);
    dv3[0] = fmin(fmax(tar[0], dv1[0]), dv2[0]);
    dv3[1] = fmin(fmax(tar[1], dv1[1]), dv2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&plus[0], r);
    c_tmp_tmp = previous[2] - step;
    minus[2] = c_tmp_tmp;
    d = previous[2] + step;
    plus[2] = fmin(fmax(tar[2], c_tmp_tmp), d);
    r = _mm_loadu_pd(&linearRequest[0]);
    r1 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_sub_pd(r, r1));
    pacAllocatePressure(
        b_linearRequest, b_dv, bias, cfg->common_bar, cfg->sTrim_bar,
        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar, cfg->maxDeltaS_bar,
        cfg->maxDeltaA_bar, cfg->activation_bar, cfg->forcePerBar_N,
        cfg->enableCommonAdaptation, tar, move, a);
    dv3[0] = fmin(fmax(tar[0], minus[0]), varargin_2[0]);
    dv3[1] = fmin(fmax(tar[1], minus[1]), varargin_2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&minus[0], r);
    r1 = _mm_loadu_pd(&plus[0]);
    r3 = _mm_set1_pd(cost);
    _mm_storeu_pd(&seeds[0], _mm_div_pd(_mm_sub_pd(r1, r), r3));
    seeds[2] = (plus[2] - fmin(fmax(tar[2], minus[2]), d)) / cost;
    e[0] = 0.0;
    e[1] = options[5];
    r = _mm_loadu_pd(&linearRequest[0]);
    r1 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_add_pd(r, r1));
    pacAllocatePressure(
        b_linearRequest, b_dv, bias, cfg->common_bar, cfg->sTrim_bar,
        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar, cfg->maxDeltaS_bar,
        cfg->maxDeltaA_bar, cfg->activation_bar, cfg->forcePerBar_N,
        cfg->enableCommonAdaptation, tar, move, a);
    r = _mm_loadu_pd(&previous[0]);
    r1 = _mm_sub_pd(r, r2);
    _mm_storeu_pd(&minus[0], r1);
    _mm_storeu_pd(&dv1[0], r1);
    r = _mm_add_pd(r, r2);
    _mm_storeu_pd(&varargin_2[0], r);
    _mm_storeu_pd(&dv2[0], r);
    dv3[0] = fmin(fmax(tar[0], dv1[0]), dv2[0]);
    dv3[1] = fmin(fmax(tar[1], dv1[1]), dv2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&plus[0], r);
    minus[2] = c_tmp_tmp;
    plus[2] = fmin(fmax(tar[2], c_tmp_tmp), d);
    r = _mm_loadu_pd(&linearRequest[0]);
    r1 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_sub_pd(r, r1));
    pacAllocatePressure(
        b_linearRequest, b_dv, bias, cfg->common_bar, cfg->sTrim_bar,
        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar, cfg->maxDeltaS_bar,
        cfg->maxDeltaA_bar, cfg->activation_bar, cfg->forcePerBar_N,
        cfg->enableCommonAdaptation, tar, move, a);
    dv3[0] = fmin(fmax(tar[0], minus[0]), varargin_2[0]);
    dv3[1] = fmin(fmax(tar[1], minus[1]), varargin_2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&minus[0], r);
    r1 = _mm_loadu_pd(&plus[0]);
    _mm_storeu_pd(&seeds[3], _mm_div_pd(_mm_sub_pd(r1, r), r3));
    seeds[5] = (plus[2] - fmin(fmax(tar[2], minus[2]), d)) / cost;
    memset(&b_network[0], 0, 6U * sizeof(double));
    r1 = _mm_loadu_pd(&y[0]);
    for (i = 0; i < 3; i++) {
      y_tmp = i << 1;
      for (k = 0; k < 128; k++) {
        int i1;
        r = _mm_loadu_pd(&b_network[y_tmp]);
        i1 = k + 384 * i;
        cost = next_tmp[k];
        c_tmp_tmp = b_r[k];
        _mm_storeu_pd(
            &b_network[y_tmp],
            _mm_add_pd(
                r, _mm_mul_pd(
                       _mm_loadu_pd(&network->output[k << 1]),
                       _mm_set1_pd(
                           cost * ((1.0 - g[k]) *
                                   (W[i1 + 256] +
                                    c_tmp[k] * (c_tmp_tmp * (1.0 - c_tmp_tmp) *
                                                W[i1]))) +
                           (h[k] - c[k]) * (b_g[k] * cost * W[i1 + 128])))));
      }
      y_tmp = i << 1;
      r = _mm_loadu_pd(&b_network[y_tmp]);
      _mm_storeu_pd(
          &b_y[y_tmp],
          _mm_add_pd(
              _mm_mul_pd(
                  r1,
                  _mm_add_pd(
                      r, _mm_div_pd(_mm_loadu_pd(&network->skip[(i + 4) << 1]),
                                    _mm_set1_pd(network->std[i + 4])))),
              _mm_loadu_pd(&dv[y_tmp])));
    }
    memset(&J[0], 0, sizeof(double) << 2);
    for (k = 0; k < 2; k++) {
      r = _mm_loadu_pd(&b_y[0]);
      y_tmp = k << 1;
      r1 = _mm_loadu_pd(&J[y_tmp]);
      _mm_storeu_pd(&J[y_tmp],
                    _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(seeds[3 * k]))));
      r = _mm_loadu_pd(&b_y[2]);
      r1 = _mm_loadu_pd(&J[y_tmp]);
      _mm_storeu_pd(
          &J[y_tmp],
          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(seeds[3 * k + 1]))));
      r = _mm_loadu_pd(&b_y[4]);
      r1 = _mm_loadu_pd(&J[y_tmp]);
      _mm_storeu_pd(
          &J[y_tmp],
          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(seeds[3 * k + 2]))));
    }
    b_y_tmp[0] = J[0];
    b_y_tmp[1] = J[2];
    b_y_tmp[2] = J[1];
    b_y_tmp[3] = J[3];
    memset(&y[0], 0, sizeof(double) << 1);
    cost = J[0];
    c_tmp_tmp = J[1];
    d = J[2];
    d1 = J[3];
    d2 = y[0];
    d3 = y[1];
    for (k = 0; k < 2; k++) {
      double d4;
      double d5;
      d4 = b_y_tmp[k + 2];
      J[k] = -((b_y_tmp[k] * cost + d4 * c_tmp_tmp) + c_y[k]);
      y_tmp = k << 1;
      d5 = b_error[k];
      d2 += b_y_tmp[y_tmp] * d5;
      J[k + 2] = -((b_y_tmp[k] * d + d4 * d1) + c_y[k + 2]);
      d3 += b_y_tmp[y_tmp + 1] * d5;
    }
    y[1] = d3;
    y[0] = d2;
    mldivide(J, y, move);
    cost = fmin(1.0, options[6] / fmax(b_norm(move), 1.0E-12));
    c_tmp_tmp = move[0] * cost;
    move[0] = c_tmp_tmp;
    bv[0] = ((!rtIsInf(c_tmp_tmp)) && (!rtIsNaN(c_tmp_tmp)));
    c_tmp_tmp = move[1] * cost;
    move[1] = c_tmp_tmp;
    bv[1] = ((!rtIsInf(c_tmp_tmp)) && (!rtIsNaN(c_tmp_tmp)));
    if ((!c_all(bv)) || (b_norm(move) < 1.0E-10)) {
      exitg1 = true;
    } else {
      boolean_T exitg2;
      limited = false;
      y_tmp = 0;
      exitg2 = false;
      while ((!exitg2) && (y_tmp < 4)) {
        r = _mm_loadu_pd(&move[0]);
        r1 = _mm_loadu_pd(&linearRequest[0]);
        _mm_storeu_pd(
            &a[0],
            _mm_add_pd(r1,
                       _mm_mul_pd(_mm_set1_pd(rt_powd_snf(0.5, y_tmp)), r)));
        candidate(a, request, b_dv, bias, previous, step, measured, featureDt,
                  h, network, cfg, shift, minus, e, tar);
        r = _mm_loadu_pd(&e[0]);
        _mm_storeu_pd(&y[0], _mm_mul_pd(r, r));
        cost = y[0] + y[1];
        if ((!rtIsInf(cost)) && (!rtIsNaN(cost)) && (cost < best - 1.0E-14)) {
          best = cost;
          pressure[0] = minus[0];
          pressure[1] = minus[1];
          pressure[2] = minus[2];
          b_error[0] = e[0];
          b_error[1] = e[1];
          target[0] = tar[0];
          target[1] = tar[1];
          target[2] = tar[2];
          linearRequest[0] = a[0];
          linearRequest[1] = a[1];
          limited = true;
          exitg2 = true;
        } else {
          y_tmp++;
        }
      }
      if (!limited) {
        exitg1 = true;
      } else {
        it++;
      }
    }
  }
  if (b_norm(b_error) > options[2]) {
    limited = true;
  } else {
    boolean_T tar_data[10];
    tar_data[0] = (fabs(pressure[0] - target[0]) > 1.0E-9);
    tar_data[1] = (fabs(pressure[1] - target[1]) > 1.0E-9);
    tar_data[2] = (fabs(pressure[2] - target[2]) > 1.0E-9);
    if (vectorAny(tar_data, 3)) {
      limited = true;
    } else {
      limited = false;
    }
  }
  return limited;
}

/*
 * File trailer for residualAllocatePressure.c
 *
 * [EOF]
 */
