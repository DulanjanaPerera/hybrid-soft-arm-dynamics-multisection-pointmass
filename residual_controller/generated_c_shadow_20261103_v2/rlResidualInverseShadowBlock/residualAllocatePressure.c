/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualAllocatePressure.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:25
 */

/* Include Files */
#include "residualAllocatePressure.h"
#include "all.h"
#include "any.h"
#include "mldivide.h"
#include "norm.h"
#include "pacAllocatePressure.h"
#include "rlResidualInverseShadowBlock_data.h"
#include "rlResidualInverseShadowBlock_internal_types.h"
#include "rlResidualInverseShadowBlock_rtwutil.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Declarations */
static void
candidate(const double u[2], const double request[2], const double previous[3],
          double step, const double measured[4], double dt, const double h[128],
          const double network_input[3072],
          const double network_recurrent[49152], const double network_bias[384],
          const double network_output[256], const double network_skip[16],
          const double network_outputBias[2], const double network_mean[8],
          const double network_std[8], const b_struct_T *cfg, double p[3],
          double e[2], double target[3]);

/* Function Definitions */
/*
 * Arguments    : const double u[2]
 *                const double request[2]
 *                const double previous[3]
 *                double step
 *                const double measured[4]
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
 *                const b_struct_T *cfg
 *                double p[3]
 *                double e[2]
 *                double target[3]
 * Return Type  : void
 */
static void
candidate(const double u[2], const double request[2], const double previous[3],
          double step, const double measured[4], double dt, const double h[128],
          const double network_input[3072],
          const double network_recurrent[49152], const double network_bias[384],
          const double network_output[256], const double network_skip[16],
          const double network_outputBias[2], const double network_mean[8],
          const double network_std[8], const b_struct_T *cfg, double p[3],
          double e[2], double target[3])
{
  __m128d r;
  __m128d r1;
  double b_u[384];
  double w[384];
  double z[8];
  double a[2];
  double a__7[2];
  double b_dv[2];
  double b_network_skip[2];
  double dv1[2];
  double d;
  double d1;
  double d2;
  int i;
  int k;
  b_network_skip[0] = 0.0;
  b_network_skip[1] = 0.0;
  pacAllocatePressure(u, b_network_skip, cfg->common_bar, cfg->sTrim_bar,
                      cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar,
                      cfg->maxDeltaS_bar, cfg->maxDeltaA_bar,
                      cfg->activation_bar, cfg->forcePerBar_N,
                      cfg->enableCommonAdaptation, target, a, a__7);
  r = _mm_loadu_pd(&previous[0]);
  r1 = _mm_set1_pd(step);
  _mm_storeu_pd(&b_dv[0], _mm_sub_pd(r, r1));
  _mm_storeu_pd(&dv1[0], _mm_add_pd(r, r1));
  b_network_skip[0] = fmin(fmax(target[0], b_dv[0]), dv1[0]);
  b_network_skip[1] = fmin(fmax(target[1], b_dv[1]), dv1[1]);
  r = _mm_loadu_pd(&b_network_skip[0]);
  _mm_storeu_pd(&p[0], r);
  p[2] = fmin(fmax(target[2], previous[2] - step), previous[2] + step);
  /* RESIDUALFORCEMAP Activation re-expression keeps F_cfg+r_eff=B0*p+r0+deltaG.
   */
  /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
   * evaluation. */
  /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history fixed. */
  _mm_storeu_pd(&z[0], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                             _mm_loadu_pd(&network_mean[0])),
                                  _mm_loadu_pd(&network_std[0])));
  _mm_storeu_pd(&z[2], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                             _mm_loadu_pd(&network_mean[2])),
                                  _mm_loadu_pd(&network_std[2])));
  r = _mm_loadu_pd(&p[0]);
  _mm_storeu_pd(&z[4], _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&network_mean[4])),
                                  _mm_loadu_pd(&network_std[4])));
  z[6] = (p[2] - network_mean[6]) / network_std[6];
  z[7] = (dt - network_mean[7]) / network_std[7];
  for (k = 0; k < 384; k++) {
    d = 0.0;
    for (i = 0; i < 8; i++) {
      d += network_input[k + 384 * i] * z[i];
    }
    b_u[k] = d + network_bias[k];
  }
  memset(&w[0], 0, 384U * sizeof(double));
  for (k = 0; k < 128; k++) {
    for (i = 0; i <= 382; i += 2) {
      r = _mm_loadu_pd(&w[i]);
      _mm_storeu_pd(
          &w[i],
          _mm_add_pd(r,
                     _mm_mul_pd(_mm_loadu_pd(&network_recurrent[i + 384 * k]),
                                _mm_set1_pd(h[k]))));
    }
  }
  memset(&a[0], 0, sizeof(double) << 1);
  for (k = 0; k < 128; k++) {
    d = 0.5 * (tanh(0.5 * (b_u[k + 128] + w[k + 128])) + 1.0);
    d = d * h[k] + (1.0 - d) * tanh(b_u[k + 256] +
                                    0.5 * (tanh(0.5 * (b_u[k] + w[k])) + 1.0) *
                                        w[k + 256]);
    r = _mm_loadu_pd(&a[0]);
    _mm_storeu_pd(
        &a[0], _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network_output[k << 1]),
                                        _mm_set1_pd(d))));
  }
  memset(&b_network_skip[0], 0, sizeof(double) << 1);
  for (k = 0; k < 8; k++) {
    r = _mm_loadu_pd(&b_network_skip[0]);
    _mm_storeu_pd(&b_network_skip[0],
                  _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network_skip[k << 1]),
                                           _mm_set1_pd(z[k]))));
  }
  d = p[0];
  d1 = p[1];
  d2 = p[2];
  for (k = 0; k < 2; k++) {
    double d3;
    d3 = tanh(10.0 * ((a[k] + b_network_skip[k]) + network_outputBias[k]) /
              80.0);
    a[k] = d3;
    e[k] = (((-13.273228961416876 * d + dv[k + 2] * d1) + dv[k + 4] * d2) +
            80.0 * d3) -
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
 *                const double previous[3]
 *                const double lastRequest[2]
 *                const double lastEffective[2]
 *                double step
 *                const double measured[4]
 *                double featureDt
 *                const double h[128]
 *                const double network_input[3072]
 *                const double network_recurrent[49152]
 *                const double network_bias[384]
 *                const double network_output[256]
 *                const double network_skip[16]
 *                const double network_outputBias[2]
 *                const double network_mean[8]
 *                const double network_std[8]
 *                const b_struct_T *cfg
 *                const double options[8]
 *                double pressure[3]
 *                double b_error[2]
 *                unsigned char *iterations
 * Return Type  : boolean_T
 */
boolean_T residualAllocatePressure(
    const double request[2], const double previous[3],
    const double lastRequest[2], const double lastEffective[2], double step,
    const double measured[4], double featureDt, const double h[128],
    const double network_input[3072], const double network_recurrent[49152],
    const double network_bias[384], const double network_output[256],
    const double network_skip[16], const double network_outputBias[2],
    const double network_mean[8], const double network_std[8],
    const b_struct_T *cfg, const double options[8], double pressure[3],
    double b_error[2], unsigned char *iterations)
{
  static const double b_y[4] = {1.0E-6, 0.0, 0.0, 1.0E-6};
  __m128d r;
  __m128d r1;
  __m128d r2;
  double W[1152];
  double u[384];
  double w[384];
  double b_g[128];
  double c[128];
  double z[8];
  double b_network_output[6];
  double seeds[6];
  double y[6];
  double J[4];
  double ep[3];
  double minus[3];
  double plus[3];
  double target[3];
  double varargin_2[3];
  double B[2];
  double a[2];
  double b_dv[2];
  double b_linearRequest[2];
  double dv1[2];
  double dv2[2];
  double dv3[2];
  double e[2];
  double linearRequest[2];
  double trial[2];
  double best;
  double c_tmp_tmp;
  double cost;
  double d;
  double z_tmp;
  int B_tmp;
  int i;
  int it;
  int k;
  boolean_T exitg1;
  boolean_T limited;
  /* RESIDUALFORCEMAP Activation re-expression keeps F_cfg+r_eff=B0*p+r0+deltaG.
   */
  /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
   * evaluation. */
  /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history fixed. */
  _mm_storeu_pd(&z[0], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                             _mm_loadu_pd(&network_mean[0])),
                                  _mm_loadu_pd(&network_std[0])));
  _mm_storeu_pd(&z[2], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                             _mm_loadu_pd(&network_mean[2])),
                                  _mm_loadu_pd(&network_std[2])));
  _mm_storeu_pd(&z[4], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&previous[0]),
                                             _mm_loadu_pd(&network_mean[4])),
                                  _mm_loadu_pd(&network_std[4])));
  z[6] = (previous[2] - network_mean[6]) / network_std[6];
  z_tmp = (featureDt - network_mean[7]) / network_std[7];
  z[7] = z_tmp;
  for (k = 0; k < 384; k++) {
    cost = 0.0;
    for (i = 0; i < 8; i++) {
      cost += network_input[k + 384 * i] * z[i];
    }
    u[k] = cost + network_bias[k];
  }
  memset(&w[0], 0, 384U * sizeof(double));
  for (k = 0; k < 128; k++) {
    for (i = 0; i <= 382; i += 2) {
      r = _mm_loadu_pd(&w[i]);
      _mm_storeu_pd(
          &w[i],
          _mm_add_pd(r,
                     _mm_mul_pd(_mm_loadu_pd(&network_recurrent[i + 384 * k]),
                                _mm_set1_pd(h[k]))));
    }
  }
  memset(&a[0], 0, sizeof(double) << 1);
  for (k = 0; k < 128; k++) {
    cost = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
    r = _mm_loadu_pd(&a[0]);
    _mm_storeu_pd(
        &a[0],
        _mm_add_pd(
            r, _mm_mul_pd(_mm_loadu_pd(&network_output[k << 1]),
                          _mm_set1_pd(
                              cost * h[k] +
                              (1.0 - cost) *
                                  tanh(u[k + 256] +
                                       0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0) *
                                           w[k + 256])))));
  }
  memset(&e[0], 0, sizeof(double) << 1);
  for (k = 0; k < 8; k++) {
    r = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&e[0],
                  _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network_skip[k << 1]),
                                           _mm_set1_pd(z[k]))));
  }
  a[0] = tanh(10.0 * ((a[0] + e[0]) + network_outputBias[0]) / 80.0);
  a[1] = tanh(10.0 * ((a[1] + e[1]) + network_outputBias[1]) / 80.0);
  _mm_storeu_pd(
      &linearRequest[0],
      _mm_sub_pd(_mm_loadu_pd(&previous[0]), _mm_set1_pd(cfg->activation_bar)));
  B[0] = fmax(linearRequest[0], 0.0);
  B[1] = fmax(linearRequest[1], 0.0);
  r = _mm_loadu_pd(&B[0]);
  _mm_storeu_pd(&ep[0], r);
  e[0] = cfg->forcePerBar_N * (ep[1] - ep[0]);
  e[1] = cfg->forcePerBar_N *
         (fmax(previous[2] - cfg->activation_bar, 0.0) - ep[0]);
  cost = previous[0];
  c_tmp_tmp = previous[1];
  d = previous[2];
  r1 = _mm_loadu_pd(&request[0]);
  _mm_storeu_pd(&seeds[0], r1);
  r = _mm_loadu_pd(&dv[2]);
  r = _mm_mul_pd(r, _mm_set1_pd(c_tmp_tmp));
  r2 = _mm_add_pd(_mm_set1_pd(-13.273228961416876 * cost), r);
  r = _mm_loadu_pd(&dv[4]);
  r = _mm_mul_pd(r, _mm_set1_pd(d));
  r2 = _mm_add_pd(r2, r);
  r = _mm_loadu_pd(&a[0]);
  r = _mm_mul_pd(_mm_set1_pd(80.0), r);
  r = _mm_add_pd(r2, r);
  r2 = _mm_loadu_pd(&e[0]);
  r = _mm_sub_pd(r, r2);
  r = _mm_sub_pd(r1, r);
  _mm_storeu_pd(&seeds[2], r);
  r = _mm_loadu_pd(&lastRequest[0]);
  r = _mm_add_pd(r, r1);
  r2 = _mm_loadu_pd(&lastEffective[0]);
  r = _mm_sub_pd(r, r2);
  _mm_storeu_pd(&seeds[4], r);
  _mm_storeu_pd(&linearRequest[0], r1);
  candidate(&seeds[0], request, previous, step, measured, featureDt, h,
            network_input, network_recurrent, network_bias, network_output,
            network_skip, network_outputBias, network_mean, network_std, cfg,
            pressure, b_error, target);
  r = _mm_loadu_pd(&b_error[0]);
  _mm_storeu_pd(&B[0], _mm_mul_pd(r, r));
  best = B[0] + B[1];
  *iterations = 0U;
  for (k = 0; k < 2; k++) {
    B_tmp = (k + 1) << 1;
    candidate(&seeds[B_tmp], request, previous, step, measured, featureDt, h,
              network_input, network_recurrent, network_bias, network_output,
              network_skip, network_outputBias, network_mean, network_std, cfg,
              minus, e, ep);
    r = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&B[0], _mm_mul_pd(r, r));
    cost = B[0] + B[1];
    if (cost < best) {
      best = cost;
      pressure[0] = minus[0];
      pressure[1] = minus[1];
      pressure[2] = minus[2];
      b_error[0] = e[0];
      b_error[1] = e[1];
      target[0] = ep[0];
      target[1] = ep[1];
      target[2] = ep[2];
      linearRequest[0] = seeds[B_tmp];
      linearRequest[1] = seeds[B_tmp + 1];
    }
  }
  it = 0;
  exitg1 = false;
  while ((!exitg1) && ((it < 12) && (!((double)it + 1.0 > options[1])) &&
                       (!(sqrt(best) <= options[2])))) {
    double b_r[128];
    double c_tmp[128];
    double g[128];
    double next_tmp[128];
    double y_tmp[4];
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
                                               _mm_loadu_pd(&network_mean[0])),
                                    _mm_loadu_pd(&network_std[0])));
    _mm_storeu_pd(&z[2], _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                               _mm_loadu_pd(&network_mean[2])),
                                    _mm_loadu_pd(&network_std[2])));
    r = _mm_loadu_pd(&pressure[0]);
    _mm_storeu_pd(&z[4],
                  _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&network_mean[4])),
                             _mm_loadu_pd(&network_std[4])));
    z[6] = (pressure[2] - network_mean[6]) / network_std[6];
    z[7] = z_tmp;
    for (k = 0; k < 384; k++) {
      cost = 0.0;
      for (i = 0; i < 8; i++) {
        cost += network_input[k + 384 * i] * z[i];
      }
      u[k] = cost + network_bias[k];
    }
    memset(&a[0], 0, sizeof(double) << 1);
    for (k = 0; k < 128; k++) {
      cost = 0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0);
      b_r[k] = cost;
      d = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
      g[k] = d;
      c_tmp_tmp = w[k + 256];
      c_tmp[k] = c_tmp_tmp;
      cost = tanh(u[k + 256] + cost * c_tmp_tmp);
      c[k] = cost;
      next_tmp[k] = 1.0 - d;
      r = _mm_loadu_pd(&a[0]);
      _mm_storeu_pd(
          &a[0],
          _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network_output[k << 1]),
                                   _mm_set1_pd(d * h[k] + (1.0 - d) * cost))));
    }
    memset(&e[0], 0, sizeof(double) << 1);
    for (k = 0; k < 8; k++) {
      r = _mm_loadu_pd(&e[0]);
      _mm_storeu_pd(
          &e[0], _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network_skip[k << 1]),
                                          _mm_set1_pd(z[k]))));
    }
    a[0] = tanh(10.0 * ((a[0] + e[0]) + network_outputBias[0]) / 80.0);
    a[1] = tanh(10.0 * ((a[1] + e[1]) + network_outputBias[1]) / 80.0);
    for (k = 0; k < 3; k++) {
      for (i = 0; i <= 382; i += 2) {
        _mm_storeu_pd(
            &W[i + 384 * k],
            _mm_div_pd(_mm_loadu_pd(&network_input[i + 384 * (k + 4)]),
                       _mm_set1_pd(network_std[k + 4])));
      }
    }
    for (k = 0; k <= 126; k += 2) {
      r = _mm_loadu_pd(&c[k]);
      _mm_storeu_pd(&b_g[k], _mm_mul_pd(r, r));
    }
    __m128d r3;
    __m128d r4;
    r = _mm_loadu_pd(&a[0]);
    _mm_storeu_pd(&B[0],
                  _mm_mul_pd(_mm_set1_pd(10.0),
                             _mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(r, r))));
    cost = 2.0 * options[5];
    e[1] = 0.0;
    e[0] = options[5];
    r = _mm_loadu_pd(&linearRequest[0]);
    r2 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_add_pd(r, r2));
    r1 = _mm_set1_pd(0.0);
    _mm_storeu_pd(&b_dv[0], r1);
    pacAllocatePressure(b_linearRequest, b_dv, cfg->common_bar, cfg->sTrim_bar,
                        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar,
                        cfg->maxDeltaS_bar, cfg->maxDeltaA_bar,
                        cfg->activation_bar, cfg->forcePerBar_N,
                        cfg->enableCommonAdaptation, ep, a, trial);
    r = _mm_loadu_pd(&previous[0]);
    r3 = _mm_set1_pd(step);
    r2 = _mm_sub_pd(r, r3);
    _mm_storeu_pd(&minus[0], r2);
    _mm_storeu_pd(&dv1[0], r2);
    r = _mm_add_pd(r, r3);
    _mm_storeu_pd(&varargin_2[0], r);
    _mm_storeu_pd(&dv2[0], r);
    dv3[0] = fmin(fmax(ep[0], dv1[0]), dv2[0]);
    dv3[1] = fmin(fmax(ep[1], dv1[1]), dv2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&plus[0], r);
    d = previous[2] - step;
    minus[2] = d;
    c_tmp_tmp = previous[2] + step;
    plus[2] = fmin(fmax(ep[2], d), c_tmp_tmp);
    r = _mm_loadu_pd(&linearRequest[0]);
    r2 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_sub_pd(r, r2));
    _mm_storeu_pd(&b_dv[0], r1);
    pacAllocatePressure(b_linearRequest, b_dv, cfg->common_bar, cfg->sTrim_bar,
                        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar,
                        cfg->maxDeltaS_bar, cfg->maxDeltaA_bar,
                        cfg->activation_bar, cfg->forcePerBar_N,
                        cfg->enableCommonAdaptation, ep, a, trial);
    dv3[0] = fmin(fmax(ep[0], minus[0]), varargin_2[0]);
    dv3[1] = fmin(fmax(ep[1], minus[1]), varargin_2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&minus[0], r);
    r2 = _mm_loadu_pd(&plus[0]);
    r4 = _mm_set1_pd(cost);
    _mm_storeu_pd(&seeds[0], _mm_div_pd(_mm_sub_pd(r2, r), r4));
    seeds[2] = (plus[2] - fmin(fmax(ep[2], minus[2]), c_tmp_tmp)) / cost;
    e[0] = 0.0;
    e[1] = options[5];
    r = _mm_loadu_pd(&linearRequest[0]);
    r2 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_add_pd(r, r2));
    _mm_storeu_pd(&b_dv[0], r1);
    pacAllocatePressure(b_linearRequest, b_dv, cfg->common_bar, cfg->sTrim_bar,
                        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar,
                        cfg->maxDeltaS_bar, cfg->maxDeltaA_bar,
                        cfg->activation_bar, cfg->forcePerBar_N,
                        cfg->enableCommonAdaptation, ep, a, trial);
    r = _mm_loadu_pd(&previous[0]);
    r2 = _mm_sub_pd(r, r3);
    _mm_storeu_pd(&minus[0], r2);
    _mm_storeu_pd(&dv1[0], r2);
    r = _mm_add_pd(r, r3);
    _mm_storeu_pd(&varargin_2[0], r);
    _mm_storeu_pd(&dv2[0], r);
    dv3[0] = fmin(fmax(ep[0], dv1[0]), dv2[0]);
    dv3[1] = fmin(fmax(ep[1], dv1[1]), dv2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&plus[0], r);
    minus[2] = d;
    plus[2] = fmin(fmax(ep[2], d), c_tmp_tmp);
    r = _mm_loadu_pd(&linearRequest[0]);
    r2 = _mm_loadu_pd(&e[0]);
    _mm_storeu_pd(&b_linearRequest[0], _mm_sub_pd(r, r2));
    _mm_storeu_pd(&b_dv[0], r1);
    pacAllocatePressure(b_linearRequest, b_dv, cfg->common_bar, cfg->sTrim_bar,
                        cfg->aTrim_bar, cfg->minimum_bar, cfg->maximum_bar,
                        cfg->maxDeltaS_bar, cfg->maxDeltaA_bar,
                        cfg->activation_bar, cfg->forcePerBar_N,
                        cfg->enableCommonAdaptation, ep, a, trial);
    dv3[0] = fmin(fmax(ep[0], minus[0]), varargin_2[0]);
    dv3[1] = fmin(fmax(ep[1], minus[1]), varargin_2[1]);
    r = _mm_loadu_pd(&dv3[0]);
    _mm_storeu_pd(&minus[0], r);
    r2 = _mm_loadu_pd(&plus[0]);
    _mm_storeu_pd(&seeds[3], _mm_div_pd(_mm_sub_pd(r2, r), r4));
    seeds[5] = (plus[2] - fmin(fmax(ep[2], minus[2]), c_tmp_tmp)) / cost;
    memset(&b_network_output[0], 0, 6U * sizeof(double));
    r2 = _mm_loadu_pd(&B[0]);
    for (i = 0; i < 3; i++) {
      B_tmp = i << 1;
      for (k = 0; k < 128; k++) {
        int i1;
        r = _mm_loadu_pd(&b_network_output[B_tmp]);
        i1 = k + 384 * i;
        cost = next_tmp[k];
        c_tmp_tmp = b_r[k];
        _mm_storeu_pd(
            &b_network_output[B_tmp],
            _mm_add_pd(
                r, _mm_mul_pd(
                       _mm_loadu_pd(&network_output[k << 1]),
                       _mm_set1_pd(
                           cost * ((1.0 - b_g[k]) *
                                   (W[i1 + 256] +
                                    c_tmp[k] * (c_tmp_tmp * (1.0 - c_tmp_tmp) *
                                                W[i1]))) +
                           (h[k] - c[k]) * (g[k] * cost * W[i1 + 128])))));
      }
      B_tmp = i << 1;
      r = _mm_loadu_pd(&b_network_output[B_tmp]);
      _mm_storeu_pd(
          &y[B_tmp],
          _mm_add_pd(
              _mm_mul_pd(
                  r2,
                  _mm_add_pd(
                      r, _mm_div_pd(_mm_loadu_pd(&network_skip[(i + 4) << 1]),
                                    _mm_set1_pd(network_std[i + 4])))),
              _mm_loadu_pd(&dv[B_tmp])));
    }
    memset(&J[0], 0, sizeof(double) << 2);
    for (k = 0; k < 2; k++) {
      r = _mm_loadu_pd(&y[0]);
      B_tmp = k << 1;
      r2 = _mm_loadu_pd(&J[B_tmp]);
      _mm_storeu_pd(&J[B_tmp],
                    _mm_add_pd(r2, _mm_mul_pd(r, _mm_set1_pd(seeds[3 * k]))));
      r = _mm_loadu_pd(&y[2]);
      r2 = _mm_loadu_pd(&J[B_tmp]);
      _mm_storeu_pd(
          &J[B_tmp],
          _mm_add_pd(r2, _mm_mul_pd(r, _mm_set1_pd(seeds[3 * k + 1]))));
      r = _mm_loadu_pd(&y[4]);
      r2 = _mm_loadu_pd(&J[B_tmp]);
      _mm_storeu_pd(
          &J[B_tmp],
          _mm_add_pd(r2, _mm_mul_pd(r, _mm_set1_pd(seeds[3 * k + 2]))));
    }
    y_tmp[0] = J[0];
    y_tmp[1] = J[2];
    y_tmp[2] = J[1];
    y_tmp[3] = J[3];
    memset(&B[0], 0, sizeof(double) << 1);
    cost = J[0];
    d = J[1];
    c_tmp_tmp = J[2];
    d1 = J[3];
    d2 = B[0];
    d3 = B[1];
    for (k = 0; k < 2; k++) {
      double d4;
      double d5;
      d4 = y_tmp[k + 2];
      J[k] = -((y_tmp[k] * cost + d4 * d) + b_y[k]);
      B_tmp = k << 1;
      d5 = b_error[k];
      d2 += y_tmp[B_tmp] * d5;
      J[k + 2] = -((y_tmp[k] * c_tmp_tmp + d4 * d1) + b_y[k + 2]);
      d3 += y_tmp[B_tmp + 1] * d5;
    }
    B[1] = d3;
    B[0] = d2;
    mldivide(J, B, a);
    cost = fmin(1.0, options[6] / fmax(b_norm(a), 1.0E-12));
    c_tmp_tmp = a[0] * cost;
    a[0] = c_tmp_tmp;
    bv[0] = ((!rtIsInf(c_tmp_tmp)) && (!rtIsNaN(c_tmp_tmp)));
    c_tmp_tmp = a[1] * cost;
    a[1] = c_tmp_tmp;
    bv[1] = ((!rtIsInf(c_tmp_tmp)) && (!rtIsNaN(c_tmp_tmp)));
    if ((!all(bv)) || (b_norm(a) < 1.0E-10)) {
      exitg1 = true;
    } else {
      boolean_T exitg2;
      limited = false;
      B_tmp = 0;
      exitg2 = false;
      while ((!exitg2) && (B_tmp < 4)) {
        r = _mm_loadu_pd(&a[0]);
        r2 = _mm_loadu_pd(&linearRequest[0]);
        _mm_storeu_pd(
            &trial[0],
            _mm_add_pd(r2,
                       _mm_mul_pd(_mm_set1_pd(rt_powd_snf(0.5, B_tmp)), r)));
        candidate(trial, request, previous, step, measured, featureDt, h,
                  network_input, network_recurrent, network_bias,
                  network_output, network_skip, network_outputBias,
                  network_mean, network_std, cfg, minus, e, ep);
        r = _mm_loadu_pd(&e[0]);
        _mm_storeu_pd(&B[0], _mm_mul_pd(r, r));
        cost = B[0] + B[1];
        if ((!rtIsInf(cost)) && (!rtIsNaN(cost)) && (cost < best - 1.0E-14)) {
          best = cost;
          pressure[0] = minus[0];
          pressure[1] = minus[1];
          pressure[2] = minus[2];
          b_error[0] = e[0];
          b_error[1] = e[1];
          target[0] = ep[0];
          target[1] = ep[1];
          target[2] = ep[2];
          linearRequest[0] = trial[0];
          linearRequest[1] = trial[1];
          limited = true;
          exitg2 = true;
        } else {
          B_tmp++;
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
    boolean_T b_minus[3];
    b_minus[0] = (fabs(pressure[0] - target[0]) > 1.0E-9);
    b_minus[1] = (fabs(pressure[1] - target[1]) > 1.0E-9);
    b_minus[2] = (fabs(pressure[2] - target[2]) > 1.0E-9);
    if (any(b_minus)) {
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
