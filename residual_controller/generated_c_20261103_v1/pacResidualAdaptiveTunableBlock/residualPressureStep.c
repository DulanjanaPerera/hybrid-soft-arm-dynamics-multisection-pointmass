/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualPressureStep.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "residualPressureStep.h"
#include "all.h"
#include "allOrAny.h"
#include "norm.h"
#include "pacPressureStep.h"
#include "pacResidualAdaptiveTunableBlock_data.h"
#include "pacResidualAdaptiveTunableBlock_internal_types.h"
#include "pacResidualAdaptiveTunableBlock_types.h"
#include "residualAllocatePressure.h"
#include "residualForceMap.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * RESIDUALPRESSURESTEP Existing PAC supervisor with model-aware inverse/anchor.
 *  Off delegates EXACTLY to pacPressureStep. On reports augmented nominal
 * force. Short missing samples hold memory (no invented observations); >maxGap
 * resets.
 *
 * Arguments    : const double requested[2]
 *                const double staticForce[2]
 *                const double measured[6]
 *                boolean_T valid
 *                double sampleTime
 *                double t
 *                boolean_T enable
 *                const e_struct_T *cfg
 *                b_struct_T *state
 *                c_struct_T *b_ns
 *                const f_struct_T *network
 *                const double options[8]
 *                boolean_T residualEnabled
 *                const double shift[2]
 *                double pressure[3]
 *                double equivalent[2]
 *                double b_error[2]
 *                unsigned char *status
 *                boolean_T *active
 *                double effective[2]
 *                struct0_T *d
 * Return Type  : boolean_T
 */
boolean_T residualPressureStep(
    const double requested[2], const double staticForce[2],
    const double measured[6], boolean_T valid, double sampleTime, double t,
    boolean_T enable, const e_struct_T *cfg, b_struct_T *state,
    c_struct_T *b_ns, const f_struct_T *network, const double options[8],
    boolean_T residualEnabled, const double shift[2], double pressure[3],
    double equivalent[2], double b_error[2], unsigned char *status,
    boolean_T *active, double effective[2], struct0_T *d)
{
  __m128d r;
  double u[384];
  double w[384];
  double x[130];
  double g[128];
  double z[8];
  double b_x[5];
  double trim[3];
  double b_a[2];
  double b_dv[2];
  double b_network[2];
  double dv1[2];
  double dv2[2];
  double dv3[2];
  double b_d;
  double d1;
  double d3;
  double dt;
  double step;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int k;
  boolean_T a[2];
  boolean_T limited;
  boolean_T newSample;
  boolean_T physical;
  d->advanced = false;
  d->contractValid = true;
  d->rawForce_N[0] = 0.0;
  d->effectiveResidual_N[0] = 0.0;
  d->inverseError_N[0] = 0.0;
  d->rawForce_N[1] = 0.0;
  d->effectiveResidual_N[1] = 0.0;
  d->inverseError_N[1] = 0.0;
  d->iterations = 0U;
  d->converged = false;
  memset(&d->parameterCorrection[0], 0, 10U * sizeof(double));
  d->enabled = residualEnabled;
  d->units = network->units;
  d->gravityReexpression_N[0] = shift[0];
  d->gravityReexpression_N[1] = shift[1];
  newSample = false;
  d->historyGapHeld = false;
  if (!b_ns->initialized) {
    b_ns->initialized = true;
    b_ns->units = network->units;
    b_ns->mode = residualEnabled;
  } else if ((b_ns->units != network->units) ||
             (b_ns->mode != residualEnabled)) {
    if (enable && state->anchorReady) {
      state->fault = true;
    } else {
      memset(&b_ns->h[0], 0, 128U * sizeof(double));
      memset(&b_ns->pre[0], 0, 128U * sizeof(double));
      b_ns->lastSample = -1.0;
      b_ns->featureDt = 0.05;
      b_ns->warm_s = 0.0;
      b_ns->initialized = true;
      b_ns->units = network->units;
      b_ns->mode = residualEnabled;
    }
  }
  if ((!rtIsInf(sampleTime)) && (!rtIsNaN(sampleTime)) &&
      ((!rtIsInf(t)) && (!rtIsNaN(t))) && (sampleTime <= t + 1.0E-9) &&
      (t - sampleTime <= cfg->maxSampleAge_s + 1.0E-9)) {
    physical = true;
  } else {
    physical = false;
  }
  if (valid && physical) {
    boolean_T bv[4];
    bv[0] = ((!rtIsInf(measured[0])) && (!rtIsNaN(measured[0])));
    bv[1] = ((!rtIsInf(measured[1])) && (!rtIsNaN(measured[1])));
    bv[2] = ((!rtIsInf(measured[2])) && (!rtIsNaN(measured[2])));
    bv[3] = ((!rtIsInf(measured[3])) && (!rtIsNaN(measured[3])));
    if (e_all(bv)) {
      a[0] = ((!rtIsInf(shift[0])) && (!rtIsNaN(shift[0])));
      a[1] = ((!rtIsInf(shift[1])) && (!rtIsNaN(shift[1])));
      if (c_all(a)) {
        physical = true;
      } else {
        physical = false;
      }
    } else {
      physical = false;
    }
  } else {
    physical = false;
  }
  if (physical && (sampleTime > b_ns->lastSample)) {
    newSample = true;
    dt = options[7];
    if (b_ns->lastSample >= 0.0) {
      dt = sampleTime - b_ns->lastSample;
    }
    if ((b_ns->lastSample >= 0.0) && (dt > options[4] + 1.0E-9)) {
      memset(&b_ns->h[0], 0, 128U * sizeof(double));
      b_ns->warm_s = 0.0;
    } else if (b_ns->lastSample >= 0.0) {
      b_ns->warm_s += dt;
      d->historyGapHeld = (dt > 1.5 * options[7]);
    }
    memcpy(&b_ns->pre[0], &b_ns->h[0], 128U * sizeof(double));
    b_ns->featureDt = dt;
    b_ns->lastSample = sampleTime;
  }
  d->ready = (b_ns->warm_s >= options[3] - 1.0E-9);
  d->warm_s = b_ns->warm_s;
  if (physical && residualEnabled) {
    boolean_T bv1[130];
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
    r = _mm_loadu_pd(&state->pressure[0]);
    _mm_storeu_pd(&z[4],
                  _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&network->mean[4])),
                             _mm_loadu_pd(&network->std[4])));
    z[6] = (state->pressure[2] - network->mean[6]) / network->std[6];
    z[7] = (b_ns->featureDt - network->mean[7]) / network->std[7];
#pragma omp parallel for num_threads(omp_get_max_threads()) private(b_d, i1)

    for (i = 0; i < 384; i++) {
      b_d = 0.0;
      for (i1 = 0; i1 < 8; i1++) {
        b_d += network->input[i + 384 * i1] * z[i1];
      }
      u[i] = b_d + network->bias[i];
    }
    memset(&w[0], 0, 384U * sizeof(double));
    for (k = 0; k < 128; k++) {
      step = b_ns->pre[k];
      for (i2 = 0; i2 <= 382; i2 += 2) {
        r = _mm_loadu_pd(&w[i2]);
        _mm_storeu_pd(
            &w[i2],
            _mm_add_pd(
                r, _mm_mul_pd(_mm_loadu_pd(&network->recurrent[i2 + 384 * k]),
                              _mm_set1_pd(step))));
      }
    }
    memset(&b_a[0], 0, sizeof(double) << 1);
    for (k = 0; k < 128; k++) {
      step = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
      step = step * b_ns->pre[k] +
             (1.0 - step) *
                 tanh(u[k + 256] +
                      0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0) * w[k + 256]);
      g[k] = step;
      r = _mm_loadu_pd(&b_a[0]);
      _mm_storeu_pd(
          &b_a[0],
          _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->output[k << 1]),
                                   _mm_set1_pd(step))));
    }
    memset(&b_network[0], 0, sizeof(double) << 1);
    for (k = 0; k < 8; k++) {
      r = _mm_loadu_pd(&b_network[0]);
      _mm_storeu_pd(
          &b_network[0],
          _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->skip[k << 1]),
                                   _mm_set1_pd(z[k]))));
    }
    b_a[0] =
        tanh(10.0 * ((b_a[0] + b_network[0]) + network->outputBias[0]) / 80.0);
    b_a[1] =
        tanh(10.0 * ((b_a[1] + b_network[1]) + network->outputBias[1]) / 80.0);
    memcpy(&x[0], &g[0], 128U * sizeof(double));
    r = _mm_loadu_pd(&b_a[0]);
    _mm_storeu_pd(&x[128], _mm_mul_pd(_mm_set1_pd(80.0), r));
    for (k = 0; k < 130; k++) {
      step = x[k];
      bv1[k] = ((!rtIsInf(step)) && (!rtIsNaN(step)));
    }
    if (!f_all(bv1)) {
      physical = false;
      newSample = false;
      d->contractValid = false;
    }
  }
  if (!residualEnabled) {
    limited = pacPressureStep(requested, staticForce, &measured[2], valid,
                              sampleTime, t, enable, cfg, state, pressure,
                              equivalent, b_error, status, active, effective);
    if (newSample) {
      boolean_T y[3];
      y[0] = ((!rtIsInf(pressure[0])) && (!rtIsNaN(pressure[0])));
      y[1] = ((!rtIsInf(pressure[1])) && (!rtIsNaN(pressure[1])));
      y[2] = ((!rtIsInf(pressure[2])) && (!rtIsNaN(pressure[2])));
      if (all(y)) {
        boolean_T bv2[128];
        /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
         * evaluation. */
        /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history
         * fixed. */
        _mm_storeu_pd(&z[0],
                      _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                            _mm_loadu_pd(&network->mean[0])),
                                 _mm_loadu_pd(&network->std[0])));
        _mm_storeu_pd(&z[2],
                      _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                            _mm_loadu_pd(&network->mean[2])),
                                 _mm_loadu_pd(&network->std[2])));
        r = _mm_loadu_pd(&pressure[0]);
        _mm_storeu_pd(&z[4],
                      _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&network->mean[4])),
                                 _mm_loadu_pd(&network->std[4])));
        z[6] = (pressure[2] - network->mean[6]) / network->std[6];
        z[7] = (b_ns->featureDt - network->mean[7]) / network->std[7];
#pragma omp parallel for num_threads(omp_get_max_threads()) private(d1, i4)

        for (i3 = 0; i3 < 384; i3++) {
          d1 = 0.0;
          for (i4 = 0; i4 < 8; i4++) {
            d1 += network->input[i3 + 384 * i4] * z[i4];
          }
          u[i3] = d1 + network->bias[i3];
        }
        memset(&w[0], 0, 384U * sizeof(double));
        for (k = 0; k < 128; k++) {
          step = b_ns->pre[k];
          for (i2 = 0; i2 <= 382; i2 += 2) {
            r = _mm_loadu_pd(&w[i2]);
            _mm_storeu_pd(
                &w[i2],
                _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(
                                             &network->recurrent[i2 + 384 * k]),
                                         _mm_set1_pd(step))));
          }
        }
        memset(&b_a[0], 0, sizeof(double) << 1);
        for (k = 0; k < 128; k++) {
          step = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
          step = step * b_ns->pre[k] +
                 (1.0 - step) *
                     tanh(u[k + 256] +
                          0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0) * w[k + 256]);
          g[k] = step;
          r = _mm_loadu_pd(&b_a[0]);
          _mm_storeu_pd(
              &b_a[0],
              _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->output[k << 1]),
                                       _mm_set1_pd(step))));
        }
        memset(&b_network[0], 0, sizeof(double) << 1);
        for (k = 0; k < 8; k++) {
          r = _mm_loadu_pd(&b_network[0]);
          _mm_storeu_pd(
              &b_network[0],
              _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->skip[k << 1]),
                                       _mm_set1_pd(z[k]))));
        }
        d->rawForce_N[0] =
            80.0 *
            tanh(10.0 * ((b_a[0] + b_network[0]) + network->outputBias[0]) /
                 80.0);
        d->rawForce_N[1] =
            80.0 *
            tanh(10.0 * ((b_a[1] + b_network[1]) + network->outputBias[1]) /
                 80.0);
        for (k = 0; k < 128; k++) {
          step = g[k];
          bv2[k] = ((!rtIsInf(step)) && (!rtIsNaN(step)));
        }
        if (g_all(bv2)) {
          memcpy(&b_ns->h[0], &g[0], 128U * sizeof(double));
          d->advanced = true;
        } else {
          memset(&b_ns->h[0], 0, 128U * sizeof(double));
          b_ns->warm_s = 0.0;
        }
      }
    }
    d->inverseError_N[0] = b_error[0];
    a[0] = ((!rtIsInf(b_error[0])) && (!rtIsNaN(b_error[0])));
    d->inverseError_N[1] = b_error[1];
    a[1] = ((!rtIsInf(b_error[1])) && (!rtIsNaN(b_error[1])));
    if (c_all(a) && (b_norm(b_error) <= options[2])) {
      d->converged = true;
    } else {
      d->converged = false;
    }
  } else {
    equivalent[0] = 0.0;
    b_error[0] = 0.0;
    effective[0] = 0.0;
    equivalent[1] = 0.0;
    b_error[1] = 0.0;
    effective[1] = 0.0;
    limited = false;
    *active = false;
    if (rtIsInf(t) || rtIsNaN(t) || (t < 0.0)) {
      state->fault = true;
      pressure[0] = 0.0;
      state->pressure[0] = 0.0;
      pressure[1] = 0.0;
      state->pressure[1] = 0.0;
      pressure[2] = 0.0;
      state->pressure[2] = 0.0;
      state->lastTime = -1.0;
      *status = 5U;
    } else {
      __m128d r1;
      double target[3];
      double d2;
      boolean_T y[3];
      boolean_T usable;
      if ((state->lastTime >= 0.0) && (t < state->lastTime)) {
        state->fault = true;
        state->lastTime = t;
      }
      if (state->lastTime < 0.0) {
        state->lastTime = t;
      }
      dt = t - state->lastTime;
      if (d->contractValid && (fabs(cfg->forcePerBar_N - 13.273228961416876) <=
                               1.3273228961416876E-9)) {
        d->contractValid = true;
      } else {
        d->contractValid = false;
      }
      if (physical && newSample) {
        a[0] = ((!rtIsInf(requested[0])) && (!rtIsNaN(requested[0])));
        a[1] = ((!rtIsInf(requested[1])) && (!rtIsNaN(requested[1])));
        if (c_all(a)) {
          a[0] = ((!rtIsInf(staticForce[0])) && (!rtIsNaN(staticForce[0])));
          a[1] = ((!rtIsInf(staticForce[1])) && (!rtIsNaN(staticForce[1])));
          if (c_all(a) && d->contractValid) {
            usable = true;
          } else {
            usable = false;
          }
        } else {
          usable = false;
        }
      } else {
        usable = false;
      }
      if (state->anchorReady && (!d->ready)) {
        usable = false;
      }
      trim[0] = cfg->common_bar - 2.0 * cfg->sTrim_bar / 3.0;
      step = cfg->common_bar + cfg->sTrim_bar / 3.0;
      trim[1] = step - cfg->aTrim_bar;
      trim[2] = step + cfg->aTrim_bar;
      target[0] = 0.0;
      target[1] = 0.0;
      target[2] = 0.0;
      if (state->fault) {
        *status = 5U;
      } else if ((!enable) || (t < cfg->pressureStart_s) ||
                 (t >= cfg->stopTime_s)) {
        *status = 0U;
        state->lastGood = -1.0;
        state->anchorReady = false;
        state->engageTime = -1.0;
      } else if (!usable) {
        if ((state->lastGood >= 0.0) &&
            (t - state->lastGood <= cfg->invalidGrace_s + 1.0E-9)) {
          target[0] = state->pressure[0];
          target[1] = state->pressure[1];
          target[2] = state->pressure[2];
          *status = 4U;
        } else {
          state->fault = true;
          *status = 5U;
        }
      } else {
        state->lastGood = t;
        target[0] = trim[0];
        target[1] = trim[1];
        target[2] = trim[2];
        *status = 1U;
        if ((!state->anchorReady) && d->ready && (t >= cfg->feedbackStart_s)) {
          y[0] = (fabs(state->pressure[0] - trim[0]) < 1.0E-6);
          y[1] = (fabs(state->pressure[1] - trim[1]) < 1.0E-6);
          y[2] = (fabs(state->pressure[2] - trim[2]) < 1.0E-6);
          if (all(y)) {
            a[0] = (fabs(measured[2]) <= cfg->maxAnchorSpeed_mPerS);
            a[1] = (fabs(measured[3]) <= cfg->maxAnchorSpeed_mPerS);
            if (c_all(a)) {
              residualForceMap(trim, measured, b_ns->featureDt, b_ns->pre,
                               network->input, network->recurrent,
                               network->bias, network->output, network->skip,
                               network->outputBias, network->mean, network->std,
                               cfg->activation_bar, cfg->forcePerBar_N, shift,
                               b_a, state->residualAnchor);
              state->anchorReady = true;
              state->engageTime = t;
              r1 = _mm_loadu_pd(&staticForce[0]);
              _mm_storeu_pd(&state->anchorForce[0], r1);
              r = _mm_loadu_pd(&b_a[0]);
              _mm_storeu_pd(&state->forceBias[0], _mm_sub_pd(r, r1));
              r = _mm_loadu_pd(&state->residualAnchor[0]);
              _mm_storeu_pd(&state->linearRequest[0], _mm_sub_pd(r1, r));
              _mm_storeu_pd(&state->lastEffective[0], r1);
            }
          }
        }
        if (state->anchorReady) {
          boolean_T bv3[5];
          step =
              fmin(fmax((t - state->engageTime) / cfg->blendTime_s, 0.0), 1.0);
          r = _mm_loadu_pd(&state->anchorForce[0]);
          _mm_storeu_pd(
              &effective[0],
              _mm_add_pd(
                  r, _mm_mul_pd(_mm_set1_pd(step * step * (3.0 - 2.0 * step)),
                                _mm_sub_pd(_mm_loadu_pd(&requested[0]), r))));
          for (k = 0; k < 2; k++) {
            b_a[k] = state->linearRequest[k];
          }
          limited = residualAllocatePressure(
              effective, state->anchorForce, state->residualAnchor,
              state->forceBias, state->pressure, b_a, state->lastEffective,
              cfg->slew_barPerS * fmax(dt, 0.0), measured, b_ns->featureDt,
              b_ns->pre, network, cfg, shift, options, target, b_error,
              &d->iterations, state->linearRequest, trim);
          state->lastEffective[0] = effective[0];
          state->lastEffective[1] = effective[1];
          b_x[0] = target[0];
          b_x[1] = target[1];
          b_x[2] = target[2];
          b_x[3] = b_error[0];
          b_x[4] = b_error[1];
          for (k = 0; k < 5; k++) {
            d2 = b_x[k];
            bv3[k] = ((!rtIsInf(d2)) && (!rtIsNaN(d2)));
          }
          if (!h_all(bv3)) {
            state->fault = true;
            *status = 5U;
            target[0] = 0.0;
            target[1] = 0.0;
            target[2] = 0.0;
          } else if (step < 1.0) {
            *status = 2U;
          } else {
            *status = 3U;
          }
        }
      }
      if (dt > 0.0) {
        step = cfg->slew_barPerS * dt;
        r = _mm_loadu_pd(&state->pressure[0]);
        r1 = _mm_set1_pd(step);
        _mm_storeu_pd(&b_dv[0], _mm_sub_pd(r, r1));
        _mm_storeu_pd(&dv1[0], _mm_add_pd(r, r1));
        dv2[0] = fmin(fmax(target[0], b_dv[0]), dv1[0]);
        dv2[1] = fmin(fmax(target[1], b_dv[1]), dv1[1]);
        r = _mm_loadu_pd(&dv2[0]);
        _mm_storeu_pd(&pressure[0], r);
        _mm_storeu_pd(&state->pressure[0], r);
        d2 = state->pressure[2];
        step = fmin(fmax(target[2], d2 - step), d2 + step);
        pressure[2] = step;
        state->pressure[2] = step;
        state->lastTime = t;
      } else {
        pressure[0] = state->pressure[0];
        pressure[1] = state->pressure[1];
        pressure[2] = state->pressure[2];
      }
      if (limited) {
        limited = true;
      } else {
        boolean_T y_data[10];
        y_data[0] = (fabs(pressure[0] - target[0]) > 1.0E-9);
        y_data[1] = (fabs(pressure[1] - target[1]) > 1.0E-9);
        y_data[2] = (fabs(pressure[2] - target[2]) > 1.0E-9);
        if (vectorAny(y_data, 3)) {
          limited = true;
        }
      }
      y[0] = (pressure[0] > cfg->activation_bar);
      y[1] = (pressure[1] > cfg->activation_bar);
      y[2] = (pressure[2] > cfg->activation_bar);
      *active = all(y);
      if (physical) {
        /* RESIDUALFORCEMAP Activation re-expression keeps
         * F_cfg+r_eff=B0*p+r0+deltaG. */
        /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
         * evaluation. */
        /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history
         * fixed. */
        _mm_storeu_pd(&z[0],
                      _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                            _mm_loadu_pd(&network->mean[0])),
                                 _mm_loadu_pd(&network->std[0])));
        _mm_storeu_pd(&z[2],
                      _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&measured[2]),
                                            _mm_loadu_pd(&network->mean[2])),
                                 _mm_loadu_pd(&network->std[2])));
        r = _mm_loadu_pd(&pressure[0]);
        _mm_storeu_pd(&z[4],
                      _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&network->mean[4])),
                                 _mm_loadu_pd(&network->std[4])));
        z[6] = (pressure[2] - network->mean[6]) / network->std[6];
        z[7] = (b_ns->featureDt - network->mean[7]) / network->std[7];
#pragma omp parallel for num_threads(omp_get_max_threads()) private(d3, i6)

        for (i5 = 0; i5 < 384; i5++) {
          d3 = 0.0;
          for (i6 = 0; i6 < 8; i6++) {
            d3 += network->input[i5 + 384 * i6] * z[i6];
          }
          u[i5] = d3 + network->bias[i5];
        }
        memset(&w[0], 0, 384U * sizeof(double));
        for (k = 0; k < 128; k++) {
          step = b_ns->pre[k];
          for (i2 = 0; i2 <= 382; i2 += 2) {
            r = _mm_loadu_pd(&w[i2]);
            _mm_storeu_pd(
                &w[i2],
                _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(
                                             &network->recurrent[i2 + 384 * k]),
                                         _mm_set1_pd(step))));
          }
        }
        memset(&b_a[0], 0, sizeof(double) << 1);
        for (k = 0; k < 128; k++) {
          step = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
          step = step * b_ns->pre[k] +
                 (1.0 - step) *
                     tanh(u[k + 256] +
                          0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0) * w[k + 256]);
          g[k] = step;
          r = _mm_loadu_pd(&b_a[0]);
          _mm_storeu_pd(
              &b_a[0],
              _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->output[k << 1]),
                                       _mm_set1_pd(step))));
        }
        memset(&b_network[0], 0, sizeof(double) << 1);
        for (k = 0; k < 8; k++) {
          r = _mm_loadu_pd(&b_network[0]);
          _mm_storeu_pd(
              &b_network[0],
              _mm_add_pd(r, _mm_mul_pd(_mm_loadu_pd(&network->skip[k << 1]),
                                       _mm_set1_pd(z[k]))));
        }
        __m128d r2;
        b_a[0] =
            80.0 *
            tanh(10.0 * ((b_a[0] + b_network[0]) + network->outputBias[0]) /
                 80.0);
        b_a[1] =
            80.0 *
            tanh(10.0 * ((b_a[1] + b_network[1]) + network->outputBias[1]) /
                 80.0);
        r = _mm_loadu_pd(&pressure[0]);
        _mm_storeu_pd(&dv3[0], _mm_sub_pd(r, _mm_set1_pd(cfg->activation_bar)));
        dv2[0] = fmax(dv3[0], 0.0);
        dv2[1] = fmax(dv3[1], 0.0);
        r = _mm_loadu_pd(&dv2[0]);
        _mm_storeu_pd(&trim[0], r);
        step = pressure[0];
        d2 = pressure[1];
        dt = pressure[2];
        r = _mm_loadu_pd(&dv[2]);
        r = _mm_mul_pd(r, _mm_set1_pd(d2));
        r1 = _mm_add_pd(_mm_set1_pd(-13.273228961416876 * step), r);
        r = _mm_loadu_pd(&dv[4]);
        r = _mm_mul_pd(r, _mm_set1_pd(dt));
        r = _mm_add_pd(r1, r);
        r2 = _mm_loadu_pd(&b_a[0]);
        r = _mm_add_pd(r, r2);
        r1 = _mm_loadu_pd(&shift[0]);
        r = _mm_add_pd(r, r1);
        _mm_storeu_pd(&equivalent[0], r);
        _mm_storeu_pd(&d->rawForce_N[0], r2);
        d->effectiveResidual_N[0] =
            equivalent[0] - cfg->forcePerBar_N * (trim[1] - trim[0]);
        d->effectiveResidual_N[1] =
            equivalent[1] -
            cfg->forcePerBar_N *
                (fmax(pressure[2] - cfg->activation_bar, 0.0) - trim[0]);
        if (newSample) {
          boolean_T bv2[128];
          for (k = 0; k < 128; k++) {
            step = g[k];
            bv2[k] = ((!rtIsInf(step)) && (!rtIsNaN(step)));
          }
          if (g_all(bv2)) {
            memcpy(&b_ns->h[0], &g[0], 128U * sizeof(double));
            d->advanced = true;
          }
        }
      } else {
        equivalent[0] = rtNaN;
        equivalent[1] = rtNaN;
      }
      if (state->anchorReady) {
        r = _mm_loadu_pd(&equivalent[0]);
        r1 = _mm_loadu_pd(&state->forceBias[0]);
        _mm_storeu_pd(&equivalent[0], _mm_sub_pd(r, r1));
        if ((*status == 2) || (*status == 3)) {
          r = _mm_loadu_pd(&equivalent[0]);
          r1 = _mm_loadu_pd(&effective[0]);
          _mm_storeu_pd(&b_error[0], _mm_sub_pd(r, r1));
        } else {
          b_error[0] = rtNaN;
          b_error[1] = rtNaN;
        }
      } else {
        equivalent[0] = rtNaN;
        b_error[0] = rtNaN;
        equivalent[1] = rtNaN;
        b_error[1] = rtNaN;
      }
      d->inverseError_N[0] = b_error[0];
      a[0] = ((!rtIsInf(b_error[0])) && (!rtIsNaN(b_error[0])));
      d->inverseError_N[1] = b_error[1];
      a[1] = ((!rtIsInf(b_error[1])) && (!rtIsNaN(b_error[1])));
      if (c_all(a) && (b_norm(b_error) <= options[2])) {
        d->converged = true;
      } else {
        d->converged = false;
      }
    }
  }
  return limited;
}

/*
 * File trailer for residualPressureStep.c
 *
 * [EOF]
 */
