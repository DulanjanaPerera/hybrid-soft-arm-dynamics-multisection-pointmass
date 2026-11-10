/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacPressureStep.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

/* Include Files */
#include "pacPressureStep.h"
#include "all.h"
#include "any.h"
#include "pacAllocatePressure.h"
#include "pdfbResidualInverseBlock_internal_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
/*
 * PACPRESSURESTEP Discrete command supervisor, no I/O. Call at 0.05 s.
 *  status: 0 off, 1 trim/wait anchor, 2 blend, 3 control, 4 invalid hold, 5
 * fault. equivalentForce/residual use NOMINAL activation model and calibrated
 * bias. Invalid grace holds pressure; sustained invalid data latches then slews
 * zero. enable=false slews zero; external physical Kill remains
 * separate/primary.
 *
 * Arguments    : const double requested[2]
 *                const double staticForce[2]
 *                const double speed[2]
 *                boolean_T valid
 *                double sampleTime
 *                double t
 *                boolean_T enable
 *                const struct_T *cfg
 *                b_struct_T *state
 *                double pressure[3]
 *                double equivalentForce[2]
 *                double residual[2]
 *                unsigned char *status
 *                boolean_T *thresholdActive
 *                double effectiveRequest[2]
 * Return Type  : boolean_T
 */
boolean_T pacPressureStep(const double requested[2],
                          const double staticForce[2], const double speed[2],
                          boolean_T valid, double sampleTime, double t,
                          boolean_T enable, const struct_T *cfg,
                          b_struct_T *state, double pressure[3],
                          double equivalentForce[2], double residual[2],
                          unsigned char *status, boolean_T *thresholdActive,
                          double effectiveRequest[2])
{
  double ep[3];
  double trim[3];
  double b_dv[2];
  double b_y[2];
  double dv1[2];
  boolean_T limited;
  equivalentForce[0] = 0.0;
  residual[0] = 0.0;
  effectiveRequest[0] = 0.0;
  equivalentForce[1] = 0.0;
  residual[1] = 0.0;
  effectiveRequest[1] = 0.0;
  limited = false;
  *thresholdActive = false;
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
    __m128d r;
    __m128d r1;
    double target[3];
    double dt;
    double step;
    boolean_T b_ep[3];
    boolean_T y[2];
    boolean_T fresh;
    if ((state->lastTime >= 0.0) && (t < state->lastTime)) {
      /*  Clock reversal is a fault, not automatic runtime re-arming. */
      state->fault = true;
      state->lastTime = t;
    }
    if (state->lastTime < 0.0) {
      state->lastTime = t;
    }
    dt = t - state->lastTime;
    if ((!rtIsInf(sampleTime)) && (!rtIsNaN(sampleTime)) &&
        (sampleTime <= t + 1.0E-9) &&
        (t - sampleTime <= cfg->maxSampleAge_s + 1.0E-9)) {
      fresh = true;
    } else {
      fresh = false;
    }
    if (valid && fresh) {
      y[0] = ((!rtIsInf(requested[0])) && (!rtIsNaN(requested[0])));
      y[1] = ((!rtIsInf(requested[1])) && (!rtIsNaN(requested[1])));
      if (d_all(y)) {
        y[0] = ((!rtIsInf(staticForce[0])) && (!rtIsNaN(staticForce[0])));
        y[1] = ((!rtIsInf(staticForce[1])) && (!rtIsNaN(staticForce[1])));
        if (d_all(y)) {
          y[0] = ((!rtIsInf(speed[0])) && (!rtIsNaN(speed[0])));
          y[1] = ((!rtIsInf(speed[1])) && (!rtIsNaN(speed[1])));
          if (d_all(y)) {
            fresh = true;
          } else {
            fresh = false;
          }
        } else {
          fresh = false;
        }
      } else {
        fresh = false;
      }
    } else {
      fresh = false;
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
      /*  A later enable starts a new trim calibration, not an old anchor. */
      state->anchorReady = false;
      state->engageTime = -1.0;
    } else if (!fresh) {
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
      if ((!state->anchorReady) && (t >= cfg->feedbackStart_s)) {
        b_ep[0] = (fabs(state->pressure[0] - trim[0]) < 1.0E-6);
        b_ep[1] = (fabs(state->pressure[1] - trim[1]) < 1.0E-6);
        b_ep[2] = (fabs(state->pressure[2] - trim[2]) < 1.0E-6);
        if (all(b_ep)) {
          y[0] = (fabs(speed[0]) <= cfg->maxAnchorSpeed_mPerS);
          y[1] = (fabs(speed[1]) <= cfg->maxAnchorSpeed_mPerS);
          if (d_all(y)) {
            state->anchorForce[0] = staticForce[0];
            state->anchorForce[1] = staticForce[1];
            r = _mm_loadu_pd(&trim[0]);
            r = _mm_sub_pd(r, _mm_set1_pd(cfg->activation_bar));
            _mm_storeu_pd(&trim[0], r);
            _mm_storeu_pd(&b_y[0], r);
            b_y[0] = fmax(b_y[0], 0.0);
            b_y[1] = fmax(b_y[1], 0.0);
            r = _mm_loadu_pd(&b_y[0]);
            _mm_storeu_pd(&ep[0], r);
            state->forceBias[0] =
                cfg->forcePerBar_N * (ep[1] - ep[0]) - staticForce[0];
            state->forceBias[1] =
                cfg->forcePerBar_N *
                    (fmax(trim[2] - cfg->activation_bar, 0.0) - ep[0]) -
                staticForce[1];
            state->anchorReady = true;
            state->engageTime = t;
          }
        }
      }
      if (state->anchorReady) {
        double a__2[2];
        step = fmin(fmax((t - state->engageTime) / cfg->blendTime_s, 0.0), 1.0);
        r = _mm_loadu_pd(&state->anchorForce[0]);
        _mm_storeu_pd(
            &effectiveRequest[0],
            _mm_add_pd(r,
                       _mm_mul_pd(_mm_set1_pd(step * step * (3.0 - 2.0 * step)),
                                  _mm_sub_pd(_mm_loadu_pd(&requested[0]), r))));
        limited = pacAllocatePressure(
            effectiveRequest, state->anchorForce, state->forceBias,
            cfg->common_bar, cfg->sTrim_bar, cfg->aTrim_bar, cfg->minimum_bar,
            cfg->maximum_bar, cfg->maxDeltaS_bar, cfg->maxDeltaA_bar,
            cfg->activation_bar, cfg->forcePerBar_N,
            cfg->enableCommonAdaptation, target, b_y, a__2);
        if (step < 1.0) {
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
      b_y[0] = fmin(fmax(target[0], b_dv[0]), dv1[0]);
      b_y[1] = fmin(fmax(target[1], b_dv[1]), dv1[1]);
      r = _mm_loadu_pd(&b_y[0]);
      _mm_storeu_pd(&pressure[0], r);
      _mm_storeu_pd(&state->pressure[0], r);
      dt = state->pressure[2];
      step = fmin(fmax(target[2], dt - step), dt + step);
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
      b_ep[0] = (fabs(pressure[0] - target[0]) > 1.0E-9);
      b_ep[1] = (fabs(pressure[1] - target[1]) > 1.0E-9);
      b_ep[2] = (fabs(pressure[2] - target[2]) > 1.0E-9);
      if (any(b_ep)) {
        limited = true;
      }
    }
    b_ep[0] = (pressure[0] > cfg->activation_bar);
    ep[0] = fmax(pressure[0] - cfg->activation_bar, 0.0);
    b_ep[1] = (pressure[1] > cfg->activation_bar);
    b_ep[2] = (pressure[2] > cfg->activation_bar);
    *thresholdActive = all(b_ep);
    if (state->anchorReady) {
      equivalentForce[0] =
          cfg->forcePerBar_N *
              (fmax(pressure[1] - cfg->activation_bar, 0.0) - ep[0]) -
          state->forceBias[0];
      equivalentForce[1] =
          cfg->forcePerBar_N *
              (fmax(pressure[2] - cfg->activation_bar, 0.0) - ep[0]) -
          state->forceBias[1];
      if ((*status == 2) || (*status == 3)) {
        r = _mm_loadu_pd(&equivalentForce[0]);
        r1 = _mm_loadu_pd(&effectiveRequest[0]);
        _mm_storeu_pd(&residual[0], _mm_sub_pd(r, r1));
      } else {
        /*  No valid current request to compare during hold/fault/disabled. */
        residual[0] = rtNaN;
        residual[1] = rtNaN;
      }
    } else {
      equivalentForce[0] = rtNaN;
      residual[0] = rtNaN;
      equivalentForce[1] = rtNaN;
      residual[1] = rtNaN;
    }
  }
  return limited;
}

/*
 * File trailer for pacPressureStep.c
 *
 * [EOF]
 */
