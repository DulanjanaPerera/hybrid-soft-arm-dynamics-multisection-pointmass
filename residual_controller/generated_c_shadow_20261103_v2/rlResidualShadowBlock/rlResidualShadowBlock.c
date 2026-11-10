/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualShadowBlock.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

/* Include Files */
#include "rlResidualShadowBlock.h"
#include "all.h"
#include "allOrAny.h"
#include "minOrMax.h"
#include "residualShadowFlow.h"
#include "rlResidualShadowBlock_data.h"
#include "rlResidualShadowBlock_initialize.h"
#include "rlResidualShadowBlock_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  double options[11];
  boolean_T configured;
  boolean_T segmentActive;
  double lastCall_sim_s;
  double lastHost_s;
  double lastFrame;
  double segmentStart_s;
  double lastQ[2];
  double velocity[2];
  double history[128];
  double neuralStart_s;
  unsigned int segment;
  boolean_T pending;
  double pendingQ[2];
  double pendingV[2];
  double pendingP[3];
  double pendingR[2];
  double pendingPrediction[6];
  double pendingHost_s;
} struct_T;
#endif /* typedef_struct_T */

/* Variable Definitions */
static struct_T state;

/* Function Definitions */
/*
 * RLRESIDUALSHADOWBLOCK Paste into a separate discrete0.05s MATLAB Function.
 *  Observation ONLY. Never connect any output to pressure/control/Kill inputs.
 *  Output columns: augmented, nominal, constant velocity. See
 * SHADOW_VALIDATION.md.
 *
 * Arguments    : const double lengths3_m[3]
 *                boolean_T orientationValid
 *                boolean_T baseValid
 *                double pollStart_s
 *                double tipFrame
 *                const double pressure_bar[3]
 *                double commandTime_sim_s
 *                double t_sim_s
 *                boolean_T reset
 *                const double shadowOptions[11]
 *                double predictionLength_m[6]
 *                double predictionAtMeasurement_m[6]
 *                double retimedPrediction_m[6]
 *                struct0_T *shadowDiagnostic
 *                struct1_T *inverseContext
 * Return Type  : void
 */
void rlResidualShadowBlock(
    const double lengths3_m[3], boolean_T orientationValid, boolean_T baseValid,
    double pollStart_s, double tipFrame, const double pressure_bar[3],
    double commandTime_sim_s, double t_sim_s, boolean_T reset,
    const double shadowOptions[11], double predictionLength_m[6],
    double predictionAtMeasurement_m[6], double retimedPrediction_m[6],
    struct0_T *shadowDiagnostic, struct1_T *inverseContext)
{
  static double network_recurrent[49152];
  static double t0_errorAtMeasurement_m[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  static double t0_retimedError_m[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double network_input[3072];
  double w[384];
  double network_output[256];
  double x[130];
  double g[128];
  double pre[128];
  double network_skip[16];
  double network_mean[8];
  double network_std[8];
  double z[8];
  double b_y[2];
  double q[2];
  double network_outputBias_idx_0;
  double network_outputBias_idx_1;
  double network_units;
  int i;
  int k;
  signed char indicesForNan[12] = {0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 2, 1};
  if (!isInitialized_rlResidualShadowBlock) {
    rlResidualShadowBlock_initialize();
  }
  for (k = 0; k < 6; k++) {
    int t0_retimedError_m_tmp;
    t0_retimedError_m_tmp = k << 1;
    t0_retimedError_m_tmp = indicesForNan[t0_retimedError_m_tmp + 1] +
                            (indicesForNan[t0_retimedError_m_tmp] << 1);
    t0_retimedError_m[t0_retimedError_m_tmp] = rtNaN;
    t0_errorAtMeasurement_m[t0_retimedError_m_tmp] = rtNaN;
  }
  /* RESIDUALCONTROLLERNETWORK Fixed-size code-generation constants; run packer
   * first. */
  /*  Set units64/128 before a fresh run. Never replace weights during a trial.
   */
  if (shadowOptions[0] == 64.0) {
    memcpy(&network_input[0], &n64.input[0], 3072U * sizeof(double));
    memcpy(&network_recurrent[0], &n64.recurrent[0], 49152U * sizeof(double));
    memcpy(&w[0], &n64.bias[0], 384U * sizeof(double));
    memcpy(&network_output[0], &n64.output[0], 256U * sizeof(double));
    memcpy(&network_skip[0], &n64.skip[0], 16U * sizeof(double));
    network_outputBias_idx_0 = n64.outputBias[0];
    network_outputBias_idx_1 = n64.outputBias[1];
    memcpy(&network_mean[0], &n64.mean[0], 8U * sizeof(double));
    memcpy(&network_std[0], &n64.std[0], 8U * sizeof(double));
    network_units = n64.units;
  } else {
    memcpy(&network_input[0], &n128.input[0], 3072U * sizeof(double));
    memcpy(&network_recurrent[0], &n128.recurrent[0], 49152U * sizeof(double));
    memcpy(&w[0], &n128.bias[0], 384U * sizeof(double));
    memcpy(&network_output[0], &n128.output[0], 256U * sizeof(double));
    memcpy(&network_skip[0], &n128.skip[0], 16U * sizeof(double));
    network_outputBias_idx_0 = n128.outputBias[0];
    network_outputBias_idx_1 = n128.outputBias[1];
    memcpy(&network_mean[0], &n128.mean[0], 8U * sizeof(double));
    memcpy(&network_std[0], &n128.std[0], 8U * sizeof(double));
    network_units = n128.units;
  }
  /* RESIDUALSHADOWSTEP Pure observer state machine. NO actuator/controller
   * outputs. */
  /*  Predictions columns: augmented, nominal, constant velocity. All m. */
  /*  No interpolation, no future measurements/commands; memory commits once. */
  /*  Missing/rejected observations split segments, matching training
   * preprocessing. */
  shadowDiagnostic->status = 0U;
  shadowDiagnostic->accepted = false;
  shadowDiagnostic->advanced = false;
  shadowDiagnostic->ready = false;
  shadowDiagnostic->historyReset = false;
  shadowDiagnostic->predictionValid = false;
  shadowDiagnostic->scoreValid = false;
  shadowDiagnostic->retimedScoreValid = false;
  shadowDiagnostic->inverseContextValid = false;
  shadowDiagnostic->qMeasured_m[0] = rtNaN;
  shadowDiagnostic->velocity_mPerS[0] = rtNaN;
  shadowDiagnostic->rawForce_N[0] = rtNaN;
  shadowDiagnostic->qMeasured_m[1] = rtNaN;
  shadowDiagnostic->velocity_mPerS[1] = rtNaN;
  shadowDiagnostic->rawForce_N[1] = rtNaN;
  for (k = 0; k < 6; k++) {
    predictionLength_m[k] = rtNaN;
    predictionAtMeasurement_m[k] = rtNaN;
    retimedPrediction_m[k] = rtNaN;
    shadowDiagnostic->errorAtMeasurement_m[k] = t0_errorAtMeasurement_m[k];
    shadowDiagnostic->retimedError_m[k] = t0_retimedError_m[k];
  }
  shadowDiagnostic->predictionTarget_host_s = rtNaN;
  shadowDiagnostic->scoreOrigin_host_s = rtNaN;
  shadowDiagnostic->dtPhysical_s = rtNaN;
  shadowDiagnostic->commandAge_sim_s = rtNaN;
  shadowDiagnostic->horizonMismatch_s = rtNaN;
  shadowDiagnostic->pressureChangeSinceOrigin_bar[0] = rtNaN;
  shadowDiagnostic->pressureChangeSinceOrigin_bar[1] = rtNaN;
  shadowDiagnostic->pressureChangeSinceOrigin_bar[2] = rtNaN;
  shadowDiagnostic->velocityWarm_s = 0.0;
  shadowDiagnostic->neuralWarm_s = 0.0;
  shadowDiagnostic->maximumTrainingZ = rtNaN;
  shadowDiagnostic->segment = 0U;
  inverseContext->valid = false;
  inverseContext->measured[0] = 0.0;
  inverseContext->measured[1] = 0.0;
  inverseContext->measured[2] = 0.0;
  inverseContext->measured[3] = 0.0;
  memset(&inverseContext->historyPre[0], 0, 128U * sizeof(double));
  inverseContext->pressure_bar[0] = 0.0;
  inverseContext->pressure_bar[1] = 0.0;
  inverseContext->pressure_bar[2] = 0.0;
  inverseContext->featureDt_s = 0.05;
  inverseContext->horizon_s = 0.05;
  inverseContext->sampleTime_host_s = -1.0;
  inverseContext->units = 128.0;
  shadowDiagnostic->units = network_units;
  shadowDiagnostic->sampleTime_host_s = pollStart_s;
  if (reset) {
    memset(&state.options[0], 0, 11U * sizeof(double));
    state.configured = false;
    state.segmentActive = false;
    state.lastCall_sim_s = -1.0;
    state.lastHost_s = -1.0;
    state.lastFrame = -1.0;
    state.segmentStart_s = -1.0;
    state.lastQ[0] = 0.0;
    state.velocity[0] = 0.0;
    state.lastQ[1] = 0.0;
    state.velocity[1] = 0.0;
    memset(&state.history[0], 0, 128U * sizeof(double));
    state.neuralStart_s = -1.0;
    state.segment = 0U;
    state.pending = false;
    state.pendingQ[0] = 0.0;
    state.pendingV[0] = 0.0;
    state.pendingQ[1] = 0.0;
    state.pendingV[1] = 0.0;
    state.pendingP[0] = 0.0;
    state.pendingP[1] = 0.0;
    state.pendingP[2] = 0.0;
    state.pendingR[0] = 0.0;
    state.pendingR[1] = 0.0;
    for (k = 0; k < 6; k++) {
      state.pendingPrediction[k] = 0.0;
    }
    state.pendingHost_s = -1.0;
  } else {
    boolean_T y_data[11];
    boolean_T guard1;
    guard1 = false;
    if (!state.configured) {
      guard1 = true;
    } else {
      for (k = 0; k < 11; k++) {
        y_data[k] = (state.options[k] != shadowOptions[k]);
      }
      if (vectorAny(y_data, 11)) {
        guard1 = true;
      }
    }
    if (guard1) {
      state.segmentActive = false;
      state.lastCall_sim_s = -1.0;
      state.lastHost_s = -1.0;
      state.lastFrame = -1.0;
      state.segmentStart_s = -1.0;
      state.lastQ[0] = 0.0;
      state.velocity[0] = 0.0;
      state.lastQ[1] = 0.0;
      state.velocity[1] = 0.0;
      memset(&state.history[0], 0, 128U * sizeof(double));
      state.neuralStart_s = -1.0;
      state.segment = 0U;
      state.pending = false;
      state.pendingQ[0] = 0.0;
      state.pendingV[0] = 0.0;
      state.pendingQ[1] = 0.0;
      state.pendingV[1] = 0.0;
      state.pendingP[0] = 0.0;
      state.pendingP[1] = 0.0;
      state.pendingP[2] = 0.0;
      state.pendingR[0] = 0.0;
      state.pendingR[1] = 0.0;
      for (k = 0; k < 6; k++) {
        state.pendingPrediction[k] = 0.0;
      }
      state.pendingHost_s = -1.0;
      state.configured = true;
      memcpy(&state.options[0], &shadowOptions[0], 11U * sizeof(double));
      shadowDiagnostic->historyReset = true;
    }
    if (rtIsInf(t_sim_s) || rtIsNaN(t_sim_s) ||
        (rtIsInf(pollStart_s) || rtIsNaN(pollStart_s)) ||
        (rtIsInf(commandTime_sim_s) || rtIsNaN(commandTime_sim_s)) ||
        (rtIsInf(tipFrame) || rtIsNaN(tipFrame)) || (tipFrame < 0.0)) {
      state.segmentActive = false;
      state.pending = false;
      memset(&state.history[0], 0, 128U * sizeof(double));
      state.velocity[0] = 0.0;
      state.velocity[1] = 0.0;
      state.neuralStart_s = -1.0;
      state.segmentStart_s = -1.0;
      shadowDiagnostic->status = 1U;
      shadowDiagnostic->historyReset = true;

      /*  A repeated execution of the SAME sensor/simulation tick cannot score
       * or */
      /*  advance history a second time. A repeated frame at a NEW tick is a
       * dropout. */
    } else if ((t_sim_s == state.lastCall_sim_s) &&
               (pollStart_s == state.lastHost_s) &&
               (tipFrame == state.lastFrame)) {
      shadowDiagnostic->status = 2U;
      shadowDiagnostic->segment = state.segment;
    } else {
      double physicalDt;
      double simDt;
      boolean_T clockBack;
      boolean_T guard2;
      boolean_T sameFrame;
      simDt = t_sim_s - state.lastCall_sim_s;
      physicalDt = pollStart_s - state.lastHost_s;
      if (state.lastHost_s < 0.0) {
        physicalDt = rtNaN;
      }
      if ((state.lastCall_sim_s >= 0.0) &&
          ((simDt <= 0.0) || (physicalDt <= 0.0) ||
           (tipFrame < state.lastFrame))) {
        clockBack = true;
      } else {
        clockBack = false;
      }
      if ((state.lastFrame >= 0.0) && (tipFrame == state.lastFrame)) {
        sameFrame = true;
      } else {
        sameFrame = false;
      }
      state.lastCall_sim_s = t_sim_s;
      shadowDiagnostic->dtPhysical_s = physicalDt;
      shadowDiagnostic->commandAge_sim_s = t_sim_s - commandTime_sim_s;
      q[0] = lengths3_m[1];
      shadowDiagnostic->qMeasured_m[0] = lengths3_m[1];
      q[1] = lengths3_m[2];
      shadowDiagnostic->qMeasured_m[1] = lengths3_m[2];
      guard1 = false;
      guard2 = false;
      if (orientationValid && baseValid) {
        boolean_T b_pressure_bar[3];
        b_pressure_bar[0] =
            ((!rtIsInf(lengths3_m[0])) && (!rtIsNaN(lengths3_m[0])));
        b_pressure_bar[1] =
            ((!rtIsInf(lengths3_m[1])) && (!rtIsNaN(lengths3_m[1])));
        b_pressure_bar[2] =
            ((!rtIsInf(lengths3_m[2])) && (!rtIsNaN(lengths3_m[2])));
        if (all(b_pressure_bar)) {
          b_pressure_bar[0] =
              ((!rtIsInf(pressure_bar[0])) && (!rtIsNaN(pressure_bar[0])));
          b_pressure_bar[1] =
              ((!rtIsInf(pressure_bar[1])) && (!rtIsNaN(pressure_bar[1])));
          b_pressure_bar[2] =
              ((!rtIsInf(pressure_bar[2])) && (!rtIsNaN(pressure_bar[2])));
          if (all(b_pressure_bar)) {
            b_pressure_bar[0] =
                ((pressure_bar[0] >= -1.0E-9) &&
                 (pressure_bar[0] <= shadowOptions[10] + 1.0E-9));
            b_pressure_bar[1] =
                ((pressure_bar[1] >= -1.0E-9) &&
                 (pressure_bar[1] <= shadowOptions[10] + 1.0E-9));
            b_pressure_bar[2] =
                ((pressure_bar[2] >= -1.0E-9) &&
                 (pressure_bar[2] <= shadowOptions[10] + 1.0E-9));
            if (all(b_pressure_bar) &&
                (shadowDiagnostic->commandAge_sim_s >= -1.0E-9) &&
                (shadowDiagnostic->commandAge_sim_s <=
                 shadowOptions[7] + 1.0E-9)) {
              boolean_T y[2];
              y[0] = (fabs(lengths3_m[1]) <= shadowOptions[9]);
              y[1] = (fabs(lengths3_m[2]) <= shadowOptions[9]);
              if (b_all(y) && (!sameFrame)) {
                __m128d r;
                __m128d r1;
                boolean_T c_y[6];
                if (state.segmentActive && (!clockBack) && (physicalDt > 0.0) &&
                    (physicalDt <= shadowOptions[5] + 1.0E-9) &&
                    (simDt > 0.0) && (simDt <= shadowOptions[6] + 1.0E-9)) {
                  clockBack = true;
                  r = _mm_loadu_pd(&state.lastQ[0]);
                  r1 = _mm_loadu_pd(&state.velocity[0]);
                  _mm_storeu_pd(
                      &state.velocity[0],
                      _mm_add_pd(
                          r1, _mm_mul_pd(
                                  _mm_set1_pd(physicalDt /
                                              (shadowOptions[2] + physicalDt)),
                                  _mm_sub_pd(
                                      _mm_div_pd(
                                          _mm_sub_pd(
                                              _mm_loadu_pd(&lengths3_m[1]), r),
                                          _mm_set1_pd(physicalDt)),
                                      r1))));
                } else {
                  unsigned int qY;
                  clockBack = false;
                  state.pending = false;
                  memset(&state.history[0], 0, 128U * sizeof(double));
                  state.velocity[0] = 0.0;
                  state.velocity[1] = 0.0;
                  state.neuralStart_s = -1.0;
                  qY = state.segment + 1U;
                  if (state.segment + 1U < state.segment) {
                    qY = MAX_uint32_T;
                  }
                  state.segment = qY;
                  state.segmentActive = true;
                  state.segmentStart_s = pollStart_s;
                  shadowDiagnostic->historyReset = true;
                  shadowDiagnostic->status = 6U;
                }
                state.lastHost_s = pollStart_s;
                state.lastFrame = tipFrame;
                shadowDiagnostic->accepted = true;
                shadowDiagnostic->segment = state.segment;
                state.lastQ[0] = lengths3_m[1];
                shadowDiagnostic->velocity_mPerS[0] = state.velocity[0];
                state.lastQ[1] = lengths3_m[2];
                shadowDiagnostic->velocity_mPerS[1] = state.velocity[1];
                shadowDiagnostic->velocityWarm_s =
                    pollStart_s - state.segmentStart_s;
                /*  Score only a prior *issued* warm prediction on an adjacent
                 * accepted sample. */
                /*  Retimed replay uses only the elapsed interval, not the new
                 * q/v/p as inputs. */
                if (state.pending && clockBack) {
                  boolean_T guard3;
                  simDt = pollStart_s - state.pendingHost_s;
                  shadowDiagnostic->scoreOrigin_host_s = state.pendingHost_s;
                  shadowDiagnostic->horizonMismatch_s =
                      simDt - shadowOptions[1];
                  r = _mm_loadu_pd(&state.pendingP[0]);
                  _mm_storeu_pd(
                      &shadowDiagnostic->pressureChangeSinceOrigin_bar[0],
                      _mm_sub_pd(_mm_loadu_pd(&pressure_bar[0]), r));
                  shadowDiagnostic->pressureChangeSinceOrigin_bar[2] =
                      pressure_bar[2] - state.pendingP[2];
                  if (fabs(shadowDiagnostic->horizonMismatch_s) <=
                      shadowOptions[8] + 1.0E-9) {
                    for (k = 0; k < 6; k++) {
                      predictionAtMeasurement_m[k] = state.pendingPrediction[k];
                    }
                    r = _mm_loadu_pd(&state.pendingPrediction[0]);
                    r1 = _mm_loadu_pd(&q[0]);
                    _mm_storeu_pd(&shadowDiagnostic->errorAtMeasurement_m[0],
                                  _mm_sub_pd(r, r1));
                    r = _mm_loadu_pd(&state.pendingPrediction[2]);
                    r1 = _mm_loadu_pd(&q[0]);
                    _mm_storeu_pd(&shadowDiagnostic->errorAtMeasurement_m[2],
                                  _mm_sub_pd(r, r1));
                    r = _mm_loadu_pd(&state.pendingPrediction[4]);
                    r1 = _mm_loadu_pd(&q[0]);
                    _mm_storeu_pd(&shadowDiagnostic->errorAtMeasurement_m[4],
                                  _mm_sub_pd(r, r1));
                    shadowDiagnostic->scoreValid = true;
                  }
                  residualShadowFlow(state.pendingQ, state.pendingV,
                                     state.pendingP, state.pendingR, simDt,
                                     retimedPrediction_m);
                  for (k = 0; k < 6; k++) {
                    simDt = retimedPrediction_m[k];
                    c_y[k] = ((!rtIsInf(simDt)) && (!rtIsNaN(simDt)));
                  }
                  guard3 = false;
                  if (c_all(c_y)) {
                    for (k = 0; k < 6; k++) {
                      c_y[k] =
                          (fabs(retimedPrediction_m[k]) <= shadowOptions[9]);
                    }
                    if (c_all(c_y)) {
                      shadowDiagnostic->retimedScoreValid = true;
                      r = _mm_loadu_pd(&retimedPrediction_m[0]);
                      r1 = _mm_loadu_pd(&q[0]);
                      _mm_storeu_pd(&shadowDiagnostic->retimedError_m[0],
                                    _mm_sub_pd(r, r1));
                      r = _mm_loadu_pd(&retimedPrediction_m[2]);
                      r1 = _mm_loadu_pd(&q[0]);
                      _mm_storeu_pd(&shadowDiagnostic->retimedError_m[2],
                                    _mm_sub_pd(r, r1));
                      r = _mm_loadu_pd(&retimedPrediction_m[4]);
                      r1 = _mm_loadu_pd(&q[0]);
                      _mm_storeu_pd(&shadowDiagnostic->retimedError_m[4],
                                    _mm_sub_pd(r, r1));
                    } else {
                      guard3 = true;
                    }
                  } else {
                    guard3 = true;
                  }
                  if (guard3) {
                    shadowDiagnostic->retimedScoreValid = false;
                    for (k = 0; k < 6; k++) {
                      retimedPrediction_m[k] = rtNaN;
                    }
                  }
                }
                state.pending = false;
                if (shadowDiagnostic->velocityWarm_s <
                    shadowOptions[3] - 1.0E-9) {
                  if (shadowDiagnostic->status != 6) {
                    shadowDiagnostic->status = 3U;
                  }
                } else {
                  __m128d r2;
                  double u[384];
                  boolean_T bv[130];
                  if (state.neuralStart_s < 0.0) {
                    state.neuralStart_s = pollStart_s;
                  }
                  shadowDiagnostic->neuralWarm_s =
                      pollStart_s - state.neuralStart_s;
                  if (!clockBack) {
                    physicalDt = shadowOptions[1];
                  }
                  memcpy(&pre[0], &state.history[0], 128U * sizeof(double));
                  /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure
                   * candidate evaluation. */
                  /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds
                   * history fixed. */
                  r = _mm_loadu_pd(&network_mean[0]);
                  r1 = _mm_loadu_pd(&network_std[0]);
                  _mm_storeu_pd(
                      &z[0],
                      _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&lengths3_m[1]), r),
                                 r1));
                  r = _mm_loadu_pd(&state.velocity[0]);
                  r1 = _mm_loadu_pd(&network_mean[2]);
                  r2 = _mm_loadu_pd(&network_std[2]);
                  _mm_storeu_pd(&z[2], _mm_div_pd(_mm_sub_pd(r, r1), r2));
                  r = _mm_loadu_pd(&network_mean[4]);
                  r1 = _mm_loadu_pd(&network_std[4]);
                  _mm_storeu_pd(
                      &z[4],
                      _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&pressure_bar[0]), r),
                                 r1));
                  z[6] = (pressure_bar[2] - network_mean[6]) / network_std[6];
                  z[7] = (physicalDt - network_mean[7]) / network_std[7];
                  for (k = 0; k < 384; k++) {
                    simDt = 0.0;
                    for (i = 0; i < 8; i++) {
                      simDt += network_input[k + 384 * i] * z[i];
                    }
                    u[k] = simDt + w[k];
                  }
                  memset(&w[0], 0, 384U * sizeof(double));
                  for (k = 0; k < 128; k++) {
                    for (i = 0; i <= 382; i += 2) {
                      r = _mm_loadu_pd(&network_recurrent[i + 384 * k]);
                      r1 = _mm_loadu_pd(&w[i]);
                      _mm_storeu_pd(
                          &w[i],
                          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(
                                                           state.history[k]))));
                    }
                  }
                  memset(&q[0], 0, sizeof(double) << 1);
                  for (k = 0; k < 128; k++) {
                    simDt = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
                    simDt = simDt * state.history[k] +
                            (1.0 - simDt) *
                                tanh(u[k + 256] +
                                     0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0) *
                                         w[k + 256]);
                    g[k] = simDt;
                    r = _mm_loadu_pd(&network_output[k << 1]);
                    r1 = _mm_loadu_pd(&q[0]);
                    _mm_storeu_pd(
                        &q[0],
                        _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(simDt))));
                  }
                  memset(&b_y[0], 0, sizeof(double) << 1);
                  for (k = 0; k < 8; k++) {
                    r = _mm_loadu_pd(&network_skip[k << 1]);
                    r1 = _mm_loadu_pd(&b_y[0]);
                    _mm_storeu_pd(
                        &b_y[0],
                        _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(z[k]))));
                  }
                  q[0] =
                      80.0 *
                      tanh(10.0 * ((q[0] + b_y[0]) + network_outputBias_idx_0) /
                           80.0);
                  q[1] =
                      80.0 *
                      tanh(10.0 * ((q[1] + b_y[1]) + network_outputBias_idx_1) /
                           80.0);
                  memcpy(&x[0], &g[0], 128U * sizeof(double));
                  x[128] = q[0];
                  x[129] = q[1];
                  for (k = 0; k < 130; k++) {
                    simDt = x[k];
                    bv[k] = ((!rtIsInf(simDt)) && (!rtIsNaN(simDt)));
                  }
                  if (!d_all(bv)) {
                    state.segmentActive = false;
                    state.pending = false;
                    memset(&state.history[0], 0, 128U * sizeof(double));
                    state.velocity[0] = 0.0;
                    state.velocity[1] = 0.0;
                    state.neuralStart_s = -1.0;
                    state.segmentStart_s = -1.0;
                    shadowDiagnostic->status = 7U;
                    shadowDiagnostic->historyReset = true;
                  } else {
                    memcpy(&state.history[0], &g[0], 128U * sizeof(double));
                    shadowDiagnostic->advanced = true;
                    shadowDiagnostic->rawForce_N[0] = q[0];
                    shadowDiagnostic->rawForce_N[1] = q[1];
                    for (k = 0; k < 8; k++) {
                      network_mean[k] = fabs(z[k]);
                    }
                    shadowDiagnostic->maximumTrainingZ =
                        b_maximum(network_mean);
                    shadowDiagnostic->ready = (shadowDiagnostic->neuralWarm_s >=
                                               shadowOptions[4] - 1.0E-9);
                    if (!shadowDiagnostic->ready) {
                      shadowDiagnostic->status = 4U;
                    } else {
                      residualShadowFlow(&lengths3_m[1], state.velocity,
                                         pressure_bar, q, shadowOptions[1],
                                         predictionLength_m);
                      for (k = 0; k < 6; k++) {
                        simDt = predictionLength_m[k];
                        c_y[k] = ((!rtIsInf(simDt)) && (!rtIsNaN(simDt)));
                      }
                      if (!c_all(c_y)) {
                        guard2 = true;
                      } else {
                        for (k = 0; k < 6; k++) {
                          y_data[k] =
                              (fabs(predictionLength_m[k]) > shadowOptions[9]);
                        }
                        if (vectorAny(y_data, 6)) {
                          guard2 = true;
                        } else {
                          shadowDiagnostic->status = 5U;
                          shadowDiagnostic->predictionValid = true;
                          shadowDiagnostic->predictionTarget_host_s =
                              pollStart_s + shadowOptions[1];
                          state.pending = true;
                          state.pendingQ[0] = lengths3_m[1];
                          state.pendingV[0] = state.velocity[0];
                          state.pendingQ[1] = lengths3_m[2];
                          state.pendingV[1] = state.velocity[1];
                          state.pendingP[0] = pressure_bar[0];
                          state.pendingP[1] = pressure_bar[1];
                          state.pendingP[2] = pressure_bar[2];
                          state.pendingR[0] = q[0];
                          state.pendingR[1] = q[1];
                          for (k = 0; k < 6; k++) {
                            state.pendingPrediction[k] = predictionLength_m[k];
                          }
                          state.pendingHost_s = pollStart_s;
                          inverseContext->valid = true;
                          inverseContext->measured[0] = lengths3_m[1];
                          inverseContext->measured[2] = state.velocity[0];
                          inverseContext->measured[1] = lengths3_m[2];
                          inverseContext->measured[3] = state.velocity[1];
                          memcpy(&inverseContext->historyPre[0], &pre[0],
                                 128U * sizeof(double));
                          inverseContext->pressure_bar[0] = pressure_bar[0];
                          inverseContext->pressure_bar[1] = pressure_bar[1];
                          inverseContext->pressure_bar[2] = pressure_bar[2];
                          inverseContext->featureDt_s = physicalDt;
                          inverseContext->horizon_s = shadowOptions[1];
                          inverseContext->sampleTime_host_s = pollStart_s;
                          inverseContext->units = network_units;
                          shadowDiagnostic->inverseContextValid = true;
                        }
                      }
                    }
                  }
                }
              } else {
                guard1 = true;
              }
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
      if (guard2) {
        for (k = 0; k < 6; k++) {
          predictionLength_m[k] = rtNaN;
        }
        state.segmentActive = false;
        state.pending = false;
        memset(&state.history[0], 0, 128U * sizeof(double));
        state.velocity[0] = 0.0;
        state.velocity[1] = 0.0;
        state.neuralStart_s = -1.0;
        state.segmentStart_s = -1.0;
        shadowDiagnostic->ready = false;
        shadowDiagnostic->status = 7U;
        shadowDiagnostic->historyReset = true;
      }
      if (guard1) {
        state.segmentActive = false;
        state.pending = false;
        memset(&state.history[0], 0, 128U * sizeof(double));
        state.velocity[0] = 0.0;
        state.velocity[1] = 0.0;
        state.neuralStart_s = -1.0;
        state.segmentStart_s = -1.0;
        state.lastHost_s = pollStart_s;
        state.lastFrame = tipFrame;
        shadowDiagnostic->historyReset = true;
        if (sameFrame) {
          shadowDiagnostic->status = 2U;
        } else {
          shadowDiagnostic->status = 1U;
        }
        shadowDiagnostic->segment = state.segment;
      }
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void rlResidualShadowBlock_init(void)
{
  static const struct_T r = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, /* options */
      false,                                                   /* configured */
      false,      /* segmentActive */
      -1.0,       /* lastCall_sim_s */
      -1.0,       /* lastHost_s */
      -1.0,       /* lastFrame */
      -1.0,       /* segmentStart_s */
      {0.0, 0.0}, /* lastQ */
      {0.0, 0.0}, /* velocity */
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, /* history */
      -1.0,                           /* neuralStart_s */
      0U,                             /* segment */
      false,                          /* pending */
      {0.0, 0.0},                     /* pendingQ */
      {0.0, 0.0},                     /* pendingV */
      {0.0, 0.0, 0.0},                /* pendingP */
      {0.0, 0.0},                     /* pendingR */
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, /* pendingPrediction */
      -1.0                            /* pendingHost_s */
  };
  state = r;
}

/*
 * File trailer for rlResidualShadowBlock.c
 *
 * [EOF]
 */
