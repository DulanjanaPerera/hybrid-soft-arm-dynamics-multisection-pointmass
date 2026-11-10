/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualInverseShadowBlock.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:25
 */

/* Include Files */
#include "rlResidualInverseShadowBlock.h"
#include "all.h"
#include "any.h"
#include "norm.h"
#include "pacAllocatePressure.h"
#include "residualAllocatePressure.h"
#include "rlResidualInverseShadowBlock_data.h"
#include "rlResidualInverseShadowBlock_initialize.h"
#include "rlResidualInverseShadowBlock_internal_types.h"
#include "rlResidualInverseShadowBlock_rtwutil.h"
#include "rlResidualInverseShadowBlock_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * RLRESIDUALINVERSESHADOWBLOCK Stateless diagnostic inverse. NO valve
 * connection. Optional second block. History comes from forward observer; never
 * commits it. Targets the fixed trained nominal LHS at measured q/v and
 * supplied ddq. This is not PAC's adaptive/anchored control demand or a
 * controller replacement.
 *
 * Arguments    : const struct0_T *inverseContext
 *                const double desiredAcceleration_mPerS2[2]
 *                const double shadowPressureOptions[10]
 *                boolean_T variableCommon
 *                const double shadowInverseOptions[8]
 *                double residualPressurePrediction_bar[3]
 *                double nominalPressurePrediction_bar[3]
 *                struct1_T *inverseDiagnostic
 * Return Type  : void
 */
void rlResidualInverseShadowBlock(const struct0_T *inverseContext,
                                  const double desiredAcceleration_mPerS2[2],
                                  const double shadowPressureOptions[10],
                                  boolean_T variableCommon,
                                  const double shadowInverseOptions[8],
                                  double residualPressurePrediction_bar[3],
                                  double nominalPressurePrediction_bar[3],
                                  struct1_T *inverseDiagnostic)
{
  static const double b_dv[11] = {0.05,
                                  -0.0059523809523809521,
                                  0.00034722222222222224,
                                  -1.2776575276575277E-5,
                                  3.2856803690137022E-7,
                                  -6.26302709636043E-9,
                                  9.1445458654842411E-11,
                                  -1.024652859680243E-12,
                                  8.54742507061871E-15,
                                  -5.0026735363152986E-17,
                                  1.743355929321998E-19};
  static const double dv1[11] = {
      0.003968253968253968,    -0.00030864197530864197, 1.2325637325637326E-5,
      -3.2503504725726946E-7,  6.2431444706576982E-9,   -1.0013286157077007E-10,
      2.648557146298413E-12,   -8.04166027354068E-14,   1.5203929502675726E-15,
      -1.7701767897731058E-17, 1.3754287410824441E-19};
  static const double dv2[7] = {0.5,
                                -0.041666666666666664,
                                0.0013888888888888889,
                                -2.48015873015873E-5,
                                2.7557319223985888E-7,
                                -2.08767569878681E-9,
                                1.1470745597729725E-11};
  static double network_recurrent[49152];
  b_struct_T expl_temp;
  double network_input[3072];
  double network_bias[384];
  double w[384];
  double network_output[256];
  double network_skip[16];
  double P[9];
  double network_mean[8];
  double network_std[8];
  double z[8];
  double x[7];
  double b_P[4];
  double c_P[4];
  double l[3];
  double b_a[2];
  double d_P[2];
  double dv3[2];
  double dv4[2];
  double e_P[2];
  double target[2];
  int i;
  int k;
  unsigned char iterations;
  if (!isInitialized_rlResidualInverseShadowBlock) {
    rlResidualInverseShadowBlock_initialize();
  }
  residualPressurePrediction_bar[0] = rtNaN;
  nominalPressurePrediction_bar[0] = rtNaN;
  residualPressurePrediction_bar[1] = rtNaN;
  nominalPressurePrediction_bar[1] = rtNaN;
  residualPressurePrediction_bar[2] = rtNaN;
  nominalPressurePrediction_bar[2] = rtNaN;
  inverseDiagnostic->valid = false;
  inverseDiagnostic->converged = false;
  inverseDiagnostic->limited = false;
  inverseDiagnostic->nominalLimited = false;
  inverseDiagnostic->iterations = 0U;
  inverseDiagnostic->targetForce_N[0] = rtNaN;
  inverseDiagnostic->predictedForce_N[0] = rtNaN;
  inverseDiagnostic->inverseError_N[0] = rtNaN;
  inverseDiagnostic->nominalInverseError_N[0] = rtNaN;
  inverseDiagnostic->targetForce_N[1] = rtNaN;
  inverseDiagnostic->predictedForce_N[1] = rtNaN;
  inverseDiagnostic->inverseError_N[1] = rtNaN;
  inverseDiagnostic->nominalInverseError_N[1] = rtNaN;
  inverseDiagnostic->pressureDifference_bar[0] = rtNaN;
  inverseDiagnostic->pressureDifference_bar[1] = rtNaN;
  inverseDiagnostic->pressureDifference_bar[2] = rtNaN;
  inverseDiagnostic->sampleTime_host_s = rtNaN;
  if (inverseContext->valid) {
    boolean_T bv[2];
    bv[0] = ((!rtIsInf(desiredAcceleration_mPerS2[0])) &&
             (!rtIsNaN(desiredAcceleration_mPerS2[0])));
    bv[1] = ((!rtIsInf(desiredAcceleration_mPerS2[1])) &&
             (!rtIsNaN(desiredAcceleration_mPerS2[1])));
    if (all(bv)) {
      __m128d r;
      __m128d r1;
      __m128d r2;
      double u[384];
      double nominal[3];
      double network_outputBias[2];
      double a;
      double a__1;
      double a_tmp;
      double b;
      double bs;
      double s;
      double step;
      double w1;
      double w2;
      boolean_T bv1[7];
      boolean_T limited;
      boolean_T nomLimited;
      /* RESIDUALCONTROLLERNETWORK Fixed-size code-generation constants; run
       * packer first. */
      /*  Set units64/128 before a fresh run. Never replace weights during a
       * trial. */
      if (inverseContext->units == 64.0) {
        memcpy(&network_input[0], &n64.input[0], 3072U * sizeof(double));
        memcpy(&network_recurrent[0], &n64.recurrent[0],
               49152U * sizeof(double));
        memcpy(&network_bias[0], &n64.bias[0], 384U * sizeof(double));
        memcpy(&network_output[0], &n64.output[0], 256U * sizeof(double));
        memcpy(&network_skip[0], &n64.skip[0], 16U * sizeof(double));
        network_outputBias[0] = n64.outputBias[0];
        network_outputBias[1] = n64.outputBias[1];
        memcpy(&network_mean[0], &n64.mean[0], 8U * sizeof(double));
        memcpy(&network_std[0], &n64.std[0], 8U * sizeof(double));
      } else {
        memcpy(&network_input[0], &n128.input[0], 3072U * sizeof(double));
        memcpy(&network_recurrent[0], &n128.recurrent[0],
               49152U * sizeof(double));
        memcpy(&network_bias[0], &n128.bias[0], 384U * sizeof(double));
        memcpy(&network_output[0], &n128.output[0], 256U * sizeof(double));
        memcpy(&network_skip[0], &n128.skip[0], 16U * sizeof(double));
        network_outputBias[0] = n128.outputBias[0];
        network_outputBias[1] = n128.outputBias[1];
        memcpy(&network_mean[0], &n128.mean[0], 8U * sizeof(double));
        memcpy(&network_std[0], &n128.std[0], 8U * sizeof(double));
      }
      /*  Shared configuration validation only; no supervisor/anchor state is
       * run. */
      /* RESIDUALPRESSURECONFIG Same vector contract as PAC tuning, no feedback
       * gains. */
      /* FORCERESIDUALPHYSICS Batched differentiable EXACT nominal M,C,G
       * polynomials. */
      /*  q,v:2xbatch, numeric or unformatted dlarray. P rows M11,M12,M22, */
      /*  C11,C12,C21,C22,G1,G2. Fixed nominal basis, no fitted approximation.
       */
      /*  Generated rational coefficients; see generate_force_polynomials.py. */
      w1 = 7889.5463510848131 * inverseContext->measured[0] +
           3944.7731755424065 * inverseContext->measured[1];
      w2 = 3944.7731755424065 * inverseContext->measured[0] +
           7889.5463510848131 * inverseContext->measured[1];
      s = inverseContext->measured[0] * w1 + inverseContext->measured[1] * w2;
      a_tmp = 0.0 * s;
      a = a_tmp + 1.743355929321998E-19;
      step = a_tmp;
      b = a_tmp + 1.3754287410824441E-19;
      bs = a_tmp;
      for (k = 0; k < 10; k++) {
        step = step * s + a;
        a = a * s + b_dv[9 - k];
        bs = bs * s + b;
        b = b * s + dv1[9 - k];
      }
      a__1 = a_tmp + 1.1470745597729725E-11;
      for (k = 0; k < 6; k++) {
        a_tmp = a_tmp * s + a__1;
        a__1 = a__1 * s + dv2[5 - k];
      }
      double a2;
      double a2_tmp;
      double b_a2_tmp;
      double b_c1_tmp;
      double c1_tmp;
      double c_a2_tmp;
      double c_c1_tmp;
      double d_a2_tmp;
      double e_a2_tmp;
      step *= 2.0;
      a2_tmp = 2.0 * bs;
      s = 2.0 * b;
      b_a2_tmp = w1 * w1;
      c_a2_tmp = s * w1;
      d_a2_tmp = step * w2;
      e_a2_tmp = d_a2_tmp * 7889.5463510848131;
      a2 = 0.0030314292099999996 * ((e_a2_tmp + a2_tmp * w2 * b_a2_tmp) +
                                    c_a2_tmp * 3944.7731755424065);
      c1_tmp = w2 * w2;
      bs = step * w1;
      b_c1_tmp = bs * 7889.5463510848131;
      c_c1_tmp = a2_tmp * w1 * c1_tmp;
      step = s * w2;
      s = 0.0030314292099999996 *
          ((b_c1_tmp + c_c1_tmp) + step * 3944.7731755424065);
      a__1 = a * 7889.5463510848131;
      P[0] = 0.0030314292099999996 * (a__1 + b * b_a2_tmp);
      P[1] = 0.0030314292099999996 * (a * 3944.7731755424065 + b * w1 * w2);
      P[2] = 0.0030314292099999996 * (a__1 + b * c1_tmp);
      P[3] = 0.5 * (0.0030314292099999996 *
                        ((b_c1_tmp + a2_tmp * rt_powd_snf(w1, 3.0)) +
                         c_a2_tmp * 7889.5463510848131) *
                        inverseContext->measured[2] +
                    a2 * inverseContext->measured[3]);
      P[4] =
          0.5 *
          (a2 * inverseContext->measured[2] +
           (2.0 * (0.0030314292099999996 *
                   ((d_a2_tmp * 3944.7731755424065 + c_c1_tmp) +
                    b * (3944.7731755424065 * w2 + 7889.5463510848131 * w1))) -
            s) *
               inverseContext->measured[3]);
      P[5] =
          0.5 *
          ((2.0 * (0.0030314292099999996 *
                   ((bs * 3944.7731755424065 + a2_tmp * b_a2_tmp * w2) +
                    b * (7889.5463510848131 * w2 + 3944.7731755424065 * w1))) -
            a2) *
               inverseContext->measured[2] +
           s * inverseContext->measured[3]);
      P[6] = 0.5 * (s * inverseContext->measured[2] +
                    0.0030314292099999996 *
                        ((e_a2_tmp + a2_tmp * rt_powd_snf(w2, 3.0)) +
                         step * 7889.5463510848131) *
                        inverseContext->measured[3]);
      s = -0.34160382 * a_tmp;
      P[7] = s * w1;
      P[8] = s * w2;
      /* ARMS_SINGLE_BOUNDFORCE Soft actuator-length restoring force in
       * q=[l2;l3]. */
      /*  Default mode1: all l=[-l2-l3;l2;l3] use the same bounds, in metres. */
      /*  Mode0 reproduces historical coordinate-only l2/l3 bounds for old
       * analyses. */
      /*  coefficient*q equals force_N; coefficient is NOT the tangent
       * stiffness. */
      l[0] = -inverseContext->measured[0] - inverseContext->measured[1];
      s = 500000.0 *
          ((tanh(2000.0 * (l[0] - 0.035)) + 2.0) -
           tanh(2000.0 * (l[0] - -0.035))) *
          l[0];
      /*  J=dl/dq=[-1 -1;1 0;0 1]; generalized force is J'*actuatorForce. */
      b_P[0] = P[0];
      b_P[2] = P[1];
      b_P[1] = P[1];
      b_P[3] = P[2];
      c_P[0] = P[3] + 40.0;
      c_P[2] = P[4];
      c_P[1] = P[5];
      c_P[3] = P[6] + 40.0;
      memset(&d_P[0], 0, sizeof(double) << 1);
      memset(&e_P[0], 0, sizeof(double) << 1);
      r = _mm_loadu_pd(&b_P[0]);
      r1 = _mm_loadu_pd(&d_P[0]);
      _mm_storeu_pd(
          &d_P[0],
          _mm_add_pd(
              r1, _mm_mul_pd(r, _mm_set1_pd(desiredAcceleration_mPerS2[0]))));
      r = _mm_loadu_pd(&c_P[0]);
      r1 = _mm_loadu_pd(&e_P[0]);
      _mm_storeu_pd(
          &e_P[0],
          _mm_add_pd(r1,
                     _mm_mul_pd(r, _mm_set1_pd(inverseContext->measured[2]))));
      r = _mm_loadu_pd(&b_P[2]);
      r1 = _mm_loadu_pd(&d_P[0]);
      _mm_storeu_pd(
          &d_P[0],
          _mm_add_pd(
              r1, _mm_mul_pd(r, _mm_set1_pd(desiredAcceleration_mPerS2[1]))));
      r = _mm_loadu_pd(&c_P[2]);
      r1 = _mm_loadu_pd(&e_P[0]);
      _mm_storeu_pd(
          &e_P[0],
          _mm_add_pd(r1,
                     _mm_mul_pd(r, _mm_set1_pd(inverseContext->measured[3]))));
      target[0] = (d_P[0] + e_P[0]) + (1300.0 * inverseContext->measured[0] +
                                       650.0 * inverseContext->measured[1]);
      target[1] = (d_P[1] + e_P[1]) + (650.0 * inverseContext->measured[0] +
                                       1300.0 * inverseContext->measured[1]);
      d_P[0] =
          500000.0 *
              ((tanh(2000.0 * (inverseContext->measured[0] - 0.035)) + 2.0) -
               tanh(2000.0 * (inverseContext->measured[0] - -0.035))) *
              inverseContext->measured[0] -
          s;
      d_P[1] =
          500000.0 *
              ((tanh(2000.0 * (inverseContext->measured[1] - 0.035)) + 2.0) -
               tanh(2000.0 * (inverseContext->measured[1] - -0.035))) *
              inverseContext->measured[1] -
          s;
      step = shadowPressureOptions[7] * inverseContext->horizon_s;
      r = _mm_loadu_pd(&target[0]);
      r1 = _mm_loadu_pd(&d_P[0]);
      r2 = _mm_loadu_pd(&P[7]);
      _mm_storeu_pd(&target[0], _mm_add_pd(_mm_add_pd(r, r1), r2));
      _mm_storeu_pd(&d_P[0], _mm_set1_pd(0.0));
      nomLimited = pacAllocatePressure(
          target, d_P, shadowPressureOptions[0], shadowPressureOptions[1],
          shadowPressureOptions[2], shadowPressureOptions[3],
          shadowPressureOptions[4], shadowPressureOptions[5],
          shadowPressureOptions[6], shadowPressureOptions[8],
          shadowPressureOptions[9], variableCommon, nominal, b_a, e_P);
      r = _mm_loadu_pd(&inverseContext->pressure_bar[0]);
      r1 = _mm_set1_pd(step);
      _mm_storeu_pd(&dv3[0], _mm_sub_pd(r, r1));
      _mm_storeu_pd(&dv4[0], _mm_add_pd(r, r1));
      b_a[0] = fmin(fmax(nominal[0], dv3[0]), dv4[0]);
      b_a[1] = fmin(fmax(nominal[1], dv3[1]), dv4[1]);
      r = _mm_loadu_pd(&b_a[0]);
      _mm_storeu_pd(&nominalPressurePrediction_bar[0], r);
      nominalPressurePrediction_bar[2] =
          fmin(fmax(nominal[2], inverseContext->pressure_bar[2] - step),
               inverseContext->pressure_bar[2] + step);
      /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
       * evaluation. */
      /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history fixed.
       */
      r = _mm_loadu_pd(&network_mean[0]);
      r1 = _mm_loadu_pd(&network_std[0]);
      _mm_storeu_pd(
          &z[0],
          _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&inverseContext->measured[0]), r),
                     r1));
      r = _mm_loadu_pd(&network_mean[2]);
      r1 = _mm_loadu_pd(&network_std[2]);
      _mm_storeu_pd(
          &z[2],
          _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&inverseContext->measured[2]), r),
                     r1));
      r = _mm_loadu_pd(&network_mean[4]);
      r1 = _mm_loadu_pd(&network_std[4]);
      _mm_storeu_pd(
          &z[4],
          _mm_div_pd(
              _mm_sub_pd(_mm_loadu_pd(&inverseContext->pressure_bar[0]), r),
              r1));
      z[6] =
          (inverseContext->pressure_bar[2] - network_mean[6]) / network_std[6];
      z[7] = (inverseContext->featureDt_s - network_mean[7]) / network_std[7];
      for (k = 0; k < 384; k++) {
        s = 0.0;
        for (i = 0; i < 8; i++) {
          s += network_input[k + 384 * i] * z[i];
        }
        u[k] = s + network_bias[k];
      }
      memset(&w[0], 0, 384U * sizeof(double));
      for (k = 0; k < 128; k++) {
        for (i = 0; i <= 382; i += 2) {
          r = _mm_loadu_pd(&network_recurrent[i + 384 * k]);
          r1 = _mm_loadu_pd(&w[i]);
          _mm_storeu_pd(
              &w[i],
              _mm_add_pd(
                  r1,
                  _mm_mul_pd(r, _mm_set1_pd(inverseContext->historyPre[k]))));
        }
      }
      memset(&b_a[0], 0, sizeof(double) << 1);
      for (k = 0; k < 128; k++) {
        s = 0.5 * (tanh(0.5 * (u[k + 128] + w[k + 128])) + 1.0);
        s = s * inverseContext->historyPre[k] +
            (1.0 - s) *
                tanh(u[k + 256] +
                     0.5 * (tanh(0.5 * (u[k] + w[k])) + 1.0) * w[k + 256]);
        r = _mm_loadu_pd(&network_output[k << 1]);
        r1 = _mm_loadu_pd(&b_a[0]);
        _mm_storeu_pd(&b_a[0], _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(s))));
      }
      memset(&d_P[0], 0, sizeof(double) << 1);
      for (k = 0; k < 8; k++) {
        r = _mm_loadu_pd(&network_skip[k << 1]);
        r1 = _mm_loadu_pd(&d_P[0]);
        _mm_storeu_pd(&d_P[0],
                      _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(z[k]))));
      }
      expl_temp.enableCommonAdaptation = variableCommon;
      expl_temp.forcePerBar_N = shadowPressureOptions[9];
      expl_temp.activation_bar = shadowPressureOptions[8];
      expl_temp.maxDeltaA_bar = shadowPressureOptions[6];
      expl_temp.maxDeltaS_bar = shadowPressureOptions[5];
      expl_temp.maximum_bar = shadowPressureOptions[4];
      expl_temp.minimum_bar = shadowPressureOptions[3];
      expl_temp.aTrim_bar = shadowPressureOptions[2];
      expl_temp.sTrim_bar = shadowPressureOptions[1];
      expl_temp.common_bar = shadowPressureOptions[0];
      b_a[0] = target[0] -
               80.0 * tanh(10.0 * ((b_a[0] + d_P[0]) + network_outputBias[0]) /
                           80.0);
      b_a[1] = target[1] -
               80.0 * tanh(10.0 * ((b_a[1] + d_P[1]) + network_outputBias[1]) /
                           80.0);
      limited = residualAllocatePressure(
          target, inverseContext->pressure_bar, b_a, target, step,
          inverseContext->measured, inverseContext->featureDt_s,
          inverseContext->historyPre, network_input, network_recurrent,
          network_bias, network_output, network_skip, network_outputBias,
          network_mean, network_std, &expl_temp, shadowInverseOptions, l, e_P,
          &iterations);
      x[0] = l[0];
      x[1] = l[1];
      x[2] = l[2];
      x[3] = e_P[0];
      x[5] = target[0];
      x[4] = e_P[1];
      x[6] = target[1];
      for (k = 0; k < 7; k++) {
        s = x[k];
        bv1[k] = ((!rtIsInf(s)) && (!rtIsNaN(s)));
      }
      if (d_all(bv1)) {
        residualPressurePrediction_bar[0] = l[0];
        residualPressurePrediction_bar[1] = l[1];
        residualPressurePrediction_bar[2] = l[2];
        inverseDiagnostic->valid = true;
        inverseDiagnostic->converged = (b_norm(e_P) <= shadowInverseOptions[2]);
        inverseDiagnostic->limited = limited;
        if (nomLimited) {
          inverseDiagnostic->nominalLimited = true;
        } else {
          boolean_T y[3];
          y[0] = (fabs(nominal[0] - nominalPressurePrediction_bar[0]) > 1.0E-9);
          y[1] = (fabs(nominal[1] - nominalPressurePrediction_bar[1]) > 1.0E-9);
          y[2] = (fabs(nominal[2] - nominalPressurePrediction_bar[2]) > 1.0E-9);
          if (any(y)) {
            inverseDiagnostic->nominalLimited = true;
          } else {
            inverseDiagnostic->nominalLimited = false;
          }
        }
        inverseDiagnostic->iterations = iterations;
        s = nominalPressurePrediction_bar[0];
        step = nominalPressurePrediction_bar[1];
        a__1 = nominalPressurePrediction_bar[2];
        r2 = _mm_loadu_pd(&target[0]);
        _mm_storeu_pd(&inverseDiagnostic->targetForce_N[0], r2);
        r = _mm_loadu_pd(&e_P[0]);
        r1 = _mm_add_pd(r2, r);
        _mm_storeu_pd(&inverseDiagnostic->predictedForce_N[0], r1);
        _mm_storeu_pd(&inverseDiagnostic->inverseError_N[0], r);
        r = _mm_loadu_pd(&dv[2]);
        r = _mm_mul_pd(r, _mm_set1_pd(step));
        r1 = _mm_add_pd(_mm_set1_pd(-13.273228961416876 * s), r);
        r = _mm_loadu_pd(&dv[4]);
        r = _mm_mul_pd(r, _mm_set1_pd(a__1));
        r = _mm_add_pd(r1, r);
        r = _mm_sub_pd(r, r2);
        _mm_storeu_pd(&inverseDiagnostic->nominalInverseError_N[0], r);
        r = _mm_loadu_pd(&l[0]);
        r1 = _mm_loadu_pd(&nominalPressurePrediction_bar[0]);
        _mm_storeu_pd(&inverseDiagnostic->pressureDifference_bar[0],
                      _mm_sub_pd(r, r1));
        inverseDiagnostic->pressureDifference_bar[2] =
            l[2] - nominalPressurePrediction_bar[2];
        inverseDiagnostic->sampleTime_host_s =
            inverseContext->sampleTime_host_s;
      }
    }
  }
}

/*
 * File trailer for rlResidualInverseShadowBlock.c
 *
 * [EOF]
 */
