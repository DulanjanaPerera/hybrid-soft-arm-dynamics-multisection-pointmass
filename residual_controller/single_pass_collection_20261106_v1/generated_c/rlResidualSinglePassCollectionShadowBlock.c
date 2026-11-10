/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualSinglePassCollectionShadowBlock.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 06-Nov-2026 18:50:26
 */

/* Include Files */
#include "rlResidualSinglePassCollectionShadowBlock.h"
#include "all.h"
#include "any.h"
#include "minOrMax.h"
#include "norm.h"
#include "rlResidualSinglePassCollectionShadowBlock_data.h"
#include "rlResidualSinglePassCollectionShadowBlock_initialize.h"
#include "rlResidualSinglePassCollectionShadowBlock_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    double d;
    y = fabs(u0);
    d = fabs(u1);
    if (rtIsInf(u1)) {
      if (y == 1.0) {
        y = 1.0;
      } else if (y > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d == 0.0) {
      y = 1.0;
    } else if (d == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

/*
 * RLRESIDUALSINGLEPASSCOLLECTIONSHADOWBLOCK Paste into discrete0.05s function.
 *  Observation ONLY: outputs to logging/Scopes, never valves/control/Kill.
 *  Same5inputs/3outputs/bus as rlResidualSinglePassShadowBlock.
 *  variableCommon is a retained compatibility input; NOT used by this policy.
 *  Pressure options used:4minimum(0),5maximum,8slew,9activation(0),10force/bar.
 *  Common pressure/trim/s-a caps(1:3,6:7) are not part of collection
 * allocation. One GRU evaluation at actual pressure; no inverse
 * iterations/state commits.
 *
 * Arguments    : const struct0_T *inverseContext
 *                const double desiredAcceleration_mPerS2[2]
 *                const double shadowPressureOptions[10]
 *                boolean_T variableCommon
 *                const double singlePassOptions[2]
 *                double c_singlePassPressurePrediction_[3]
 *                double nominalPressurePrediction_bar[3]
 *                struct1_T *singlePassDiagnostic
 * Return Type  : void
 */
void rlResidualSinglePassCollectionShadowBlock(
    const struct0_T *inverseContext, const double desiredAcceleration_mPerS2[2],
    const double shadowPressureOptions[10], boolean_T variableCommon,
    const double singlePassOptions[2],
    double c_singlePassPressurePrediction_[3],
    double nominalPressurePrediction_bar[3], struct1_T *singlePassDiagnostic)
{
  static const double dv[11] = {0.05,
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
  static const signed char b_b[6] = {-1, -1, 1, 0, 0, 1};
  double network_input[3072];
  double w[384];
  double network_output[256];
  double x[140];
  double b_x[16];
  double network_skip[16];
  double P[9];
  double network_mean[8];
  double network_std[8];
  double z[8];
  double B[6];
  double b_P[4];
  double c_P[4];
  double corrected[3];
  double p[3];
  double pn[3];
  double y[3];
  double b_a[2];
  double dv3[2];
  double dv4[2];
  double forceBar[2];
  double request[2];
  double target[2];
  int i;
  int k;
  (void)variableCommon;
  if (!isInitialized_rlResidualSinglePassCollectionShadowBlock) {
    rlResidualSinglePassCollectionShadowBlock_initialize();
  }
  c_singlePassPressurePrediction_[0] = rtNaN;
  nominalPressurePrediction_bar[0] = rtNaN;
  c_singlePassPressurePrediction_[1] = rtNaN;
  nominalPressurePrediction_bar[1] = rtNaN;
  c_singlePassPressurePrediction_[2] = rtNaN;
  nominalPressurePrediction_bar[2] = rtNaN;
  singlePassDiagnostic->valid = false;
  singlePassDiagnostic->networkEvaluations = 0U;
  singlePassDiagnostic->allocationLimited = false;
  singlePassDiagnostic->slewLimited = false;
  singlePassDiagnostic->nominalAllocationLimited = false;
  singlePassDiagnostic->nominalSlewLimited = false;
  singlePassDiagnostic->frozenBalanceSatisfied = false;
  singlePassDiagnostic->targetForce_N[0] = rtNaN;
  singlePassDiagnostic->residualForce_N[0] = rtNaN;
  singlePassDiagnostic->correctedForceRequest_N[0] = rtNaN;
  singlePassDiagnostic->frozenTotalForce_N[0] = rtNaN;
  singlePassDiagnostic->allocationError_N[0] = rtNaN;
  singlePassDiagnostic->slewError_N[0] = rtNaN;
  singlePassDiagnostic->frozenBalanceError_N[0] = rtNaN;
  singlePassDiagnostic->nominalBalanceError_N[0] = rtNaN;
  singlePassDiagnostic->targetForce_N[1] = rtNaN;
  singlePassDiagnostic->residualForce_N[1] = rtNaN;
  singlePassDiagnostic->correctedForceRequest_N[1] = rtNaN;
  singlePassDiagnostic->frozenTotalForce_N[1] = rtNaN;
  singlePassDiagnostic->allocationError_N[1] = rtNaN;
  singlePassDiagnostic->slewError_N[1] = rtNaN;
  singlePassDiagnostic->frozenBalanceError_N[1] = rtNaN;
  singlePassDiagnostic->nominalBalanceError_N[1] = rtNaN;
  singlePassDiagnostic->residualInputPressure_bar[0] = rtNaN;
  singlePassDiagnostic->correctedPressureBeforeSlew_bar[0] = rtNaN;
  singlePassDiagnostic->nominalPressureBeforeSlew_bar[0] = rtNaN;
  singlePassDiagnostic->c_pressureDifferenceBeforeSlew_[0] = rtNaN;
  singlePassDiagnostic->pressureDifference_bar[0] = rtNaN;
  singlePassDiagnostic->residualInputPressure_bar[1] = rtNaN;
  singlePassDiagnostic->correctedPressureBeforeSlew_bar[1] = rtNaN;
  singlePassDiagnostic->nominalPressureBeforeSlew_bar[1] = rtNaN;
  singlePassDiagnostic->c_pressureDifferenceBeforeSlew_[1] = rtNaN;
  singlePassDiagnostic->pressureDifference_bar[1] = rtNaN;
  singlePassDiagnostic->residualInputPressure_bar[2] = rtNaN;
  singlePassDiagnostic->correctedPressureBeforeSlew_bar[2] = rtNaN;
  singlePassDiagnostic->nominalPressureBeforeSlew_bar[2] = rtNaN;
  singlePassDiagnostic->c_pressureDifferenceBeforeSlew_[2] = rtNaN;
  singlePassDiagnostic->pressureDifference_bar[2] = rtNaN;
  singlePassDiagnostic->sampleTime_host_s = rtNaN;
  singlePassDiagnostic->horizon_s = rtNaN;
  if (inverseContext->valid) {
    double s;
    boolean_T bv[140];
    x[0] = inverseContext->measured[0];
    x[1] = inverseContext->measured[1];
    x[2] = inverseContext->measured[2];
    x[3] = inverseContext->measured[3];
    x[4] = inverseContext->pressure_bar[0];
    x[5] = inverseContext->pressure_bar[1];
    x[6] = inverseContext->pressure_bar[2];
    memcpy(&x[7], &inverseContext->historyPre[0], 128U * sizeof(double));
    x[135] = inverseContext->featureDt_s;
    x[136] = inverseContext->horizon_s;
    x[137] = inverseContext->sampleTime_host_s;
    x[138] = desiredAcceleration_mPerS2[0];
    x[139] = desiredAcceleration_mPerS2[1];
    for (k = 0; k < 140; k++) {
      s = x[k];
      bv[k] = ((!rtIsInf(s)) && (!rtIsNaN(s)));
    }
    if (all(bv) && (!(inverseContext->featureDt_s <= 0.0)) &&
        (!(inverseContext->horizon_s <= 0.0))) {
      double network_outputBias_idx_0;
      double network_outputBias_idx_1;
      boolean_T b_inverseContext[3];
      /* RESIDUALCONTROLLERNETWORK Fixed-size code-generation constants; run
       * packer first. */
      /*  Set units64/128 before a fresh run. Never replace weights during a
       * trial. */
      if (singlePassOptions[0] == 64.0) {
        memcpy(&network_input[0], &n64.input[0], 3072U * sizeof(double));
        memcpy(&network_recurrent[0], &n64.recurrent[0],
               49152U * sizeof(double));
        memcpy(&w[0], &n64.bias[0], 384U * sizeof(double));
        memcpy(&network_output[0], &n64.output[0], 256U * sizeof(double));
        memcpy(&network_skip[0], &n64.skip[0], 16U * sizeof(double));
        network_outputBias_idx_0 = n64.outputBias[0];
        network_outputBias_idx_1 = n64.outputBias[1];
        memcpy(&network_mean[0], &n64.mean[0], 8U * sizeof(double));
        memcpy(&network_std[0], &n64.std[0], 8U * sizeof(double));
      } else {
        memcpy(&network_input[0], &n128.input[0], 3072U * sizeof(double));
        memcpy(&network_recurrent[0], &n128.recurrent[0],
               49152U * sizeof(double));
        memcpy(&w[0], &n128.bias[0], 384U * sizeof(double));
        memcpy(&network_output[0], &n128.output[0], 256U * sizeof(double));
        memcpy(&network_skip[0], &n128.skip[0], 16U * sizeof(double));
        network_outputBias_idx_0 = n128.outputBias[0];
        network_outputBias_idx_1 = n128.outputBias[1];
        memcpy(&network_mean[0], &n128.mean[0], 8U * sizeof(double));
        memcpy(&network_std[0], &n128.std[0], 8U * sizeof(double));
      }
      b_inverseContext[0] = (inverseContext->pressure_bar[0] < -1.0E-9);
      b_inverseContext[1] = (inverseContext->pressure_bar[1] < -1.0E-9);
      b_inverseContext[2] = (inverseContext->pressure_bar[2] < -1.0E-9);
      if (!any(b_inverseContext)) {
        b_inverseContext[0] = (inverseContext->pressure_bar[0] >
                               shadowPressureOptions[4] + 1.0E-9);
        b_inverseContext[1] = (inverseContext->pressure_bar[1] >
                               shadowPressureOptions[4] + 1.0E-9);
        b_inverseContext[2] = (inverseContext->pressure_bar[2] >
                               shadowPressureOptions[4] + 1.0E-9);
        if (!any(b_inverseContext)) {
          __m128d r;
          __m128d r1;
          __m128d r2;
          double u[384];
          double P_tmp;
          double a;
          double a2_tmp;
          double a_tmp;
          double as;
          double b;
          double bs;
          double c1;
          double c1_tmp_tmp;
          double l_idx_0;
          double l_idx_1;
          double w1;
          double w2;
          boolean_T bv1[16];
          boolean_T allocationLimited;
          boolean_T nominalLimited;
          /* FORCERESIDUALPHYSICS Batched differentiable EXACT nominal M,C,G
           * polynomials. */
          /*  q,v:2xbatch, numeric or unformatted dlarray. P rows M11,M12,M22,
           */
          /*  C11,C12,C21,C22,G1,G2. Fixed nominal basis, no fitted
           * approximation. */
          /*  Generated rational coefficients; see
           * generate_force_polynomials.py. */
          w1 = 7889.5463510848131 * inverseContext->measured[0] +
               3944.7731755424065 * inverseContext->measured[1];
          w2 = 3944.7731755424065 * inverseContext->measured[0] +
               7889.5463510848131 * inverseContext->measured[1];
          s = inverseContext->measured[0] * w1 +
              inverseContext->measured[1] * w2;
          a_tmp = 0.0 * s;
          a = a_tmp + 1.743355929321998E-19;
          as = a_tmp;
          b = a_tmp + 1.3754287410824441E-19;
          bs = a_tmp;
          for (k = 0; k < 10; k++) {
            as = as * s + a;
            a = a * s + dv[9 - k];
            bs = bs * s + b;
            b = b * s + dv1[9 - k];
          }
          l_idx_1 = a_tmp + 1.1470745597729725E-11;
          for (k = 0; k < 6; k++) {
            a_tmp = a_tmp * s + l_idx_1;
            l_idx_1 = l_idx_1 * s + dv2[5 - k];
          }
          double a2;
          double b_a2_tmp;
          double c1_tmp;
          double c_a2_tmp;
          double d_a2_tmp;
          l_idx_1 = 2.0 * as;
          l_idx_0 = 2.0 * bs;
          a2_tmp = 2.0 * b;
          bs = w1 * w1;
          b_a2_tmp = a2_tmp * w1;
          c_a2_tmp = l_idx_1 * w2;
          d_a2_tmp = c_a2_tmp * 7889.5463510848131;
          a2 = 0.0030314292099999996 *
               ((d_a2_tmp + l_idx_0 * w2 * bs) + b_a2_tmp * 3944.7731755424065);
          c1_tmp = w2 * w2;
          c1_tmp_tmp = l_idx_1 * w1;
          s = c1_tmp_tmp * 7889.5463510848131;
          as = l_idx_0 * w1 * c1_tmp;
          l_idx_1 = a2_tmp * w2;
          c1 =
              0.0030314292099999996 * ((s + as) + l_idx_1 * 3944.7731755424065);
          P_tmp = a * 7889.5463510848131;
          P[0] = 0.0030314292099999996 * (P_tmp + b * bs);
          P[1] = 0.0030314292099999996 * (a * 3944.7731755424065 + b * w1 * w2);
          P[2] = 0.0030314292099999996 * (P_tmp + b * c1_tmp);
          P[3] = 0.5 * (0.0030314292099999996 *
                            ((s + l_idx_0 * rt_powd_snf(w1, 3.0)) +
                             b_a2_tmp * 7889.5463510848131) *
                            inverseContext->measured[2] +
                        a2 * inverseContext->measured[3]);
          P[4] = 0.5 * (a2 * inverseContext->measured[2] +
                        (2.0 * (0.0030314292099999996 *
                                ((c_a2_tmp * 3944.7731755424065 + as) +
                                 b * (3944.7731755424065 * w2 +
                                      7889.5463510848131 * w1))) -
                         c1) *
                            inverseContext->measured[3]);
          P[5] =
              0.5 *
              ((2.0 * (0.0030314292099999996 *
                       ((c1_tmp_tmp * 3944.7731755424065 + l_idx_0 * bs * w2) +
                        b * (7889.5463510848131 * w2 +
                             3944.7731755424065 * w1))) -
                a2) *
                   inverseContext->measured[2] +
               c1 * inverseContext->measured[3]);
          P[6] = 0.5 * (c1 * inverseContext->measured[2] +
                        0.0030314292099999996 *
                            ((d_a2_tmp + l_idx_0 * rt_powd_snf(w2, 3.0)) +
                             l_idx_1 * 7889.5463510848131) *
                            inverseContext->measured[3]);
          s = -0.34160382 * a_tmp;
          P[7] = s * w1;
          P[8] = s * w2;
          /* ARMS_SINGLE_BOUNDFORCE Soft actuator-length restoring force in
           * q=[l2;l3]. */
          /*  Default mode1: all l=[-l2-l3;l2;l3] use the same bounds, in
           * metres. */
          /*  Mode0 reproduces historical coordinate-only l2/l3 bounds for old
           * analyses. */
          /*  coefficient*q equals force_N; coefficient is NOT the tangent
           * stiffness. */
          l_idx_0 = -inverseContext->measured[0] - inverseContext->measured[1];
          s = 500000.0 *
              ((tanh(2000.0 * (l_idx_0 - 0.035)) + 2.0) -
               tanh(2000.0 * (l_idx_0 - -0.035))) *
              l_idx_0;
          /*  J=dl/dq=[-1 -1;1 0;0 1]; generalized force is J'*actuatorForce. */
          b_P[0] = P[0];
          b_P[2] = P[1];
          b_P[1] = P[1];
          b_P[3] = P[2];
          c_P[0] = P[3] + 40.0;
          c_P[2] = P[4];
          c_P[1] = P[5];
          c_P[3] = P[6] + 40.0;
          memset(&forceBar[0], 0, sizeof(double) << 1);
          memset(&b_a[0], 0, sizeof(double) << 1);
          r = _mm_loadu_pd(&b_P[0]);
          r1 = _mm_loadu_pd(&forceBar[0]);
          _mm_storeu_pd(
              &forceBar[0],
              _mm_add_pd(
                  r1,
                  _mm_mul_pd(r, _mm_set1_pd(desiredAcceleration_mPerS2[0]))));
          r = _mm_loadu_pd(&c_P[0]);
          r1 = _mm_loadu_pd(&b_a[0]);
          _mm_storeu_pd(
              &b_a[0],
              _mm_add_pd(
                  r1, _mm_mul_pd(r, _mm_set1_pd(inverseContext->measured[2]))));
          r = _mm_loadu_pd(&b_P[2]);
          r1 = _mm_loadu_pd(&forceBar[0]);
          _mm_storeu_pd(
              &forceBar[0],
              _mm_add_pd(
                  r1,
                  _mm_mul_pd(r, _mm_set1_pd(desiredAcceleration_mPerS2[1]))));
          r = _mm_loadu_pd(&c_P[2]);
          r1 = _mm_loadu_pd(&b_a[0]);
          _mm_storeu_pd(
              &b_a[0],
              _mm_add_pd(
                  r1, _mm_mul_pd(r, _mm_set1_pd(inverseContext->measured[3]))));
          target[0] =
              (forceBar[0] + b_a[0]) + (1300.0 * inverseContext->measured[0] +
                                        650.0 * inverseContext->measured[1]);
          target[1] =
              (forceBar[1] + b_a[1]) + (650.0 * inverseContext->measured[0] +
                                        1300.0 * inverseContext->measured[1]);
          forceBar[0] =
              500000.0 *
                  ((tanh(2000.0 * (inverseContext->measured[0] - 0.035)) +
                    2.0) -
                   tanh(2000.0 * (inverseContext->measured[0] - -0.035))) *
                  inverseContext->measured[0] -
              s;
          forceBar[1] =
              500000.0 *
                  ((tanh(2000.0 * (inverseContext->measured[1] - 0.035)) +
                    2.0) -
                   tanh(2000.0 * (inverseContext->measured[1] - -0.035))) *
                  inverseContext->measured[1] -
              s;
          r = _mm_loadu_pd(&target[0]);
          r1 = _mm_loadu_pd(&forceBar[0]);
          r2 = _mm_loadu_pd(&P[7]);
          _mm_storeu_pd(&target[0], _mm_add_pd(_mm_add_pd(r, r1), r2));
          /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure candidate
           * evaluation. */
          /*  x=[q2;q3;v2;v3;p1;p2;p3;previous physical dt]. Jp holds history
           * fixed. */
          r = _mm_loadu_pd(&network_mean[0]);
          r1 = _mm_loadu_pd(&network_std[0]);
          _mm_storeu_pd(
              &z[0],
              _mm_div_pd(
                  _mm_sub_pd(_mm_loadu_pd(&inverseContext->measured[0]), r),
                  r1));
          r = _mm_loadu_pd(&network_mean[2]);
          r1 = _mm_loadu_pd(&network_std[2]);
          _mm_storeu_pd(
              &z[2],
              _mm_div_pd(
                  _mm_sub_pd(_mm_loadu_pd(&inverseContext->measured[2]), r),
                  r1));
          r = _mm_loadu_pd(&network_mean[4]);
          r1 = _mm_loadu_pd(&network_std[4]);
          _mm_storeu_pd(
              &z[4],
              _mm_div_pd(
                  _mm_sub_pd(_mm_loadu_pd(&inverseContext->pressure_bar[0]), r),
                  r1));
          z[6] = (inverseContext->pressure_bar[2] - network_mean[6]) /
                 network_std[6];
          z[7] =
              (inverseContext->featureDt_s - network_mean[7]) / network_std[7];
          for (k = 0; k < 384; k++) {
            s = 0.0;
            for (i = 0; i < 8; i++) {
              s += network_input[k + 384 * i] * z[i];
            }
            u[k] = s + w[k];
          }
          memset(&w[0], 0, 384U * sizeof(double));
          for (k = 0; k < 128; k++) {
            for (i = 0; i <= 382; i += 2) {
              r = _mm_loadu_pd(&network_recurrent[i + 384 * k]);
              r1 = _mm_loadu_pd(&w[i]);
              _mm_storeu_pd(
                  &w[i],
                  _mm_add_pd(
                      r1, _mm_mul_pd(
                              r, _mm_set1_pd(inverseContext->historyPre[k]))));
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
            _mm_storeu_pd(&b_a[0],
                          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(s))));
          }
          memset(&forceBar[0], 0, sizeof(double) << 1);
          for (k = 0; k < 8; k++) {
            r = _mm_loadu_pd(&network_skip[k << 1]);
            r1 = _mm_loadu_pd(&forceBar[0]);
            _mm_storeu_pd(&forceBar[0],
                          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(z[k]))));
          }
          singlePassDiagnostic->networkEvaluations = 1U;
          /*  Same minimum-nonnegative conversion as the collection inverse. No
           * preload. */
          l_idx_1 =
              80.0 *
              tanh(10.0 * ((b_a[0] + forceBar[0]) + network_outputBias_idx_0) /
                   80.0);
          b_a[0] = l_idx_1;
          s = target[0] - l_idx_1;
          request[0] = s;
          forceBar[0] = s / shadowPressureOptions[9];
          l_idx_1 =
              80.0 *
              tanh(10.0 * ((b_a[1] + forceBar[1]) + network_outputBias_idx_1) /
                   80.0);
          b_a[1] = l_idx_1;
          s = target[1] - l_idx_1;
          request[1] = s;
          forceBar[1] = s / shadowPressureOptions[9];
          y[0] = 0.0;
          y[1] = -forceBar[0];
          y[2] = -forceBar[1];
          s = maximum(y);
          corrected[0] = fmax(s, 0.0);
          corrected[1] = fmax(s + forceBar[0], 0.0);
          corrected[2] = fmax(s + forceBar[1], 0.0);
          s = maximum(corrected);
          allocationLimited = (s > shadowPressureOptions[4]);
          if (allocationLimited) {
            s = shadowPressureOptions[4] / s;
            r = _mm_loadu_pd(&corrected[0]);
            _mm_storeu_pd(&corrected[0], _mm_mul_pd(r, _mm_set1_pd(s)));
            corrected[2] *= s;
          }
          if (corrected[0] <= 1.0E-8) {
            corrected[0] = 0.0;
          }
          if (corrected[1] <= 1.0E-8) {
            corrected[1] = 0.0;
          }
          if (corrected[2] <= 1.0E-8) {
            corrected[2] = 0.0;
          }
          r = _mm_loadu_pd(&target[0]);
          _mm_storeu_pd(&forceBar[0],
                        _mm_div_pd(r, _mm_set1_pd(shadowPressureOptions[9])));
          y[0] = 0.0;
          y[1] = -forceBar[0];
          y[2] = -forceBar[1];
          s = maximum(y);
          y[0] = fmax(s, 0.0);
          y[1] = fmax(s + forceBar[0], 0.0);
          y[2] = fmax(s + forceBar[1], 0.0);
          s = maximum(y);
          nominalLimited = (s > shadowPressureOptions[4]);
          if (nominalLimited) {
            s = shadowPressureOptions[4] / s;
            r = _mm_loadu_pd(&y[0]);
            _mm_storeu_pd(&y[0], _mm_mul_pd(r, _mm_set1_pd(s)));
            y[2] *= s;
          }
          if (y[0] <= 1.0E-8) {
            y[0] = 0.0;
          }
          if (y[1] <= 1.0E-8) {
            y[1] = 0.0;
          }
          if (y[2] <= 1.0E-8) {
            y[2] = 0.0;
          }
          s = shadowPressureOptions[7] * inverseContext->horizon_s;
          r = _mm_loadu_pd(&inverseContext->pressure_bar[0]);
          r1 = _mm_set1_pd(s);
          _mm_storeu_pd(&dv3[0], _mm_sub_pd(r, r1));
          _mm_storeu_pd(&dv4[0], _mm_add_pd(r, r1));
          forceBar[0] = fmin(fmax(corrected[0], dv3[0]), dv4[0]);
          forceBar[1] = fmin(fmax(corrected[1], dv3[1]), dv4[1]);
          r = _mm_loadu_pd(&forceBar[0]);
          _mm_storeu_pd(&p[0], r);
          forceBar[0] = fmin(fmax(y[0], dv3[0]), dv4[0]);
          forceBar[1] = fmin(fmax(y[1], dv3[1]), dv4[1]);
          r = _mm_loadu_pd(&forceBar[0]);
          _mm_storeu_pd(&pn[0], r);
          as = inverseContext->pressure_bar[2] - s;
          s += inverseContext->pressure_bar[2];
          p[2] = fmin(fmax(corrected[2], as), s);
          pn[2] = fmin(fmax(y[2], as), s);
          b_x[0] = target[0];
          b_x[2] = b_a[0];
          b_x[1] = target[1];
          b_x[3] = l_idx_1;
          b_x[4] = corrected[0];
          b_x[7] = y[0];
          b_x[10] = p[0];
          b_x[13] = pn[0];
          b_x[5] = corrected[1];
          b_x[8] = y[1];
          b_x[11] = p[1];
          b_x[14] = pn[1];
          b_x[6] = corrected[2];
          b_x[9] = y[2];
          b_x[12] = p[2];
          b_x[15] = pn[2];
          for (k = 0; k < 16; k++) {
            s = b_x[k];
            bv1[k] = ((!rtIsInf(s)) && (!rtIsNaN(s)));
          }
          if (b_all(bv1)) {
            for (k = 0; k < 6; k++) {
              B[k] = shadowPressureOptions[9] * (double)b_b[k];
            }
            __m128d r3;
            __m128d r4;
            __m128d r5;
            __m128d r6;
            singlePassDiagnostic->valid = true;
            singlePassDiagnostic->allocationLimited = allocationLimited;
            singlePassDiagnostic->nominalAllocationLimited = nominalLimited;
            c_singlePassPressurePrediction_[0] = p[0];
            nominalPressurePrediction_bar[0] = pn[0];
            s = p[0] - corrected[0];
            l_idx_0 = s;
            b_inverseContext[0] = (fabs(s) > 1.0E-9);
            c_singlePassPressurePrediction_[1] = p[1];
            nominalPressurePrediction_bar[1] = pn[1];
            s = p[1] - corrected[1];
            l_idx_1 = s;
            b_inverseContext[1] = (fabs(s) > 1.0E-9);
            c_singlePassPressurePrediction_[2] = p[2];
            nominalPressurePrediction_bar[2] = pn[2];
            s = p[2] - corrected[2];
            b_inverseContext[2] = (fabs(s) > 1.0E-9);
            singlePassDiagnostic->slewLimited = any(b_inverseContext);
            b_inverseContext[0] = (fabs(pn[0] - y[0]) > 1.0E-9);
            b_inverseContext[1] = (fabs(pn[1] - y[1]) > 1.0E-9);
            b_inverseContext[2] = (fabs(pn[2] - y[2]) > 1.0E-9);
            singlePassDiagnostic->nominalSlewLimited = any(b_inverseContext);
            as = p[0];
            bs = p[1];
            a2_tmp = p[2];
            c1 = corrected[0];
            P_tmp = corrected[1];
            c1_tmp_tmp = corrected[2];
            r = _mm_loadu_pd(&target[0]);
            _mm_storeu_pd(&singlePassDiagnostic->targetForce_N[0], r);
            r3 = _mm_loadu_pd(&b_a[0]);
            _mm_storeu_pd(&singlePassDiagnostic->residualForce_N[0], r3);
            r4 = _mm_loadu_pd(&request[0]);
            _mm_storeu_pd(&singlePassDiagnostic->correctedForceRequest_N[0],
                          r4);
            r = _mm_loadu_pd(&B[0]);
            r1 = _mm_mul_pd(r, _mm_set1_pd(as));
            r5 = _mm_mul_pd(r, _mm_set1_pd(c1));
            r2 = _mm_loadu_pd(&B[2]);
            r = _mm_mul_pd(r2, _mm_set1_pd(bs));
            r6 = _mm_add_pd(r1, r);
            r = _mm_mul_pd(r2, _mm_set1_pd(P_tmp));
            r5 = _mm_add_pd(r5, r);
            r1 = _mm_loadu_pd(&B[4]);
            r = _mm_mul_pd(r1, _mm_set1_pd(a2_tmp));
            r2 = _mm_add_pd(r6, r);
            r = _mm_mul_pd(r1, _mm_set1_pd(c1_tmp_tmp));
            r1 = _mm_add_pd(r5, r);
            r = _mm_add_pd(r2, r3);
            _mm_storeu_pd(&forceBar[0], r);
            _mm_storeu_pd(&singlePassDiagnostic->frozenTotalForce_N[0], r);
            r = _mm_sub_pd(r1, r4);
            _mm_storeu_pd(&singlePassDiagnostic->allocationError_N[0], r);
            _mm_storeu_pd(&singlePassDiagnostic->slewError_N[0],
                          _mm_set1_pd(0.0));
            r = _mm_loadu_pd(&B[0]);
            r1 = _mm_loadu_pd(&singlePassDiagnostic->slewError_N[0]);
            _mm_storeu_pd(&singlePassDiagnostic->slewError_N[0],
                          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(l_idx_0))));
            r = _mm_loadu_pd(&B[2]);
            r1 = _mm_loadu_pd(&singlePassDiagnostic->slewError_N[0]);
            _mm_storeu_pd(&singlePassDiagnostic->slewError_N[0],
                          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(l_idx_1))));
            r = _mm_loadu_pd(&B[4]);
            r1 = _mm_loadu_pd(&singlePassDiagnostic->slewError_N[0]);
            _mm_storeu_pd(&singlePassDiagnostic->slewError_N[0],
                          _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(s))));
            s = pn[0];
            as = pn[1];
            l_idx_1 = pn[2];
            r = _mm_loadu_pd(&forceBar[0]);
            r2 = _mm_loadu_pd(&target[0]);
            r = _mm_sub_pd(r, r2);
            _mm_storeu_pd(&singlePassDiagnostic->frozenBalanceError_N[0], r);
            r = _mm_loadu_pd(&B[0]);
            r1 = _mm_mul_pd(r, _mm_set1_pd(s));
            r = _mm_loadu_pd(&B[2]);
            r = _mm_mul_pd(r, _mm_set1_pd(as));
            r1 = _mm_add_pd(r1, r);
            r = _mm_loadu_pd(&B[4]);
            r = _mm_mul_pd(r, _mm_set1_pd(l_idx_1));
            r = _mm_add_pd(r1, r);
            r = _mm_sub_pd(r, r2);
            _mm_storeu_pd(&singlePassDiagnostic->nominalBalanceError_N[0], r);
            singlePassDiagnostic->frozenBalanceSatisfied =
                (b_norm(singlePassDiagnostic->frozenBalanceError_N) <=
                 singlePassOptions[1]);
            _mm_storeu_pd(&singlePassDiagnostic->residualInputPressure_bar[0],
                          _mm_loadu_pd(&inverseContext->pressure_bar[0]));
            r = _mm_loadu_pd(&corrected[0]);
            _mm_storeu_pd(
                &singlePassDiagnostic->correctedPressureBeforeSlew_bar[0], r);
            r1 = _mm_loadu_pd(&y[0]);
            _mm_storeu_pd(
                &singlePassDiagnostic->nominalPressureBeforeSlew_bar[0], r1);
            _mm_storeu_pd(
                &singlePassDiagnostic->c_pressureDifferenceBeforeSlew_[0],
                _mm_sub_pd(r, r1));
            r = _mm_loadu_pd(&p[0]);
            r1 = _mm_loadu_pd(&pn[0]);
            _mm_storeu_pd(&singlePassDiagnostic->pressureDifference_bar[0],
                          _mm_sub_pd(r, r1));
            singlePassDiagnostic->residualInputPressure_bar[2] =
                inverseContext->pressure_bar[2];
            singlePassDiagnostic->correctedPressureBeforeSlew_bar[2] =
                corrected[2];
            singlePassDiagnostic->nominalPressureBeforeSlew_bar[2] = y[2];
            singlePassDiagnostic->c_pressureDifferenceBeforeSlew_[2] =
                corrected[2] - y[2];
            singlePassDiagnostic->pressureDifference_bar[2] = p[2] - pn[2];
            singlePassDiagnostic->sampleTime_host_s =
                inverseContext->sampleTime_host_s;
            singlePassDiagnostic->horizon_s = inverseContext->horizon_s;
          }
        }
      }
    }
  }
}

/*
 * File trailer for rlResidualSinglePassCollectionShadowBlock.c
 *
 * [EOF]
 */
