/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualDesiredPressureBlock.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 09-Nov-2026 11:59:05
 */

/* Include Files */
#include "rlResidualDesiredPressureBlock.h"
#include "all.h"
#include "any.h"
#include "minOrMax.h"
#include "norm.h"
#include "rlResidualDesiredPressureBlock_data.h"
#include "rlResidualDesiredPressureBlock_initialize.h"
#include "rlResidualDesiredPressureBlock_types.h"
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
 * RLRESIDUALDESIREDPRESSUREBLOCK Paste into a discrete 0.05 s MATLAB Function.
 *  desiredMotion = [l2d;l3d;dl2d;dl3d;ddl2d;ddl3d], double 6x1, SI units.
 *  Feedforward pressure for DESIRED motion, with a frozen observed residual.
 *  One GRU evaluation at actual measured q/v, actual command and causal
 * history. No recurrent-state update, inverse iterations or integration in this
 * block. This is a single-pass inverse approximation, not an exact nonlinear
 * solve. Initial installation: Scopes/logging only; existing valve/Kill wiring
 * stays. Diagnostic output type: Bus: ResidualSinglePassDiagnosticsBus
 * (unchanged). variableCommon is retained for port compatibility and is not
 * used. Pressure options: 4=min(0), 5=max, 8=slew, 9=activation(0), 10=N/bar.
 *
 * Arguments    : const struct0_T *inverseContext
 *                const double desiredMotion[6]
 *                const double shadowPressureOptions[10]
 *                boolean_T variableCommon
 *                const double singlePassOptions[2]
 *                double desiredPressurePrediction_bar[3]
 *                double nominalPressurePrediction_bar[3]
 *                struct1_T *desiredPressureDiagnostic
 * Return Type  : void
 */
void rlResidualDesiredPressureBlock(const struct0_T *inverseContext,
                                    const double desiredMotion[6],
                                    const double shadowPressureOptions[10],
                                    boolean_T variableCommon,
                                    const double singlePassOptions[2],
                                    double desiredPressurePrediction_bar[3],
                                    double nominalPressurePrediction_bar[3],
                                    struct1_T *desiredPressureDiagnostic)
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
  double x[144];
  double b_x[16];
  double network_skip[16];
  double network_mean[8];
  double network_std[8];
  double z[8];
  double B[6];
  double P[4];
  double b_P[4];
  double corrected[3];
  double p[3];
  double pn[3];
  double y[3];
  double a[2];
  double dv3[2];
  double dv4[2];
  double forceBar[2];
  double request[2];
  double target[2];
  int b_i;
  int i;
  (void)variableCommon;
  if (!isInitialized_rlResidualDesiredPressureBlock) {
    rlResidualDesiredPressureBlock_initialize();
  }
  desiredPressurePrediction_bar[0] = rtNaN;
  nominalPressurePrediction_bar[0] = rtNaN;
  desiredPressurePrediction_bar[1] = rtNaN;
  nominalPressurePrediction_bar[1] = rtNaN;
  desiredPressurePrediction_bar[2] = rtNaN;
  nominalPressurePrediction_bar[2] = rtNaN;
  desiredPressureDiagnostic->valid = false;
  desiredPressureDiagnostic->networkEvaluations = 0U;
  desiredPressureDiagnostic->allocationLimited = false;
  desiredPressureDiagnostic->slewLimited = false;
  desiredPressureDiagnostic->nominalAllocationLimited = false;
  desiredPressureDiagnostic->nominalSlewLimited = false;
  desiredPressureDiagnostic->frozenBalanceSatisfied = false;
  desiredPressureDiagnostic->targetForce_N[0] = rtNaN;
  desiredPressureDiagnostic->residualForce_N[0] = rtNaN;
  desiredPressureDiagnostic->correctedForceRequest_N[0] = rtNaN;
  desiredPressureDiagnostic->frozenTotalForce_N[0] = rtNaN;
  desiredPressureDiagnostic->allocationError_N[0] = rtNaN;
  desiredPressureDiagnostic->slewError_N[0] = rtNaN;
  desiredPressureDiagnostic->frozenBalanceError_N[0] = rtNaN;
  desiredPressureDiagnostic->nominalBalanceError_N[0] = rtNaN;
  desiredPressureDiagnostic->targetForce_N[1] = rtNaN;
  desiredPressureDiagnostic->residualForce_N[1] = rtNaN;
  desiredPressureDiagnostic->correctedForceRequest_N[1] = rtNaN;
  desiredPressureDiagnostic->frozenTotalForce_N[1] = rtNaN;
  desiredPressureDiagnostic->allocationError_N[1] = rtNaN;
  desiredPressureDiagnostic->slewError_N[1] = rtNaN;
  desiredPressureDiagnostic->frozenBalanceError_N[1] = rtNaN;
  desiredPressureDiagnostic->nominalBalanceError_N[1] = rtNaN;
  desiredPressureDiagnostic->residualInputPressure_bar[0] = rtNaN;
  desiredPressureDiagnostic->correctedPressureBeforeSlew_bar[0] = rtNaN;
  desiredPressureDiagnostic->nominalPressureBeforeSlew_bar[0] = rtNaN;
  desiredPressureDiagnostic->c_pressureDifferenceBeforeSlew_[0] = rtNaN;
  desiredPressureDiagnostic->pressureDifference_bar[0] = rtNaN;
  desiredPressureDiagnostic->residualInputPressure_bar[1] = rtNaN;
  desiredPressureDiagnostic->correctedPressureBeforeSlew_bar[1] = rtNaN;
  desiredPressureDiagnostic->nominalPressureBeforeSlew_bar[1] = rtNaN;
  desiredPressureDiagnostic->c_pressureDifferenceBeforeSlew_[1] = rtNaN;
  desiredPressureDiagnostic->pressureDifference_bar[1] = rtNaN;
  desiredPressureDiagnostic->residualInputPressure_bar[2] = rtNaN;
  desiredPressureDiagnostic->correctedPressureBeforeSlew_bar[2] = rtNaN;
  desiredPressureDiagnostic->nominalPressureBeforeSlew_bar[2] = rtNaN;
  desiredPressureDiagnostic->c_pressureDifferenceBeforeSlew_[2] = rtNaN;
  desiredPressureDiagnostic->pressureDifference_bar[2] = rtNaN;
  desiredPressureDiagnostic->sampleTime_host_s = rtNaN;
  desiredPressureDiagnostic->horizon_s = rtNaN;
  if (inverseContext->valid) {
    double s;
    boolean_T bv[144];
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
    for (i = 0; i < 6; i++) {
      x[i + 138] = desiredMotion[i];
    }
    for (i = 0; i < 144; i++) {
      s = x[i];
      bv[i] = ((!rtIsInf(s)) && (!rtIsNaN(s)));
    }
    if (b_all(bv) && (!(inverseContext->featureDt_s <= 0.0)) &&
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
          double P_idx_1;
          double a_tmp;
          double actuatorForce_N_idx_0;
          double b;
          double b_a2_tmp;
          double bs;
          double c1;
          double c1_tmp_tmp;
          double l_idx_0;
          double w1;
          double w2;
          boolean_T bv1[2];
          /*  Both branches target precisely the SAME desired trajectory. */
          /* FORCERESIDUALPHYSICS Batched differentiable EXACT nominal M,C,G
           * polynomials. */
          /*  q,v:2xbatch, numeric or unformatted dlarray. P rows M11,M12,M22,
           */
          /*  C11,C12,C21,C22,G1,G2. Fixed nominal basis, no fitted
           * approximation. */
          /*  Generated rational coefficients; see
           * generate_force_polynomials.py. */
          w1 = 7889.5463510848131 * desiredMotion[0] +
               3944.7731755424065 * desiredMotion[1];
          w2 = 3944.7731755424065 * desiredMotion[0] +
               7889.5463510848131 * desiredMotion[1];
          s = desiredMotion[0] * w1 + desiredMotion[1] * w2;
          a_tmp = 0.0 * s;
          P_idx_1 = a_tmp + 1.743355929321998E-19;
          c1 = a_tmp;
          b = a_tmp + 1.3754287410824441E-19;
          bs = a_tmp;
          for (i = 0; i < 10; i++) {
            c1 = c1 * s + P_idx_1;
            P_idx_1 = P_idx_1 * s + dv[9 - i];
            bs = bs * s + b;
            b = b * s + dv1[9 - i];
          }
          actuatorForce_N_idx_0 = a_tmp + 1.1470745597729725E-11;
          for (i = 0; i < 6; i++) {
            a_tmp = a_tmp * s + actuatorForce_N_idx_0;
            actuatorForce_N_idx_0 = actuatorForce_N_idx_0 * s + dv2[5 - i];
          }
          double a2;
          double a2_tmp;
          double b_c1_tmp;
          double c1_tmp;
          double c_a2_tmp;
          double c_c1_tmp;
          double d_a2_tmp;
          double e_a2_tmp;
          double f_a2_tmp;
          actuatorForce_N_idx_0 = 2.0 * c1;
          a2_tmp = 2.0 * bs;
          b_a2_tmp = 2.0 * b;
          c_a2_tmp = w1 * w1;
          d_a2_tmp = b_a2_tmp * w1;
          e_a2_tmp = actuatorForce_N_idx_0 * w2;
          f_a2_tmp = e_a2_tmp * 7889.5463510848131;
          a2 = 0.0030314292099999996 * ((f_a2_tmp + a2_tmp * w2 * c_a2_tmp) +
                                        d_a2_tmp * 3944.7731755424065);
          c1_tmp = w2 * w2;
          c1_tmp_tmp = actuatorForce_N_idx_0 * w1;
          b_c1_tmp = c1_tmp_tmp * 7889.5463510848131;
          c_c1_tmp = a2_tmp * w1 * c1_tmp;
          s = b_a2_tmp * w2;
          c1 = 0.0030314292099999996 *
               ((b_c1_tmp + c_c1_tmp) + s * 3944.7731755424065);
          bs = P_idx_1 * 7889.5463510848131;
          P_idx_1 = 0.0030314292099999996 *
                    (P_idx_1 * 3944.7731755424065 + b * w1 * w2);
          a_tmp *= -0.34160382;
          /* ARMS_SINGLE_BOUNDFORCE Soft actuator-length restoring force in
           * q=[l2;l3]. */
          /*  Default mode1: all l=[-l2-l3;l2;l3] use the same bounds, in
           * metres. */
          /*  Mode0 reproduces historical coordinate-only l2/l3 bounds for old
           * analyses. */
          /*  coefficient*q equals force_N; coefficient is NOT the tangent
           * stiffness. */
          l_idx_0 = -desiredMotion[0] - desiredMotion[1];
          actuatorForce_N_idx_0 = 500000.0 *
                                  ((tanh(2000.0 * (l_idx_0 - 0.035)) + 2.0) -
                                   tanh(2000.0 * (l_idx_0 - -0.035))) *
                                  l_idx_0;
          /*  J=dl/dq=[-1 -1;1 0;0 1]; generalized force is J'*actuatorForce. */
          P[0] = 0.0030314292099999996 * (bs + b * c_a2_tmp);
          P[2] = P_idx_1;
          P[1] = P_idx_1;
          P[3] = 0.0030314292099999996 * (bs + b * c1_tmp);
          b_P[0] = 0.5 * (0.0030314292099999996 *
                              ((b_c1_tmp + a2_tmp * rt_powd_snf(w1, 3.0)) +
                               d_a2_tmp * 7889.5463510848131) *
                              desiredMotion[2] +
                          a2 * desiredMotion[3]) +
                   40.0;
          b_P[2] = 0.5 * (a2 * desiredMotion[2] +
                          (2.0 * (0.0030314292099999996 *
                                  ((e_a2_tmp * 3944.7731755424065 + c_c1_tmp) +
                                   b * (3944.7731755424065 * w2 +
                                        7889.5463510848131 * w1))) -
                           c1) *
                              desiredMotion[3]);
          b_P[1] = 0.5 * ((2.0 * (0.0030314292099999996 *
                                  ((c1_tmp_tmp * 3944.7731755424065 +
                                    a2_tmp * c_a2_tmp * w2) +
                                   b * (7889.5463510848131 * w2 +
                                        3944.7731755424065 * w1))) -
                           a2) *
                              desiredMotion[2] +
                          c1 * desiredMotion[3]);
          b_P[3] = 0.5 * (c1 * desiredMotion[2] +
                          0.0030314292099999996 *
                              ((f_a2_tmp + a2_tmp * rt_powd_snf(w2, 3.0)) +
                               s * 7889.5463510848131) *
                              desiredMotion[3]) +
                   40.0;
          memset(&forceBar[0], 0, sizeof(double) << 1);
          memset(&a[0], 0, sizeof(double) << 1);
          r = _mm_loadu_pd(&P[0]);
          r1 = _mm_loadu_pd(&forceBar[0]);
          _mm_storeu_pd(
              &forceBar[0],
              _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(desiredMotion[4]))));
          r = _mm_loadu_pd(&b_P[0]);
          r1 = _mm_loadu_pd(&a[0]);
          _mm_storeu_pd(
              &a[0],
              _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(desiredMotion[2]))));
          r = _mm_loadu_pd(&P[2]);
          r1 = _mm_loadu_pd(&forceBar[0]);
          _mm_storeu_pd(
              &forceBar[0],
              _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(desiredMotion[5]))));
          r = _mm_loadu_pd(&b_P[2]);
          r1 = _mm_loadu_pd(&a[0]);
          _mm_storeu_pd(
              &a[0],
              _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(desiredMotion[3]))));
          target[0] = (forceBar[0] + a[0]) +
                      (1300.0 * desiredMotion[0] + 650.0 * desiredMotion[1]);
          target[1] = (forceBar[1] + a[1]) +
                      (650.0 * desiredMotion[0] + 1300.0 * desiredMotion[1]);
          c1_tmp_tmp =
              (target[0] +
               (500000.0 *
                    ((tanh(2000.0 * (desiredMotion[0] - 0.035)) + 2.0) -
                     tanh(2000.0 * (desiredMotion[0] - -0.035))) *
                    desiredMotion[0] -
                actuatorForce_N_idx_0)) +
              a_tmp * w1;
          target[0] = c1_tmp_tmp;
          bv1[0] = ((!rtIsInf(c1_tmp_tmp)) && (!rtIsNaN(c1_tmp_tmp)));
          c1_tmp_tmp =
              (target[1] +
               (500000.0 *
                    ((tanh(2000.0 * (desiredMotion[1] - 0.035)) + 2.0) -
                     tanh(2000.0 * (desiredMotion[1] - -0.035))) *
                    desiredMotion[1] -
                actuatorForce_N_idx_0)) +
              a_tmp * w2;
          target[1] = c1_tmp_tmp;
          bv1[1] = ((!rtIsInf(c1_tmp_tmp)) && (!rtIsNaN(c1_tmp_tmp)));
          if (all(bv1)) {
            double u[384];
            boolean_T bv2[16];
            boolean_T allocationLimited;
            boolean_T nominalLimited;
            /*  History and GRU features remain ACTUAL observations, not desired
             * motion. */
            /* RESIDUALGRUCELL Frozen single-layer GRU, padded128; pure
             * candidate evaluation. */
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
                    _mm_sub_pd(_mm_loadu_pd(&inverseContext->pressure_bar[0]),
                               r),
                    r1));
            z[6] = (inverseContext->pressure_bar[2] - network_mean[6]) /
                   network_std[6];
            z[7] = (inverseContext->featureDt_s - network_mean[7]) /
                   network_std[7];
            for (i = 0; i < 384; i++) {
              s = 0.0;
              for (b_i = 0; b_i < 8; b_i++) {
                s += network_input[i + 384 * b_i] * z[b_i];
              }
              u[i] = s + w[i];
            }
            memset(&w[0], 0, 384U * sizeof(double));
            for (i = 0; i < 128; i++) {
              for (b_i = 0; b_i <= 382; b_i += 2) {
                r = _mm_loadu_pd(&network_recurrent[b_i + 384 * i]);
                r1 = _mm_loadu_pd(&w[b_i]);
                _mm_storeu_pd(
                    &w[b_i],
                    _mm_add_pd(
                        r1, _mm_mul_pd(r, _mm_set1_pd(
                                              inverseContext->historyPre[i]))));
              }
            }
            memset(&a[0], 0, sizeof(double) << 1);
            for (i = 0; i < 128; i++) {
              s = 0.5 * (tanh(0.5 * (u[i + 128] + w[i + 128])) + 1.0);
              s = s * inverseContext->historyPre[i] +
                  (1.0 - s) * tanh(u[i + 256] +
                                   0.5 * (tanh(0.5 * (u[i] + w[i])) + 1.0) *
                                       w[i + 256]);
              r = _mm_loadu_pd(&network_output[i << 1]);
              r1 = _mm_loadu_pd(&a[0]);
              _mm_storeu_pd(&a[0],
                            _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(s))));
            }
            memset(&forceBar[0], 0, sizeof(double) << 1);
            for (i = 0; i < 8; i++) {
              r = _mm_loadu_pd(&network_skip[i << 1]);
              r1 = _mm_loadu_pd(&forceBar[0]);
              _mm_storeu_pd(&forceBar[0],
                            _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(z[i]))));
            }
            desiredPressureDiagnostic->networkEvaluations = 1U;
            actuatorForce_N_idx_0 =
                80.0 *
                tanh(10.0 * ((a[0] + forceBar[0]) + network_outputBias_idx_0) /
                     80.0);
            a[0] = actuatorForce_N_idx_0;
            s = target[0] - actuatorForce_N_idx_0;
            request[0] = s;
            forceBar[0] = s / shadowPressureOptions[9];
            actuatorForce_N_idx_0 =
                80.0 *
                tanh(10.0 * ((a[1] + forceBar[1]) + network_outputBias_idx_1) /
                     80.0);
            a[1] = actuatorForce_N_idx_0;
            s = c1_tmp_tmp - actuatorForce_N_idx_0;
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
            c1 = inverseContext->pressure_bar[2] - s;
            s += inverseContext->pressure_bar[2];
            p[2] = fmin(fmax(corrected[2], c1), s);
            pn[2] = fmin(fmax(y[2], c1), s);
            b_x[0] = target[0];
            b_x[2] = a[0];
            b_x[1] = c1_tmp_tmp;
            b_x[3] = actuatorForce_N_idx_0;
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
            for (i = 0; i < 16; i++) {
              s = b_x[i];
              bv2[i] = ((!rtIsInf(s)) && (!rtIsNaN(s)));
            }
            if (c_all(bv2)) {
              for (i = 0; i < 6; i++) {
                B[i] = shadowPressureOptions[9] * (double)b_b[i];
              }
              __m128d r2;
              __m128d r3;
              __m128d r4;
              __m128d r5;
              __m128d r6;
              desiredPressureDiagnostic->valid = true;
              desiredPressureDiagnostic->allocationLimited = allocationLimited;
              desiredPressureDiagnostic->nominalAllocationLimited =
                  nominalLimited;
              desiredPressurePrediction_bar[0] = p[0];
              nominalPressurePrediction_bar[0] = pn[0];
              s = p[0] - corrected[0];
              l_idx_0 = s;
              b_inverseContext[0] = (fabs(s) > 1.0E-9);
              desiredPressurePrediction_bar[1] = p[1];
              nominalPressurePrediction_bar[1] = pn[1];
              s = p[1] - corrected[1];
              actuatorForce_N_idx_0 = s;
              b_inverseContext[1] = (fabs(s) > 1.0E-9);
              desiredPressurePrediction_bar[2] = p[2];
              nominalPressurePrediction_bar[2] = pn[2];
              s = p[2] - corrected[2];
              b_inverseContext[2] = (fabs(s) > 1.0E-9);
              desiredPressureDiagnostic->slewLimited = any(b_inverseContext);
              b_inverseContext[0] = (fabs(pn[0] - y[0]) > 1.0E-9);
              b_inverseContext[1] = (fabs(pn[1] - y[1]) > 1.0E-9);
              b_inverseContext[2] = (fabs(pn[2] - y[2]) > 1.0E-9);
              desiredPressureDiagnostic->nominalSlewLimited =
                  any(b_inverseContext);
              c1 = p[0];
              bs = p[1];
              b_a2_tmp = p[2];
              c1_tmp_tmp = corrected[0];
              P_idx_1 = corrected[1];
              a_tmp = corrected[2];
              r = _mm_loadu_pd(&target[0]);
              _mm_storeu_pd(&desiredPressureDiagnostic->targetForce_N[0], r);
              r2 = _mm_loadu_pd(&a[0]);
              _mm_storeu_pd(&desiredPressureDiagnostic->residualForce_N[0], r2);
              r3 = _mm_loadu_pd(&request[0]);
              _mm_storeu_pd(
                  &desiredPressureDiagnostic->correctedForceRequest_N[0], r3);
              r = _mm_loadu_pd(&B[0]);
              r1 = _mm_mul_pd(r, _mm_set1_pd(c1));
              r4 = _mm_mul_pd(r, _mm_set1_pd(c1_tmp_tmp));
              r5 = _mm_loadu_pd(&B[2]);
              r = _mm_mul_pd(r5, _mm_set1_pd(bs));
              r6 = _mm_add_pd(r1, r);
              r = _mm_mul_pd(r5, _mm_set1_pd(P_idx_1));
              r4 = _mm_add_pd(r4, r);
              r1 = _mm_loadu_pd(&B[4]);
              r = _mm_mul_pd(r1, _mm_set1_pd(b_a2_tmp));
              r5 = _mm_add_pd(r6, r);
              r = _mm_mul_pd(r1, _mm_set1_pd(a_tmp));
              r1 = _mm_add_pd(r4, r);
              r = _mm_add_pd(r5, r2);
              _mm_storeu_pd(&forceBar[0], r);
              _mm_storeu_pd(&desiredPressureDiagnostic->frozenTotalForce_N[0],
                            r);
              r = _mm_sub_pd(r1, r3);
              _mm_storeu_pd(&desiredPressureDiagnostic->allocationError_N[0],
                            r);
              _mm_storeu_pd(&desiredPressureDiagnostic->slewError_N[0],
                            _mm_set1_pd(0.0));
              r = _mm_loadu_pd(&B[0]);
              r1 = _mm_loadu_pd(&desiredPressureDiagnostic->slewError_N[0]);
              _mm_storeu_pd(
                  &desiredPressureDiagnostic->slewError_N[0],
                  _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(l_idx_0))));
              r = _mm_loadu_pd(&B[2]);
              r1 = _mm_loadu_pd(&desiredPressureDiagnostic->slewError_N[0]);
              _mm_storeu_pd(
                  &desiredPressureDiagnostic->slewError_N[0],
                  _mm_add_pd(
                      r1, _mm_mul_pd(r, _mm_set1_pd(actuatorForce_N_idx_0))));
              r = _mm_loadu_pd(&B[4]);
              r1 = _mm_loadu_pd(&desiredPressureDiagnostic->slewError_N[0]);
              _mm_storeu_pd(&desiredPressureDiagnostic->slewError_N[0],
                            _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(s))));
              s = pn[0];
              c1 = pn[1];
              actuatorForce_N_idx_0 = pn[2];
              r = _mm_loadu_pd(&forceBar[0]);
              r5 = _mm_loadu_pd(&target[0]);
              r = _mm_sub_pd(r, r5);
              _mm_storeu_pd(&desiredPressureDiagnostic->frozenBalanceError_N[0],
                            r);
              r = _mm_loadu_pd(&B[0]);
              r1 = _mm_mul_pd(r, _mm_set1_pd(s));
              r = _mm_loadu_pd(&B[2]);
              r = _mm_mul_pd(r, _mm_set1_pd(c1));
              r1 = _mm_add_pd(r1, r);
              r = _mm_loadu_pd(&B[4]);
              r = _mm_mul_pd(r, _mm_set1_pd(actuatorForce_N_idx_0));
              r = _mm_add_pd(r1, r);
              r = _mm_sub_pd(r, r5);
              _mm_storeu_pd(
                  &desiredPressureDiagnostic->nominalBalanceError_N[0], r);
              desiredPressureDiagnostic->frozenBalanceSatisfied =
                  (b_norm(desiredPressureDiagnostic->frozenBalanceError_N) <=
                   singlePassOptions[1]);
              _mm_storeu_pd(
                  &desiredPressureDiagnostic->residualInputPressure_bar[0],
                  _mm_loadu_pd(&inverseContext->pressure_bar[0]));
              r = _mm_loadu_pd(&corrected[0]);
              _mm_storeu_pd(&desiredPressureDiagnostic
                                 ->correctedPressureBeforeSlew_bar[0],
                            r);
              r1 = _mm_loadu_pd(&y[0]);
              _mm_storeu_pd(
                  &desiredPressureDiagnostic->nominalPressureBeforeSlew_bar[0],
                  r1);
              _mm_storeu_pd(&desiredPressureDiagnostic
                                 ->c_pressureDifferenceBeforeSlew_[0],
                            _mm_sub_pd(r, r1));
              r = _mm_loadu_pd(&p[0]);
              r1 = _mm_loadu_pd(&pn[0]);
              _mm_storeu_pd(
                  &desiredPressureDiagnostic->pressureDifference_bar[0],
                  _mm_sub_pd(r, r1));
              desiredPressureDiagnostic->residualInputPressure_bar[2] =
                  inverseContext->pressure_bar[2];
              desiredPressureDiagnostic->correctedPressureBeforeSlew_bar[2] =
                  corrected[2];
              desiredPressureDiagnostic->nominalPressureBeforeSlew_bar[2] =
                  y[2];
              desiredPressureDiagnostic->c_pressureDifferenceBeforeSlew_[2] =
                  corrected[2] - y[2];
              desiredPressureDiagnostic->pressureDifference_bar[2] =
                  p[2] - pn[2];
              desiredPressureDiagnostic->sampleTime_host_s =
                  inverseContext->sampleTime_host_s;
              desiredPressureDiagnostic->horizon_s = inverseContext->horizon_s;
            }
          }
        }
      }
    }
  }
}

/*
 * File trailer for rlResidualDesiredPressureBlock.c
 *
 * [EOF]
 */
