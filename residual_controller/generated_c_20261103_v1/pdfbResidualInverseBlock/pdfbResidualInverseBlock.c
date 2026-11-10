/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pdfbResidualInverseBlock.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

/* Include Files */
#include "pdfbResidualInverseBlock.h"
#include "all.h"
#include "armS_single_core.h"
#include "det.h"
#include "pdfbResidualInverseBlock_data.h"
#include "pdfbResidualInverseBlock_initialize.h"
#include "pdfbResidualInverseBlock_internal_types.h"
#include "pdfbResidualInverseBlock_types.h"
#include "residualPressureStep.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Variable Definitions */
static b_struct_T ps;

static c_struct_T ns;

/* Function Definitions */
/*
 * PDFBRESIDUALINVERSEBLOCK Paste-ready model/inverse; outer PDFB/PID unchanged.
 *  First6 inputs/first8 outputs match saved FBL_commonP_circle_BRot inverse.
 *  Append physical measurement time/enable, residual switch/configuration.
 *  Fresh model-off trials in THIS wrapper are the paired nominal baseline.
 *  The old embedded bound/allocator is not claimed bit-identical to this pair.
 *
 * Arguments    : const double measured[6]
 *                const double desiredAcceleration[2]
 *                boolean_T orientationValid
 *                double t
 *                const double gravity[3]
 *                boolean_T baseOrientationValid
 *                double measurementTime_s
 *                boolean_T enablePressure
 *                boolean_T enableResidual
 *                boolean_T enableCommonPressureAdaptation
 *                const double pressureOptions[10]
 *                const double supervisorOptions[7]
 *                const double residualOptions[8]
 *                double pressure_bar[3]
 *                double tauDesired_N[2]
 *                double tauAvailable_N[2]
 *                boolean_T *limited
 *                boolean_T *modelValid
 *                double MassForce_N[2]
 *                double CoriolisForce_N[2]
 *                double GnStiff_N[2]
 *                unsigned char *status
 *                struct0_T *residualDiagnostic
 * Return Type  : void
 */
void pdfbResidualInverseBlock(
    const double measured[6], const double desiredAcceleration[2],
    boolean_T orientationValid, double t, const double gravity[3],
    boolean_T baseOrientationValid, double measurementTime_s,
    boolean_T enablePressure, boolean_T enableResidual,
    boolean_T enableCommonPressureAdaptation, const double pressureOptions[10],
    const double supervisorOptions[7], const double residualOptions[8],
    double pressure_bar[3], double tauDesired_N[2], double tauAvailable_N[2],
    boolean_T *limited, boolean_T *modelValid, double MassForce_N[2],
    double CoriolisForce_N[2], double GnStiff_N[2], unsigned char *status,
    struct0_T *residualDiagnostic)
{
  static d_struct_T network;
  static const double t2_g[3] = {0.0, 0.0, -9.81};
  struct_T expl_temp;
  double M[4];
  double G[2];
  double a__5[2];
  double shift[2];
  boolean_T bv[4];
  boolean_T b_gravity[3];
  boolean_T bv1[2];
  boolean_T good;
  if (!isInitialized_pdfbResidualInverseBlock) {
    pdfbResidualInverseBlock_initialize();
  }
  /* RESIDUALCONTROLLERNETWORK Fixed-size code-generation constants; run packer
   * first. */
  /*  Set units64/128 before a fresh run. Never replace weights during a trial.
   */
  if (residualOptions[0] == 64.0) {
    network = n64;
  } else {
    network = n128;
  }
  /* RESIDUALPRESSURECONFIG Same vector contract as PAC tuning, no feedback
   * gains. */
  tauDesired_N[0] = 0.0;
  MassForce_N[0] = 0.0;
  CoriolisForce_N[0] = 0.0;
  GnStiff_N[0] = 0.0;
  tauDesired_N[1] = 0.0;
  MassForce_N[1] = 0.0;
  CoriolisForce_N[1] = 0.0;
  GnStiff_N[1] = 0.0;
  if (orientationValid && baseOrientationValid) {
    bv[0] = ((!rtIsInf(measured[0])) && (!rtIsNaN(measured[0])));
    bv[1] = ((!rtIsInf(measured[1])) && (!rtIsNaN(measured[1])));
    bv[2] = ((!rtIsInf(measured[2])) && (!rtIsNaN(measured[2])));
    bv[3] = ((!rtIsInf(measured[3])) && (!rtIsNaN(measured[3])));
    if (c_all(bv)) {
      bv1[0] = ((!rtIsInf(desiredAcceleration[0])) &&
                (!rtIsNaN(desiredAcceleration[0])));
      bv1[1] = ((!rtIsInf(desiredAcceleration[1])) &&
                (!rtIsNaN(desiredAcceleration[1])));
      if (d_all(bv1)) {
        b_gravity[0] = ((!rtIsInf(gravity[0])) && (!rtIsNaN(gravity[0])));
        b_gravity[1] = ((!rtIsInf(gravity[1])) && (!rtIsNaN(gravity[1])));
        b_gravity[2] = ((!rtIsInf(gravity[2])) && (!rtIsNaN(gravity[2])));
        if (all(b_gravity)) {
          good = true;
        } else {
          good = false;
        }
      } else {
        good = false;
      }
    } else {
      good = false;
    }
  } else {
    good = false;
  }
  *modelValid = false;
  shift[0] = 0.0;
  shift[1] = 0.0;
  if (good) {
    __m128d r;
    __m128d r1;
    double C[4];
    double actuatorForce_N_idx_0;
    double l_idx_0;
    /* PACNOMINALPARAMETERS Initial nominal model; gains are offline starting
     * values. */
    /*  lambda: 1/s; dampingGain: N*s/m. No payload or adaptation. */
    armS_single_core(&measured[0], &measured[2], gravity, M, C, G);
    /* ARMS_SINGLE_BOUNDFORCE Soft actuator-length restoring force in q=[l2;l3].
     */
    /*  Default mode1: all l=[-l2-l3;l2;l3] use the same bounds, in metres. */
    /*  Mode0 reproduces historical coordinate-only l2/l3 bounds for old
     * analyses. */
    /*  coefficient*q equals force_N; coefficient is NOT the tangent stiffness.
     */
    l_idx_0 = -measured[0] - measured[1];
    actuatorForce_N_idx_0 = 500000.0 *
                            ((tanh(2000.0 * (l_idx_0 - 0.035)) + 2.0) -
                             tanh(2000.0 * (l_idx_0 - -0.035))) *
                            l_idx_0;
    /*  J=dl/dq=[-1 -1;1 0;0 1]; generalized force is J'*actuatorForce. */
    memset(&MassForce_N[0], 0, sizeof(double) << 1);
    r = _mm_loadu_pd(&M[0]);
    r1 = _mm_loadu_pd(&MassForce_N[0]);
    _mm_storeu_pd(
        &MassForce_N[0],
        _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(desiredAcceleration[0]))));
    r = _mm_loadu_pd(&M[2]);
    r1 = _mm_loadu_pd(&MassForce_N[0]);
    _mm_storeu_pd(
        &MassForce_N[0],
        _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(desiredAcceleration[1]))));
    C[0] += 40.0;
    C[3] += 40.0;
    memset(&CoriolisForce_N[0], 0, sizeof(double) << 1);
    CoriolisForce_N[0] += C[0] * measured[2];
    CoriolisForce_N[1] += C[1] * measured[2];
    GnStiff_N[0] = G[0] + (1300.0 * measured[0] + 650.0 * measured[1]);
    CoriolisForce_N[0] += C[2] * measured[3];
    CoriolisForce_N[1] += C[3] * measured[3];
    GnStiff_N[1] = G[1] + (650.0 * measured[0] + 1300.0 * measured[1]);
    l_idx_0 =
        GnStiff_N[0] + (500000.0 *
                            ((tanh(2000.0 * (measured[0] - 0.035)) + 2.0) -
                             tanh(2000.0 * (measured[0] - -0.035))) *
                            measured[0] -
                        actuatorForce_N_idx_0);
    GnStiff_N[0] = l_idx_0;
    l_idx_0 += MassForce_N[0] + CoriolisForce_N[0];
    tauDesired_N[0] = l_idx_0;
    bv1[0] = ((!rtIsInf(l_idx_0)) && (!rtIsNaN(l_idx_0)));
    l_idx_0 =
        GnStiff_N[1] + (500000.0 *
                            ((tanh(2000.0 * (measured[1] - 0.035)) + 2.0) -
                             tanh(2000.0 * (measured[1] - -0.035))) *
                            measured[1] -
                        actuatorForce_N_idx_0);
    GnStiff_N[1] = l_idx_0;
    l_idx_0 += MassForce_N[1] + CoriolisForce_N[1];
    tauDesired_N[1] = l_idx_0;
    bv1[1] = ((!rtIsInf(l_idx_0)) && (!rtIsNaN(l_idx_0)));
    if (d_all(bv1)) {
      bv[0] = ((!rtIsInf(M[0])) && (!rtIsNaN(M[0])));
      bv[1] = ((!rtIsInf(M[1])) && (!rtIsNaN(M[1])));
      bv[2] = ((!rtIsInf(M[2])) && (!rtIsNaN(M[2])));
      bv[3] = ((!rtIsInf(M[3])) && (!rtIsNaN(M[3])));
      if (c_all(bv) && (M[0] > 0.0) && (det(M) > 0.0)) {
        good = true;
      } else {
        good = false;
      }
    } else {
      good = false;
    }
    *modelValid = good;
    if (enableResidual && good) {
      /* RESIDUALGRAVITYSHIFT Re-express a fixed-vertical trained model, not new
       * learning. */
      /*  LHS(g)=LHS(g0)+G(g)-G(g0), hence r_g=r0+G(g)-G(g0). */
      /*  The paired model still predicts the original trained physical
       * configuration. */
      shift[0] = 0.0;
      shift[1] = 0.0;
      bv[0] = ((!rtIsInf(measured[0])) && (!rtIsNaN(measured[0])));
      bv[1] = ((!rtIsInf(measured[1])) && (!rtIsNaN(measured[1])));
      bv[2] = ((!rtIsInf(measured[2])) && (!rtIsNaN(measured[2])));
      bv[3] = ((!rtIsInf(measured[3])) && (!rtIsNaN(measured[3])));
      if (c_all(bv)) {
        b_gravity[0] = ((!rtIsInf(gravity[0])) && (!rtIsNaN(gravity[0])));
        b_gravity[1] = ((!rtIsInf(gravity[1])) && (!rtIsNaN(gravity[1])));
        b_gravity[2] = ((!rtIsInf(gravity[2])) && (!rtIsNaN(gravity[2])));
        if (all(b_gravity)) {
          b_gravity[0] = (gravity[0] == 0.0);
          b_gravity[1] = (gravity[1] == 0.0);
          b_gravity[2] = (gravity[2] == -9.81);
          if (!all(b_gravity)) {
            /* PACNOMINALPARAMETERS Initial nominal model; gains are offline
             * starting values. */
            /*  lambda: 1/s; dampingGain: N*s/m. No payload or adaptation. */
            armS_single_core(&measured[0], &measured[2], gravity, M, C, shift);
            armS_single_core(&measured[0], &measured[2], t2_g, M, C, G);
            r = _mm_loadu_pd(&shift[0]);
            r1 = _mm_loadu_pd(&G[0]);
            _mm_storeu_pd(&shift[0], _mm_sub_pd(r, r1));
          }
        }
      }
    }
  }
  if (!*modelValid) {
    tauDesired_N[0] = 0.0;
    MassForce_N[0] = 0.0;
    CoriolisForce_N[0] = 0.0;
    GnStiff_N[0] = 0.0;
    tauDesired_N[1] = 0.0;
    MassForce_N[1] = 0.0;
    CoriolisForce_N[1] = 0.0;
    GnStiff_N[1] = 0.0;
  }
  expl_temp.enableCommonAdaptation = enableCommonPressureAdaptation;
  expl_temp.maxAnchorSpeed_mPerS = supervisorOptions[5];
  expl_temp.maxSampleAge_s = supervisorOptions[4];
  expl_temp.invalidGrace_s = supervisorOptions[3];
  expl_temp.stopTime_s = supervisorOptions[6];
  expl_temp.blendTime_s = supervisorOptions[2];
  expl_temp.feedbackStart_s = supervisorOptions[1];
  expl_temp.pressureStart_s = supervisorOptions[0];
  expl_temp.forcePerBar_N = pressureOptions[9];
  expl_temp.activation_bar = pressureOptions[8];
  expl_temp.slew_barPerS = pressureOptions[7];
  expl_temp.maxDeltaA_bar = pressureOptions[6];
  expl_temp.maxDeltaS_bar = pressureOptions[5];
  expl_temp.maximum_bar = pressureOptions[4];
  expl_temp.minimum_bar = pressureOptions[3];
  expl_temp.aTrim_bar = pressureOptions[2];
  expl_temp.sTrim_bar = pressureOptions[1];
  expl_temp.common_bar = pressureOptions[0];
  *limited = residualPressureStep(
      tauDesired_N, GnStiff_N, measured, *modelValid, measurementTime_s, t,
      enablePressure, &expl_temp, &ps, &ns, &network, residualOptions,
      enableResidual, shift, pressure_bar, tauAvailable_N, G, status, &good,
      a__5, residualDiagnostic);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void pdfbResidualInverseBlock_init(void)
{
  static const b_struct_T r = {
      {0.0, 0.0, 0.0}, /* pressure */
      -1.0,            /* lastTime */
      -1.0,            /* lastGood */
      false,           /* anchorReady */
      {0.0, 0.0},      /* anchorForce */
      {0.0, 0.0},      /* forceBias */
      -1.0,            /* engageTime */
      false,           /* fault */
      {0.0, 0.0},      /* residualAnchor */
      {0.0, 0.0},      /* linearRequest */
      {0.0, 0.0}       /* lastEffective */
  };
  static const c_struct_T r1 = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, /* h */
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, /* pre */
      -1.0,                                                    /* lastSample */
      0.05,                                                    /* featureDt */
      0.0,                                                     /* warm_s */
      0.0,                                                     /* units */
      false,                                                   /* mode */
      false                                                    /* initialized */
  };
  ps = r;
  ns = r1;
}

/*
 * File trailer for pdfbResidualInverseBlock.c
 *
 * [EOF]
 */
