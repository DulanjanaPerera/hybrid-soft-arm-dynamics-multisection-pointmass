/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacResidualAdaptiveTunableBlock.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "pacResidualAdaptiveTunableBlock.h"
#include "all.h"
#include "allOrAny.h"
#include "armS_single_core.h"
#include "det.h"
#include "diag.h"
#include "eye.h"
#include "mean.h"
#include "minOrMax.h"
#include "norm.h"
#include "pacResidualAdaptiveTunableBlock_data.h"
#include "pacResidualAdaptiveTunableBlock_initialize.h"
#include "pacResidualAdaptiveTunableBlock_internal_types.h"
#include "pacResidualAdaptiveTunableBlock_types.h"
#include "pacTuningConfig.h"
#include "residualPressureStep.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_d_struct_T
#define typedef_d_struct_T
typedef struct {
  double theta[10];
  double z[3];
  double lastSample;
} d_struct_T;
#endif /* typedef_d_struct_T */

/* Variable Definitions */
static b_struct_T ps;

static boolean_T ps_not_empty;

static c_struct_T ns;

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void c_pacResidualAdaptiveTunableBlo(void)
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
  ps_not_empty = false;
}

/*
 * Existing23 inputs/28 outputs preserved. Append enableResidual and8x1 options.
 *  The existing feedback/adaptive update stays unchanged. theta-theta0 now
 *  estimates remaining correction around the learned fixed-nominal model.
 *  Restart/reset all states before each model-off/on comparison.
 *
 * Arguments    : const double measured[6]
 *                const double desired[6]
 *                const double gravity_arm[3]
 *                boolean_T measurementValid
 *                boolean_T baseValid
 *                double measurementTime_s
 *                double t
 *                boolean_T enablePressure
 *                const double lambda[2]
 *                const double dampingGain[2]
 *                boolean_T enableAdaptation
 *                boolean_T enableCommonPressureAdaptation
 *                const double adaptationGamma[10]
 *                const double parameterScale[10]
 *                const double parameterMaxRate[10]
 *                const double parameterLower[10]
 *                const double parameterUpper[10]
 *                const double parameterInitial[10]
 *                const double adaptationOptions[3]
 *                const double historyOptions[4]
 *                const double pressureOptions[10]
 *                const double supervisorOptions[7]
 *                boolean_T freezeHysteresisWeights
 *                boolean_T enableResidual
 *                const double residualOptions[8]
 *                double pressureCommand_bar[3]
 *                double tauRequested_N[2]
 *                double tauEquivalent_N[2]
 *                double allocationResidual_N[2]
 *                double s_mPerS[2]
 *                double inertiaForce_N[2]
 *                double coriolisForce_N[2]
 *                double dampingForce_N[2]
 *                double staticForce_N[2]
 *                double feedbackForce_N[2]
 *                boolean_T *modelValid
 *                boolean_T *limited
 *                unsigned char *status
 *                boolean_T *allAboveActivation
 *                boolean_T *anchorReady
 *                double forceBias_N[2]
 *                double tauBlended_N[2]
 *                double stiffnessEstimate_NPerM[4]
 *                double dampingEstimate_NsPerM[4]
 *                double hysteresisForce_N[2]
 *                double adaptiveBias_N[2]
 *                double parameterEstimate[10]
 *                double historyState[3]
 *                boolean_T *adaptationActive
 *                boolean_T *parameterAtBound
 *                double adaptiveForce_N[2]
 *                double *commonPressure_bar
 *                double *commonPressureShift_bar
 *                struct0_T *residualDiagnostic
 * Return Type  : void
 */
void pacResidualAdaptiveTunableBlock(
    const double measured[6], const double desired[6],
    const double gravity_arm[3], boolean_T measurementValid,
    boolean_T baseValid, double measurementTime_s, double t,
    boolean_T enablePressure, const double lambda[2],
    const double dampingGain[2], boolean_T enableAdaptation,
    boolean_T enableCommonPressureAdaptation, const double adaptationGamma[10],
    const double parameterScale[10], const double parameterMaxRate[10],
    const double parameterLower[10], const double parameterUpper[10],
    const double parameterInitial[10], const double adaptationOptions[3],
    const double historyOptions[4], const double pressureOptions[10],
    const double supervisorOptions[7], boolean_T freezeHysteresisWeights,
    boolean_T enableResidual, const double residualOptions[8],
    double pressureCommand_bar[3], double tauRequested_N[2],
    double tauEquivalent_N[2], double allocationResidual_N[2],
    double s_mPerS[2], double inertiaForce_N[2], double coriolisForce_N[2],
    double dampingForce_N[2], double staticForce_N[2],
    double feedbackForce_N[2], boolean_T *modelValid, boolean_T *limited,
    unsigned char *status, boolean_T *allAboveActivation,
    boolean_T *anchorReady, double forceBias_N[2], double tauBlended_N[2],
    double stiffnessEstimate_NPerM[4], double dampingEstimate_NsPerM[4],
    double hysteresisForce_N[2], double adaptiveBias_N[2],
    double parameterEstimate[10], double historyState[3],
    boolean_T *adaptationActive, boolean_T *parameterAtBound,
    double adaptiveForce_N[2], double *commonPressure_bar,
    double *commonPressureShift_bar, struct0_T *residualDiagnostic)
{
  static d_struct_T as;
  static f_struct_T network;
  static const double t3_g[3] = {0.0, 0.0, -9.81};
  static const short iv[10] = {650, 650, 650, 40, 40, 0, 0, 0, 0, 0};
  static const signed char iv1[6] = {-1, -1, 1, 0, 0, 1};
  __m128d r;
  __m128d r1;
  __m128d r2;
  __m128d r3;
  e_struct_T expl_temp;
  g_struct_T cfg;
  struct_T pc;
  double Y[20];
  double C[4];
  double M[4];
  double actuatorForce_N[3];
  double w[3];
  double shift[2];
  double v[2];
  double cfg_historyBeta;
  double cfg_historyBound;
  double cfg_historyGamma;
  double d;
  double dt;
  double h;
  int b_i;
  int i;
  int n;
  unsigned char b_status;
  boolean_T bv3[10];
  boolean_T bv4[10];
  boolean_T lambda_data[10];
  boolean_T bv[6];
  boolean_T bv2[4];
  boolean_T b_gravity_arm[3];
  boolean_T bv1[2];
  boolean_T b_limited;
  boolean_T b_modelValid;
  boolean_T good;
  if (!isInitialized_pacResidualAdaptiveTunableBlock) {
    pacResidualAdaptiveTunableBlock_initialize();
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
  pacTuningConfig(adaptationGamma, parameterScale, parameterMaxRate,
                  parameterLower, parameterUpper, parameterInitial,
                  adaptationOptions, historyOptions, pressureOptions,
                  supervisorOptions, &cfg, &pc);
  cfg_historyBeta = cfg.historyBeta;
  cfg_historyGamma = cfg.historyGamma;
  cfg_historyBound = cfg.historyBound;
  if (freezeHysteresisWeights) {
    cfg.gamma[5] = 0.0;
    cfg.gamma[6] = 0.0;
    cfg.gamma[7] = 0.0;
  }
  if (!ps_not_empty) {
    ps_not_empty = true;
    as.z[0] = 0.0;
    as.z[1] = 0.0;
    as.z[2] = 0.0;
    as.lastSample = -1.0;
    memcpy(&as.theta[0], &cfg.initial[0], 10U * sizeof(double));
  }
  if (measurementValid && baseValid &&
      ((!rtIsInf(measurementTime_s)) && (!rtIsNaN(measurementTime_s))) &&
      ((!rtIsInf(t)) && (!rtIsNaN(t))) && (measurementTime_s <= t + 1.0E-9) &&
      (t - measurementTime_s <= pc.maxSampleAge_s + 1.0E-9)) {
    good = true;
  } else {
    good = false;
  }
  dt = 0.0;
  if (good && (measurementTime_s > as.lastSample)) {
    if (as.lastSample >= 0.0) {
      dt = measurementTime_s - as.lastSample;
    }
    for (i = 0; i < 6; i++) {
      h = measured[i];
      bv[i] = ((!rtIsInf(h)) && (!rtIsNaN(h)));
    }
    if (d_all(bv) && (!ps.fault) && (!(dt <= 0.0)) &&
        (!(dt > cfg.maxUpdateInterval_s))) {
      /* PACADAPTIVEHISTORY Bounded normalized Bouc-Wen n=1 feature states. */
      /*  Internal alpha/beta/gamma fixed; force weights adapt, not hidden-state
       * law. */
      bv1[0] = ((!rtIsInf(measured[2])) && (!rtIsNaN(measured[2])));
      bv1[1] = ((!rtIsInf(measured[3])) && (!rtIsNaN(measured[3])));
      if (c_all(bv1)) {
        w[0] = (-measured[2] - measured[3]) / cfg.historyLength_m;
        w[1] = measured[2] / cfg.historyLength_m;
        w[2] = measured[3] / cfg.historyLength_m;
        /*  Stable substeps for large measured velocity; avoid explicit-Euler
         * explosions. */
        actuatorForce_N[0] = fabs(w[0]);
        actuatorForce_N[1] = fabs(w[1]);
        actuatorForce_N[2] = fabs(w[2]);
        n = (int)fmin(200.0,
                      fmax(1.0, ceil(dt * maximum(actuatorForce_N) / 0.05)));
        h = dt / (double)n;
        d = -cfg.historyBound;
        r = _mm_loadu_pd(&actuatorForce_N[0]);
        r1 = _mm_loadu_pd(&w[0]);
        for (i = 0; i < n; i++) {
          shift[0] = fabs(as.z[0]);
          shift[1] = fabs(as.z[1]);
          r2 = _mm_loadu_pd(&as.z[0]);
          r3 = _mm_loadu_pd(&shift[0]);
          _mm_storeu_pd(
              &v[0],
              _mm_add_pd(
                  r2, _mm_mul_pd(
                          _mm_set1_pd(h),
                          _mm_sub_pd(
                              _mm_sub_pd(
                                  r1, _mm_mul_pd(
                                          _mm_mul_pd(
                                              _mm_set1_pd(cfg_historyBeta), r),
                                          r2)),
                              _mm_mul_pd(
                                  _mm_mul_pd(_mm_set1_pd(cfg_historyGamma), r1),
                                  r3)))));
          shift[0] = fmin(fmax(v[0], d), cfg_historyBound);
          shift[1] = fmin(fmax(v[1], d), cfg_historyBound);
          r2 = _mm_loadu_pd(&shift[0]);
          _mm_storeu_pd(&as.z[0], r2);
          as.z[2] = fmin(
              fmax(as.z[2] + h * ((w[2] - cfg_historyBeta * actuatorForce_N[2] *
                                              as.z[2]) -
                                  cfg_historyGamma * w[2] * fabs(as.z[2])),
                   -cfg_historyBound),
              cfg_historyBound);
        }
      }
    }
    as.lastSample = measurementTime_s;
  }
  /* PACADAPTIVEFORCE Nominal M/C/G, structured adaptive K/D/history/bias. */
  tauRequested_N[0] = 0.0;
  s_mPerS[0] = 0.0;
  tauRequested_N[1] = 0.0;
  s_mPerS[1] = 0.0;
  memset(&Y[0], 0, 20U * sizeof(double));
  inertiaForce_N[0] = 0.0;
  coriolisForce_N[0] = 0.0;
  dampingForce_N[0] = 0.0;
  staticForce_N[0] = 0.0;
  feedbackForce_N[0] = 0.0;
  hysteresisForce_N[0] = 0.0;
  adaptiveBias_N[0] = 0.0;
  inertiaForce_N[1] = 0.0;
  coriolisForce_N[1] = 0.0;
  dampingForce_N[1] = 0.0;
  staticForce_N[1] = 0.0;
  feedbackForce_N[1] = 0.0;
  hysteresisForce_N[1] = 0.0;
  adaptiveBias_N[1] = 0.0;
  stiffnessEstimate_NPerM[0] = 0.0;
  dampingEstimate_NsPerM[0] = 0.0;
  stiffnessEstimate_NPerM[1] = 0.0;
  dampingEstimate_NsPerM[1] = 0.0;
  stiffnessEstimate_NPerM[2] = 0.0;
  dampingEstimate_NsPerM[2] = 0.0;
  stiffnessEstimate_NPerM[3] = 0.0;
  dampingEstimate_NsPerM[3] = 0.0;
  b_modelValid = false;
  if (good) {
    for (i = 0; i < 6; i++) {
      h = measured[i];
      bv[i] = ((!rtIsInf(h)) && (!rtIsNaN(h)));
    }
    if (d_all(bv)) {
      for (i = 0; i < 6; i++) {
        h = desired[i];
        bv[i] = ((!rtIsInf(h)) && (!rtIsNaN(h)));
      }
      if (d_all(bv)) {
        b_gravity_arm[0] =
            ((!rtIsInf(gravity_arm[0])) && (!rtIsNaN(gravity_arm[0])));
        b_gravity_arm[1] =
            ((!rtIsInf(gravity_arm[1])) && (!rtIsNaN(gravity_arm[1])));
        b_gravity_arm[2] =
            ((!rtIsInf(gravity_arm[2])) && (!rtIsNaN(gravity_arm[2])));
        if (all(b_gravity_arm)) {
          bv1[0] = ((!rtIsInf(lambda[0])) && (!rtIsNaN(lambda[0])));
          bv1[1] = ((!rtIsInf(lambda[1])) && (!rtIsNaN(lambda[1])));
          if (c_all(bv1)) {
            bv1[0] = ((!rtIsInf(dampingGain[0])) && (!rtIsNaN(dampingGain[0])));
            bv1[1] = ((!rtIsInf(dampingGain[1])) && (!rtIsNaN(dampingGain[1])));
            if (c_all(bv1)) {
              lambda_data[0] = (lambda[0] <= 0.0);
              lambda_data[1] = (lambda[1] <= 0.0);
              if (!vectorAny(lambda_data, 2)) {
                lambda_data[0] = (dampingGain[0] <= 0.0);
                lambda_data[1] = (dampingGain[1] <= 0.0);
                if (!vectorAny(lambda_data, 2)) {
                  for (i = 0; i < 10; i++) {
                    lambda_data[i] =
                        ((!rtIsInf(as.theta[i])) && (!rtIsNaN(as.theta[i])));
                  }
                  if (b_all(lambda_data)) {
                    b_gravity_arm[0] =
                        ((!rtIsInf(as.z[0])) && (!rtIsNaN(as.z[0])));
                    b_gravity_arm[1] =
                        ((!rtIsInf(as.z[1])) && (!rtIsNaN(as.z[1])));
                    b_gravity_arm[2] =
                        ((!rtIsInf(as.z[2])) && (!rtIsNaN(as.z[2])));
                    if (all(b_gravity_arm)) {
                      /* PACNOMINALPARAMETERS Initial nominal model; gains are
                       * offline starting values. */
                      /*  lambda: 1/s; dampingGain: N*s/m. No payload or
                       * adaptation. */
                      r2 = _mm_sub_pd(
                          _mm_loadu_pd(&desired[2]),
                          _mm_mul_pd(_mm_loadu_pd(&lambda[0]),
                                     _mm_sub_pd(_mm_loadu_pd(&measured[0]),
                                                _mm_loadu_pd(&desired[0]))));
                      _mm_storeu_pd(&v[0], r2);
                      _mm_storeu_pd(&s_mPerS[0],
                                    _mm_sub_pd(_mm_loadu_pd(&measured[2]), r2));
                      armS_single_core(&measured[0], &measured[2], gravity_arm,
                                       M, C, shift);
                      bv2[0] = ((!rtIsInf(M[0])) && (!rtIsNaN(M[0])));
                      bv2[1] = ((!rtIsInf(M[1])) && (!rtIsNaN(M[1])));
                      bv2[2] = ((!rtIsInf(M[2])) && (!rtIsNaN(M[2])));
                      bv2[3] = ((!rtIsInf(M[3])) && (!rtIsNaN(M[3])));
                      if (e_all(bv2)) {
                        bv2[0] = ((!rtIsInf(C[0])) && (!rtIsNaN(C[0])));
                        bv2[1] = ((!rtIsInf(C[1])) && (!rtIsNaN(C[1])));
                        bv2[2] = ((!rtIsInf(C[2])) && (!rtIsNaN(C[2])));
                        bv2[3] = ((!rtIsInf(C[3])) && (!rtIsNaN(C[3])));
                        if (e_all(bv2)) {
                          bv1[0] =
                              ((!rtIsInf(shift[0])) && (!rtIsNaN(shift[0])));
                          bv1[1] =
                              ((!rtIsInf(shift[1])) && (!rtIsNaN(shift[1])));
                          if (c_all(bv1) && (!(M[0] <= 0.0)) &&
                              (!(det(M) <= 0.0))) {
                            double b_dv[9];
                            stiffnessEstimate_NPerM[0] =
                                as.theta[0] + as.theta[2];
                            stiffnessEstimate_NPerM[2] = as.theta[2];
                            stiffnessEstimate_NPerM[1] = as.theta[2];
                            stiffnessEstimate_NPerM[3] =
                                as.theta[1] + as.theta[2];
                            diag(&as.theta[3], dampingEstimate_NsPerM);
                            /*  Project each physical muscle history force into
                             * reduced length coordinates. */
                            Y[0] = measured[0];
                            Y[1] = 0.0;
                            Y[2] = 0.0;
                            Y[3] = measured[1];
                            h = measured[0] + measured[1];
                            Y[4] = h;
                            Y[5] = h;
                            Y[6] = v[0];
                            Y[7] = 0.0;
                            Y[8] = 0.0;
                            Y[9] = v[1];
                            b_diag(as.z, b_dv);
                            /* ARMS_SINGLE_BOUNDFORCE Soft actuator-length
                             * restoring force in q=[l2;l3]. */
                            /*  Default mode1: all l=[-l2-l3;l2;l3] use the same
                             * bounds, in metres. */
                            /*  Mode0 reproduces historical coordinate-only
                             * l2/l3 bounds for old analyses. */
                            /*  coefficient*q equals force_N; coefficient is NOT
                             * the tangent stiffness. */
                            w[0] = -measured[0] - measured[1];
                            w[1] = measured[0];
                            w[2] = measured[1];
                            for (i = 0; i < 3; i++) {
                              n = (i + 5) << 1;
                              h = 0.0;
                              d = 0.0;
                              for (b_i = 0; b_i < 3; b_i++) {
                                int Y_tmp;
                                cfg_historyBeta = b_dv[b_i + 3 * i];
                                Y_tmp = b_i << 1;
                                h += (double)iv1[Y_tmp] * cfg_historyBeta;
                                d += (double)iv1[Y_tmp + 1] * cfg_historyBeta;
                              }
                              Y[n + 1] = d;
                              Y[n] = h;
                              h = w[i];
                              actuatorForce_N[i] =
                                  500000.0 *
                                  ((tanh(2000.0 * (h - 0.035)) + 2.0) -
                                   tanh(2000.0 * (h - -0.035))) *
                                  h;
                            }
                            eye(&Y[16]);
                            /*  J=dl/dq=[-1 -1;1 0;0 1]; generalized force is
                             * J'*actuatorForce. */
                            memset(&inertiaForce_N[0], 0, sizeof(double) << 1);
                            memset(&coriolisForce_N[0], 0, sizeof(double) << 1);
                            memset(&dampingForce_N[0], 0, sizeof(double) << 1);
                            h = desired[4] -
                                lambda[0] * (measured[2] - desired[2]);
                            inertiaForce_N[0] += M[0] * h;
                            coriolisForce_N[0] += C[0] * v[0];
                            dampingForce_N[0] +=
                                dampingEstimate_NsPerM[0] * v[0];
                            inertiaForce_N[1] += M[1] * h;
                            coriolisForce_N[1] += v[0] * C[1];
                            dampingForce_N[1] +=
                                v[0] * dampingEstimate_NsPerM[1];
                            staticForce_N[0] =
                                shift[0] +
                                (stiffnessEstimate_NPerM[0] * measured[0] +
                                 measured[1] * as.theta[2]);
                            h = desired[5] -
                                lambda[1] * (measured[3] - desired[3]);
                            inertiaForce_N[0] += M[2] * h;
                            coriolisForce_N[0] += v[1] * C[2];
                            dampingForce_N[0] +=
                                v[1] * dampingEstimate_NsPerM[2];
                            inertiaForce_N[1] += M[3] * h;
                            coriolisForce_N[1] += v[1] * C[3];
                            dampingForce_N[1] +=
                                v[1] * dampingEstimate_NsPerM[3];
                            staticForce_N[1] =
                                shift[1] +
                                (measured[0] * as.theta[2] +
                                 measured[1] * stiffnessEstimate_NPerM[3]);
                            v[0] = actuatorForce_N[1] - actuatorForce_N[0];
                            v[1] = actuatorForce_N[2] - actuatorForce_N[0];
                            r2 = _mm_loadu_pd(&staticForce_N[0]);
                            r3 = _mm_loadu_pd(&v[0]);
                            _mm_storeu_pd(&staticForce_N[0],
                                          _mm_add_pd(r2, r3));
                            r2 = _mm_loadu_pd(&s_mPerS[0]);
                            _mm_storeu_pd(
                                &feedbackForce_N[0],
                                _mm_mul_pd(
                                    _mm_mul_pd(_mm_loadu_pd(&dampingGain[0]),
                                               _mm_set1_pd(-1.0)),
                                    r2));
                            memset(&hysteresisForce_N[0], 0,
                                   sizeof(double) << 1);
                            r2 = _mm_loadu_pd(&Y[10]);
                            r3 = _mm_loadu_pd(&hysteresisForce_N[0]);
                            _mm_storeu_pd(
                                &hysteresisForce_N[0],
                                _mm_add_pd(
                                    r3,
                                    _mm_mul_pd(r2, _mm_set1_pd(as.theta[5]))));
                            r2 = _mm_loadu_pd(&Y[12]);
                            r3 = _mm_loadu_pd(&hysteresisForce_N[0]);
                            _mm_storeu_pd(
                                &hysteresisForce_N[0],
                                _mm_add_pd(
                                    r3,
                                    _mm_mul_pd(r2, _mm_set1_pd(as.theta[6]))));
                            r2 = _mm_loadu_pd(&Y[14]);
                            r3 = _mm_loadu_pd(&hysteresisForce_N[0]);
                            _mm_storeu_pd(
                                &hysteresisForce_N[0],
                                _mm_add_pd(
                                    r3,
                                    _mm_mul_pd(r2, _mm_set1_pd(as.theta[7]))));
                            adaptiveBias_N[0] = as.theta[8];
                            h = (((((inertiaForce_N[0] + coriolisForce_N[0]) +
                                    dampingForce_N[0]) +
                                   staticForce_N[0]) +
                                  feedbackForce_N[0]) +
                                 hysteresisForce_N[0]) +
                                as.theta[8];
                            tauRequested_N[0] = h;
                            bv1[0] = ((!rtIsInf(h)) && (!rtIsNaN(h)));
                            adaptiveBias_N[1] = as.theta[9];
                            h = (((((inertiaForce_N[1] + coriolisForce_N[1]) +
                                    dampingForce_N[1]) +
                                   staticForce_N[1]) +
                                  feedbackForce_N[1]) +
                                 hysteresisForce_N[1]) +
                                as.theta[9];
                            tauRequested_N[1] = h;
                            bv1[1] = ((!rtIsInf(h)) && (!rtIsNaN(h)));
                            b_modelValid = c_all(bv1);
                            if (!b_modelValid) {
                              tauRequested_N[0] = 0.0;
                              tauRequested_N[1] = 0.0;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  *modelValid = b_modelValid;
  /*  Anchor includes the same residual/activation/gravity definition as
   * inversion. */
  shift[0] = 0.0;
  shift[1] = 0.0;
  if (enableResidual && b_modelValid) {
    /* RESIDUALGRAVITYSHIFT Re-express a fixed-vertical trained model, not new
     * learning. */
    /*  LHS(g)=LHS(g0)+G(g)-G(g0), hence r_g=r0+G(g)-G(g0). */
    /*  The paired model still predicts the original trained physical
     * configuration. */
    shift[0] = 0.0;
    shift[1] = 0.0;
    bv2[0] = ((!rtIsInf(measured[0])) && (!rtIsNaN(measured[0])));
    bv2[1] = ((!rtIsInf(measured[1])) && (!rtIsNaN(measured[1])));
    bv2[2] = ((!rtIsInf(measured[2])) && (!rtIsNaN(measured[2])));
    bv2[3] = ((!rtIsInf(measured[3])) && (!rtIsNaN(measured[3])));
    if (e_all(bv2)) {
      b_gravity_arm[0] =
          ((!rtIsInf(gravity_arm[0])) && (!rtIsNaN(gravity_arm[0])));
      b_gravity_arm[1] =
          ((!rtIsInf(gravity_arm[1])) && (!rtIsNaN(gravity_arm[1])));
      b_gravity_arm[2] =
          ((!rtIsInf(gravity_arm[2])) && (!rtIsNaN(gravity_arm[2])));
      if (all(b_gravity_arm)) {
        b_gravity_arm[0] = (gravity_arm[0] == 0.0);
        b_gravity_arm[1] = (gravity_arm[1] == 0.0);
        b_gravity_arm[2] = (gravity_arm[2] == -9.81);
        if (!all(b_gravity_arm)) {
          /* PACNOMINALPARAMETERS Initial nominal model; gains are offline
           * starting values. */
          /*  lambda: 1/s; dampingGain: N*s/m. No payload or adaptation. */
          armS_single_core(&measured[0], &measured[2], gravity_arm, M, C,
                           shift);
          armS_single_core(&measured[0], &measured[2], t3_g, M, C, v);
          r2 = _mm_loadu_pd(&shift[0]);
          r3 = _mm_loadu_pd(&v[0]);
          _mm_storeu_pd(&shift[0], _mm_sub_pd(r2, r3));
        }
      }
    }
  }
  expl_temp.enableCommonAdaptation = enableCommonPressureAdaptation;
  expl_temp.maxAnchorSpeed_mPerS = pc.maxAnchorSpeed_mPerS;
  expl_temp.maxSampleAge_s = pc.maxSampleAge_s;
  expl_temp.invalidGrace_s = pc.invalidGrace_s;
  expl_temp.stopTime_s = pc.stopTime_s;
  expl_temp.blendTime_s = pc.blendTime_s;
  expl_temp.feedbackStart_s = pc.feedbackStart_s;
  expl_temp.pressureStart_s = pc.pressureStart_s;
  expl_temp.forcePerBar_N = pc.forcePerBar_N;
  expl_temp.activation_bar = pc.activation_bar;
  expl_temp.slew_barPerS = pc.slew_barPerS;
  expl_temp.maxDeltaA_bar = pc.maxDeltaA_bar;
  expl_temp.maxDeltaS_bar = pc.maxDeltaS_bar;
  expl_temp.maximum_bar = pc.maximum_bar;
  expl_temp.minimum_bar = pc.minimum_bar;
  expl_temp.aTrim_bar = pc.aTrim_bar;
  expl_temp.sTrim_bar = pc.sTrim_bar;
  expl_temp.common_bar = pc.common_bar;
  r2 = _mm_loadu_pd(&staticForce_N[0]);
  r3 = _mm_loadu_pd(&hysteresisForce_N[0]);
  r = _mm_loadu_pd(&adaptiveBias_N[0]);
  _mm_storeu_pd(&v[0], _mm_add_pd(_mm_add_pd(r2, r3), r));
  b_limited = residualPressureStep(
      tauRequested_N, v, measured, b_modelValid, measurementTime_s, t,
      enablePressure, &expl_temp, &ps, &ns, &network, residualOptions,
      enableResidual, shift, pressureCommand_bar, tauEquivalent_N,
      allocationResidual_N, &b_status, allAboveActivation, tauBlended_N,
      residualDiagnostic);
  *limited = b_limited;
  *status = b_status;
  *commonPressure_bar = mean(pressureCommand_bar);
  *commonPressureShift_bar = *commonPressure_bar - pc.common_bar;
  *anchorReady = ps.anchorReady;
  forceBias_N[0] = ps.forceBias[0];
  forceBias_N[1] = ps.forceBias[1];
  historyState[0] = as.z[0];
  historyState[1] = as.z[1];
  historyState[2] = as.z[2];
  for (i = 0; i < 10; i++) {
    h = as.theta[i];
    parameterEstimate[i] = h;
    residualDiagnostic->parameterCorrection[i] = h - (double)iv[i];
  }
  /* PACNOMINALPARAMETERS Initial nominal model; gains are offline starting
   * values. */
  /*  lambda: 1/s; dampingGain: N*s/m. No payload or adaptation. */
  M[0] = stiffnessEstimate_NPerM[0] - 1300.0;
  C[0] = dampingEstimate_NsPerM[0] - 40.0;
  M[1] = stiffnessEstimate_NPerM[1] - 650.0;
  C[1] = dampingEstimate_NsPerM[1];
  M[2] = stiffnessEstimate_NPerM[1] - 650.0;
  C[2] = dampingEstimate_NsPerM[2];
  M[3] = stiffnessEstimate_NPerM[3] - 1300.0;
  C[3] = dampingEstimate_NsPerM[3] - 40.0;
  memset(&adaptiveForce_N[0], 0, sizeof(double) << 1);
  memset(&v[0], 0, sizeof(double) << 1);
  r2 = _mm_loadu_pd(&M[0]);
  r3 = _mm_loadu_pd(&adaptiveForce_N[0]);
  _mm_storeu_pd(&adaptiveForce_N[0],
                _mm_add_pd(r3, _mm_mul_pd(r2, _mm_set1_pd(measured[0]))));
  r2 = _mm_loadu_pd(&C[0]);
  r3 = _mm_loadu_pd(&v[0]);
  _mm_storeu_pd(
      &v[0],
      _mm_add_pd(r3, _mm_mul_pd(r2, _mm_set1_pd(desired[2] -
                                                lambda[0] * (measured[0] -
                                                             desired[0])))));
  r2 = _mm_loadu_pd(&M[2]);
  r3 = _mm_loadu_pd(&adaptiveForce_N[0]);
  _mm_storeu_pd(&adaptiveForce_N[0],
                _mm_add_pd(r3, _mm_mul_pd(r2, _mm_set1_pd(measured[1]))));
  r2 = _mm_loadu_pd(&C[2]);
  r3 = _mm_loadu_pd(&v[0]);
  _mm_storeu_pd(
      &v[0],
      _mm_add_pd(r3, _mm_mul_pd(r2, _mm_set1_pd(desired[3] -
                                                lambda[1] * (measured[1] -
                                                             desired[1])))));
  r2 = _mm_loadu_pd(&adaptiveForce_N[0]);
  r3 = _mm_loadu_pd(&v[0]);
  r = _mm_loadu_pd(&hysteresisForce_N[0]);
  r1 = _mm_loadu_pd(&adaptiveBias_N[0]);
  _mm_storeu_pd(&adaptiveForce_N[0],
                _mm_add_pd(_mm_add_pd(_mm_add_pd(r2, r3), r), r1));
  if (enableAdaptation && good && b_modelValid && (b_status == 3) &&
      (!b_limited)) {
    bv1[0] = ((!rtIsInf(allocationResidual_N[0])) &&
              (!rtIsNaN(allocationResidual_N[0])));
    bv1[1] = ((!rtIsInf(allocationResidual_N[1])) &&
              (!rtIsNaN(allocationResidual_N[1])));
    if (c_all(bv1) &&
        (b_norm(allocationResidual_N) <= cfg.residualTolerance_N)) {
      good = true;
    } else {
      good = false;
    }
  } else {
    good = false;
  }
  /* PACADAPTIVEUPDATE Box-projected, rate-bounded gradient. No hidden I/O. */
  /*  Outputs are estimates for the NEXT controller sample. */
  *adaptationActive = false;
  for (i = 0; i < 10; i++) {
    h = as.theta[i];
    b_modelValid = (h <= cfg.lower[i]);
    bv3[i] = b_modelValid;
    b_limited = (h >= cfg.upper[i]);
    bv4[i] = b_limited;
    lambda_data[i] = (b_modelValid || b_limited);
  }
  *parameterAtBound = vectorAny(lambda_data, 10);
  if (good && (!(dt <= 0.0)) && (!(dt > cfg.maxUpdateInterval_s))) {
    boolean_T bv5[20];
    for (i = 0; i < 20; i++) {
      h = Y[i];
      bv5[i] = ((!rtIsInf(h)) && (!rtIsNaN(h)));
    }
    if (i_all(bv5)) {
      bv1[0] = ((!rtIsInf(s_mPerS[0])) && (!rtIsNaN(s_mPerS[0])));
      bv1[1] = ((!rtIsInf(s_mPerS[1])) && (!rtIsNaN(s_mPerS[1])));
      if (c_all(bv1) && (!(b_norm(s_mPerS) <= cfg.deadband_mPerS))) {
        /*  zeta=theta/scale; Yscaled=Y*diag(scale), zeta_dot=-Gamma*Yscaled'*s.
         */
        cfg_historyBeta = s_mPerS[0];
        cfg_historyGamma = s_mPerS[1];
        for (i = 0; i < 10; i++) {
          n = i << 1;
          h = cfg.scale[i];
          d = cfg.maxRate[i];
          h = fmin(
              fmax(-(h * h) * cfg.gamma[i] *
                       (Y[n] * cfg_historyBeta + Y[n + 1] * cfg_historyGamma),
                   -d),
              d);
          if ((bv3[i] && (h < 0.0)) || (bv4[i] && (h > 0.0))) {
            h = 0.0;
          }
          as.theta[i] =
              fmin(fmax(as.theta[i] + dt * h, cfg.lower[i]), cfg.upper[i]);
          lambda_data[i] = (h != 0.0);
        }
        *adaptationActive = vectorAny(lambda_data, 10);
        for (i = 0; i < 10; i++) {
          h = as.theta[i];
          lambda_data[i] = ((h <= cfg.lower[i]) || (h >= cfg.upper[i]));
        }
        *parameterAtBound = vectorAny(lambda_data, 10);
      }
    }
  }
}

/*
 * File trailer for pacResidualAdaptiveTunableBlock.c
 *
 * [EOF]
 */
