/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacResidualAdaptiveTunableBlock.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

#ifndef PACRESIDUALADAPTIVETUNABLEBLOCK_H
#define PACRESIDUALADAPTIVETUNABLEBLOCK_H

/* Include Files */
#include "pacResidualAdaptiveTunableBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void c_pacResidualAdaptiveTunableBlo(void);

extern void pacResidualAdaptiveTunableBlock(
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
    double *commonPressureShift_bar, struct0_T *residualDiagnostic);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for pacResidualAdaptiveTunableBlock.h
 *
 * [EOF]
 */
