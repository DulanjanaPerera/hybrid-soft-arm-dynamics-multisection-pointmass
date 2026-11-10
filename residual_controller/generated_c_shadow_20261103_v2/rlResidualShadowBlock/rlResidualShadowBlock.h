/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualShadowBlock.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

#ifndef RLRESIDUALSHADOWBLOCK_H
#define RLRESIDUALSHADOWBLOCK_H

/* Include Files */
#include "rlResidualShadowBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void rlResidualShadowBlock(
    const double lengths3_m[3], boolean_T orientationValid, boolean_T baseValid,
    double pollStart_s, double tipFrame, const double pressure_bar[3],
    double commandTime_sim_s, double t_sim_s, boolean_T reset,
    const double shadowOptions[11], double predictionLength_m[6],
    double predictionAtMeasurement_m[6], double retimedPrediction_m[6],
    struct0_T *shadowDiagnostic, struct1_T *inverseContext);

void rlResidualShadowBlock_init(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for rlResidualShadowBlock.h
 *
 * [EOF]
 */
