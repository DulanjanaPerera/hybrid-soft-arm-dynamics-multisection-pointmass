/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualDesiredPressureBlock.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 09-Nov-2026 11:59:05
 */

#ifndef RLRESIDUALDESIREDPRESSUREBLOCK_H
#define RLRESIDUALDESIREDPRESSUREBLOCK_H

/* Include Files */
#include "rlResidualDesiredPressureBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void rlResidualDesiredPressureBlock(
    const struct0_T *inverseContext, const double desiredMotion[6],
    const double shadowPressureOptions[10], boolean_T variableCommon,
    const double singlePassOptions[2], double desiredPressurePrediction_bar[3],
    double nominalPressurePrediction_bar[3],
    struct1_T *desiredPressureDiagnostic);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for rlResidualDesiredPressureBlock.h
 *
 * [EOF]
 */
