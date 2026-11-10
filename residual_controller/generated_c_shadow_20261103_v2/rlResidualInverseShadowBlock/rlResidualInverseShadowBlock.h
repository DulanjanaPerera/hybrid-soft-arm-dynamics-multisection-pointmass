/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualInverseShadowBlock.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:25
 */

#ifndef RLRESIDUALINVERSESHADOWBLOCK_H
#define RLRESIDUALINVERSESHADOWBLOCK_H

/* Include Files */
#include "rlResidualInverseShadowBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void rlResidualInverseShadowBlock(
    const struct0_T *inverseContext, const double desiredAcceleration_mPerS2[2],
    const double shadowPressureOptions[10], boolean_T variableCommon,
    const double shadowInverseOptions[8],
    double residualPressurePrediction_bar[3],
    double nominalPressurePrediction_bar[3], struct1_T *inverseDiagnostic);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for rlResidualInverseShadowBlock.h
 *
 * [EOF]
 */
