/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualSinglePassShadowBlock.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 05-Nov-2026 20:37:30
 */

#ifndef RLRESIDUALSINGLEPASSSHADOWBLOCK_H
#define RLRESIDUALSINGLEPASSSHADOWBLOCK_H

/* Include Files */
#include "rlResidualSinglePassShadowBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void rlResidualSinglePassShadowBlock(
    const struct0_T *inverseContext, const double desiredAcceleration_mPerS2[2],
    const double shadowPressureOptions[10], boolean_T variableCommon,
    const double singlePassOptions[2],
    double c_singlePassPressurePrediction_[3],
    double nominalPressurePrediction_bar[3], struct1_T *singlePassDiagnostic);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for rlResidualSinglePassShadowBlock.h
 *
 * [EOF]
 */
