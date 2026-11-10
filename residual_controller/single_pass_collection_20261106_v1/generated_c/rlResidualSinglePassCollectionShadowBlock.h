/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualSinglePassCollectionShadowBlock.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 06-Nov-2026 18:50:26
 */

#ifndef RLRESIDUALSINGLEPASSCOLLECTIONSHADOWBLOCK_H
#define RLRESIDUALSINGLEPASSCOLLECTIONSHADOWBLOCK_H

/* Include Files */
#include "rlResidualSinglePassCollectionShadowBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void rlResidualSinglePassCollectionShadowBlock(
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
 * File trailer for rlResidualSinglePassCollectionShadowBlock.h
 *
 * [EOF]
 */
