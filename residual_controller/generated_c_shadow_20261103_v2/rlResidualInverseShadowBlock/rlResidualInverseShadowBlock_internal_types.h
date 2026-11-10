/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualInverseShadowBlock_internal_types.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:25
 */

#ifndef RLRESIDUALINVERSESHADOWBLOCK_INTERNAL_TYPES_H
#define RLRESIDUALINVERSESHADOWBLOCK_INTERNAL_TYPES_H

/* Include Files */
#include "rlResidualInverseShadowBlock_types.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct {
  double common_bar;
  double sTrim_bar;
  double aTrim_bar;
  double minimum_bar;
  double maximum_bar;
  double maxDeltaS_bar;
  double maxDeltaA_bar;
  double activation_bar;
  double forcePerBar_N;
  boolean_T enableCommonAdaptation;
} b_struct_T;
#endif /* typedef_b_struct_T */

#endif
/*
 * File trailer for rlResidualInverseShadowBlock_internal_types.h
 *
 * [EOF]
 */
