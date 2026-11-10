/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacAllocateVariableCommon.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

#ifndef PACALLOCATEVARIABLECOMMON_H
#define PACALLOCATEVARIABLECOMMON_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
boolean_T pacAllocateVariableCommon(
    const double requested[2], const double bias[2], double cfg_common_bar,
    double cfg_sTrim_bar, double cfg_aTrim_bar, double cfg_minimum_bar,
    double cfg_maximum_bar, double cfg_maxDeltaS_bar, double cfg_maxDeltaA_bar,
    double cfg_activation_bar, double cfg_forcePerBar_N, double target[3],
    double equivalent[2], double residual[2]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for pacAllocateVariableCommon.h
 *
 * [EOF]
 */
