/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacAllocatePressure.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 05-Nov-2026 20:37:30
 */

#ifndef PACALLOCATEPRESSURE_H
#define PACALLOCATEPRESSURE_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
boolean_T
pacAllocatePressure(const double requested[2], double cfg_common_bar,
                    double cfg_sTrim_bar, double cfg_aTrim_bar,
                    double cfg_minimum_bar, double cfg_maximum_bar,
                    double cfg_maxDeltaS_bar, double cfg_maxDeltaA_bar,
                    double cfg_activation_bar, double cfg_forcePerBar_N,
                    boolean_T cfg_enableCommonAdaptation, double target[3],
                    double equivalentForce[2], double residual[2]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for pacAllocatePressure.h
 *
 * [EOF]
 */
