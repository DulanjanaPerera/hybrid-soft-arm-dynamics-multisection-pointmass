/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualInverseShadowBlock_types.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:25
 */

#ifndef RLRESIDUALINVERSESHADOWBLOCK_TYPES_H
#define RLRESIDUALINVERSESHADOWBLOCK_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  boolean_T valid;
  double measured[4];
  double historyPre[128];
  double pressure_bar[3];
  double featureDt_s;
  double horizon_s;
  double sampleTime_host_s;
  double units;
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  boolean_T valid;
  boolean_T converged;
  boolean_T limited;
  boolean_T nominalLimited;
  unsigned char iterations;
  double targetForce_N[2];
  double predictedForce_N[2];
  double inverseError_N[2];
  double nominalInverseError_N[2];
  double pressureDifference_bar[3];
  double sampleTime_host_s;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  double input[3072];
  double recurrent[49152];
  double bias[384];
  double output[256];
  double skip[16];
  double outputBias[2];
  double mean[8];
  double std[8];
} struct_T;
#endif /* typedef_struct_T */

#endif
/*
 * File trailer for rlResidualInverseShadowBlock_types.h
 *
 * [EOF]
 */
