/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualDesiredPressureBlock_types.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 09-Nov-2026 11:59:05
 */

#ifndef RLRESIDUALDESIREDPRESSUREBLOCK_TYPES_H
#define RLRESIDUALDESIREDPRESSUREBLOCK_TYPES_H

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

#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  boolean_T valid;
  unsigned char networkEvaluations;
  boolean_T allocationLimited;
  boolean_T slewLimited;
  boolean_T nominalAllocationLimited;
  boolean_T nominalSlewLimited;
  boolean_T frozenBalanceSatisfied;
  double targetForce_N[2];
  double residualForce_N[2];
  double correctedForceRequest_N[2];
  double frozenTotalForce_N[2];
  double allocationError_N[2];
  double slewError_N[2];
  double frozenBalanceError_N[2];
  double nominalBalanceError_N[2];
  double residualInputPressure_bar[3];
  double correctedPressureBeforeSlew_bar[3];
  double nominalPressureBeforeSlew_bar[3];
  double c_pressureDifferenceBeforeSlew_[3];
  double pressureDifference_bar[3];
  double sampleTime_host_s;
  double horizon_s;
} struct1_T;
#endif /* typedef_struct1_T */

#endif
/*
 * File trailer for rlResidualDesiredPressureBlock_types.h
 *
 * [EOF]
 */
