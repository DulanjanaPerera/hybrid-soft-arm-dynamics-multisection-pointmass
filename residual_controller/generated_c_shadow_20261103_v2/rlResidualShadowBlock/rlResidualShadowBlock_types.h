/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualShadowBlock_types.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

#ifndef RLRESIDUALSHADOWBLOCK_TYPES_H
#define RLRESIDUALSHADOWBLOCK_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  boolean_T valid;
  double measured[4];
  double historyPre[128];
  double pressure_bar[3];
  double featureDt_s;
  double horizon_s;
  double sampleTime_host_s;
  double units;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  unsigned char status;
  boolean_T accepted;
  boolean_T advanced;
  boolean_T ready;
  boolean_T historyReset;
  boolean_T predictionValid;
  boolean_T scoreValid;
  boolean_T retimedScoreValid;
  boolean_T inverseContextValid;
  double qMeasured_m[2];
  double velocity_mPerS[2];
  double rawForce_N[2];
  double errorAtMeasurement_m[6];
  double retimedError_m[6];
  double sampleTime_host_s;
  double predictionTarget_host_s;
  double scoreOrigin_host_s;
  double dtPhysical_s;
  double commandAge_sim_s;
  double horizonMismatch_s;
  double pressureChangeSinceOrigin_bar[3];
  double velocityWarm_s;
  double neuralWarm_s;
  double maximumTrainingZ;
  double units;
  unsigned int segment;
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct {
  double input[3072];
  double recurrent[49152];
  double bias[384];
  double output[256];
  double skip[16];
  double outputBias[2];
  double mean[8];
  double std[8];
  double units;
} b_struct_T;
#endif /* typedef_b_struct_T */

#endif
/*
 * File trailer for rlResidualShadowBlock_types.h
 *
 * [EOF]
 */
