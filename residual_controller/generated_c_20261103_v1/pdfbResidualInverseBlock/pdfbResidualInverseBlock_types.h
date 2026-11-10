/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pdfbResidualInverseBlock_types.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

#ifndef PDFBRESIDUALINVERSEBLOCK_TYPES_H
#define PDFBRESIDUALINVERSEBLOCK_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  boolean_T enabled;
  boolean_T ready;
  boolean_T advanced;
  boolean_T historyGapHeld;
  boolean_T contractValid;
  double rawForce_N[2];
  double effectiveResidual_N[2];
  double gravityReexpression_N[2];
  double inverseError_N[2];
  unsigned char iterations;
  boolean_T converged;
  double warm_s;
  double units;
  double parameterCorrection[10];
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_d_struct_T
#define typedef_d_struct_T
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
} d_struct_T;
#endif /* typedef_d_struct_T */

#endif
/*
 * File trailer for pdfbResidualInverseBlock_types.h
 *
 * [EOF]
 */
