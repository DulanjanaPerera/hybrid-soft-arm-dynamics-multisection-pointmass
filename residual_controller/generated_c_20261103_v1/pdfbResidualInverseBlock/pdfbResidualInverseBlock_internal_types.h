/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pdfbResidualInverseBlock_internal_types.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

#ifndef PDFBRESIDUALINVERSEBLOCK_INTERNAL_TYPES_H
#define PDFBRESIDUALINVERSEBLOCK_INTERNAL_TYPES_H

/* Include Files */
#include "pdfbResidualInverseBlock_types.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  double common_bar;
  double sTrim_bar;
  double aTrim_bar;
  double minimum_bar;
  double maximum_bar;
  double maxDeltaS_bar;
  double maxDeltaA_bar;
  double slew_barPerS;
  double activation_bar;
  double forcePerBar_N;
  double pressureStart_s;
  double feedbackStart_s;
  double blendTime_s;
  double stopTime_s;
  double invalidGrace_s;
  double maxSampleAge_s;
  double maxAnchorSpeed_mPerS;
  boolean_T enableCommonAdaptation;
} struct_T;
#endif /* typedef_struct_T */

#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct {
  double pressure[3];
  double lastTime;
  double lastGood;
  boolean_T anchorReady;
  double anchorForce[2];
  double forceBias[2];
  double engageTime;
  boolean_T fault;
  double residualAnchor[2];
  double linearRequest[2];
  double lastEffective[2];
} b_struct_T;
#endif /* typedef_b_struct_T */

#ifndef typedef_c_struct_T
#define typedef_c_struct_T
typedef struct {
  double h[128];
  double pre[128];
  double lastSample;
  double featureDt;
  double warm_s;
  double units;
  boolean_T mode;
  boolean_T initialized;
} c_struct_T;
#endif /* typedef_c_struct_T */

#endif
/*
 * File trailer for pdfbResidualInverseBlock_internal_types.h
 *
 * [EOF]
 */
