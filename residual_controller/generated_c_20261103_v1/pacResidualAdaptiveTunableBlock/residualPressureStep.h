/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualPressureStep.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

#ifndef RESIDUALPRESSURESTEP_H
#define RESIDUALPRESSURESTEP_H

/* Include Files */
#include "pacResidualAdaptiveTunableBlock_internal_types.h"
#include "pacResidualAdaptiveTunableBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
boolean_T residualPressureStep(
    const double requested[2], const double staticForce[2],
    const double measured[6], boolean_T valid, double sampleTime, double t,
    boolean_T enable, const e_struct_T *cfg, b_struct_T *state,
    c_struct_T *b_ns, const f_struct_T *network, const double options[8],
    boolean_T residualEnabled, const double shift[2], double pressure[3],
    double equivalent[2], double b_error[2], unsigned char *status,
    boolean_T *active, double effective[2], struct0_T *d);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for residualPressureStep.h
 *
 * [EOF]
 */
