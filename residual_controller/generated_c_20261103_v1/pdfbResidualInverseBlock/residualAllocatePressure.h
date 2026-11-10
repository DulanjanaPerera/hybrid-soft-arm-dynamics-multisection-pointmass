/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualAllocatePressure.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

#ifndef RESIDUALALLOCATEPRESSURE_H
#define RESIDUALALLOCATEPRESSURE_H

/* Include Files */
#include "pdfbResidualInverseBlock_internal_types.h"
#include "pdfbResidualInverseBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
boolean_T residualAllocatePressure(
    const double request[2], const double anchor[2],
    const double anchorResidual[2], const double bias[2],
    const double previous[3], const double lastRequest[2],
    const double lastEffective[2], double step, const double measured[6],
    double featureDt, const double h[128], const d_struct_T *network,
    const struct_T *cfg, const double shift[2], const double options[8],
    double pressure[3], double b_error[2], unsigned char *iterations,
    double linearRequest[2], double target[3]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for residualAllocatePressure.h
 *
 * [EOF]
 */
