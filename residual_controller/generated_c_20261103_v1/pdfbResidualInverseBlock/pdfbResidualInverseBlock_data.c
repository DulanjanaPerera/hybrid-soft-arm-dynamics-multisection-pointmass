/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pdfbResidualInverseBlock_data.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

/* Include Files */
#include "pdfbResidualInverseBlock_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
d_struct_T n64;

d_struct_T n128;

omp_nest_lock_t pdfbResidualInverseBlock_nestLockGlobal;

const double dv[6] = {
    -13.273228961416876, -13.273228961416876, 13.273228961416876, 0.0, 0.0,
    13.273228961416876};

boolean_T isInitialized_pdfbResidualInverseBlock = false;

/*
 * File trailer for pdfbResidualInverseBlock_data.c
 *
 * [EOF]
 */
