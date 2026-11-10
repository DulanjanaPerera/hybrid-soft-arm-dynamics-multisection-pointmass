/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pdfbResidualInverseBlock_terminate.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

/* Include Files */
#include "pdfbResidualInverseBlock_terminate.h"
#include "pdfbResidualInverseBlock_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void pdfbResidualInverseBlock_terminate(void)
{
  omp_destroy_nest_lock(&pdfbResidualInverseBlock_nestLockGlobal);
  isInitialized_pdfbResidualInverseBlock = false;
}

/*
 * File trailer for pdfbResidualInverseBlock_terminate.c
 *
 * [EOF]
 */
