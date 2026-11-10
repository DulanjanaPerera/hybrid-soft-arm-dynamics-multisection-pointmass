/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pdfbResidualInverseBlock_initialize.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

/* Include Files */
#include "pdfbResidualInverseBlock_initialize.h"
#include "pdfbResidualInverseBlock.h"
#include "pdfbResidualInverseBlock_data.h"
#include "residualControllerNetwork.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void pdfbResidualInverseBlock_initialize(void)
{
  omp_init_nest_lock(&pdfbResidualInverseBlock_nestLockGlobal);
  pdfbResidualInverseBlock_init();
  residualControllerNetwork_init();
  isInitialized_pdfbResidualInverseBlock = true;
}

/*
 * File trailer for pdfbResidualInverseBlock_initialize.c
 *
 * [EOF]
 */
