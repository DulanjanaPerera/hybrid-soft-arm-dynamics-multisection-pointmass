/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacResidualAdaptiveTunableBlock_terminate.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "pacResidualAdaptiveTunableBlock_terminate.h"
#include "pacResidualAdaptiveTunableBlock_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void pacResidualAdaptiveTunableBlock_terminate(void)
{
  omp_destroy_nest_lock(&pacResidualAdaptiveTunableBlock_nestLockGlobal);
  isInitialized_pacResidualAdaptiveTunableBlock = false;
}

/*
 * File trailer for pacResidualAdaptiveTunableBlock_terminate.c
 *
 * [EOF]
 */
