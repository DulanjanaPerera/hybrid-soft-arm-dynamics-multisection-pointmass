/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacResidualAdaptiveTunableBlock_initialize.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "pacResidualAdaptiveTunableBlock_initialize.h"
#include "pacResidualAdaptiveTunableBlock.h"
#include "pacResidualAdaptiveTunableBlock_data.h"
#include "residualControllerNetwork.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void pacResidualAdaptiveTunableBlock_initialize(void)
{
  omp_init_nest_lock(&pacResidualAdaptiveTunableBlock_nestLockGlobal);
  c_pacResidualAdaptiveTunableBlo();
  residualControllerNetwork_init();
  isInitialized_pacResidualAdaptiveTunableBlock = true;
}

/*
 * File trailer for pacResidualAdaptiveTunableBlock_initialize.c
 *
 * [EOF]
 */
