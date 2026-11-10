/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualSinglePassShadowBlock_initialize.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 05-Nov-2026 20:37:30
 */

/* Include Files */
#include "rlResidualSinglePassShadowBlock_initialize.h"
#include "residualControllerNetwork.h"
#include "rlResidualSinglePassShadowBlock_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void rlResidualSinglePassShadowBlock_initialize(void)
{
  residualControllerNetwork_init();
  isInitialized_rlResidualSinglePassShadowBlock = true;
}

/*
 * File trailer for rlResidualSinglePassShadowBlock_initialize.c
 *
 * [EOF]
 */
