/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualShadowBlock_initialize.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

/* Include Files */
#include "rlResidualShadowBlock_initialize.h"
#include "residualControllerNetwork.h"
#include "rlResidualShadowBlock.h"
#include "rlResidualShadowBlock_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void rlResidualShadowBlock_initialize(void)
{
  rlResidualShadowBlock_init();
  residualControllerNetwork_init();
  isInitialized_rlResidualShadowBlock = true;
}

/*
 * File trailer for rlResidualShadowBlock_initialize.c
 *
 * [EOF]
 */
