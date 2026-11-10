/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: rlResidualSinglePassCollectionShadowBlock_initialize.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 06-Nov-2026 18:50:26
 */

/* Include Files */
#include "rlResidualSinglePassCollectionShadowBlock_initialize.h"
#include "residualControllerNetwork.h"
#include "rlResidualSinglePassCollectionShadowBlock_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void rlResidualSinglePassCollectionShadowBlock_initialize(void)
{
  residualControllerNetwork_init();
  isInitialized_rlResidualSinglePassCollectionShadowBlock = true;
}

/*
 * File trailer for rlResidualSinglePassCollectionShadowBlock_initialize.c
 *
 * [EOF]
 */
