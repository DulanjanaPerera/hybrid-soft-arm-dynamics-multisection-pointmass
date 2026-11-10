/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacTuningConfig.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "pacTuningConfig.h"
#include "pacResidualAdaptiveTunableBlock_internal_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 * PACTUNINGCONFIG Map fixed-size external tuning vectors. See TUNABLE_BLOCK.md.
 *
 * Arguments    : const double b_gamma[10]
 *                const double scale[10]
 *                const double maxRate[10]
 *                const double lower[10]
 *                const double upper[10]
 *                const double initial[10]
 *                const double options[3]
 *                const double history[4]
 *                const double pressure[10]
 *                const double supervisor[7]
 *                g_struct_T *cfg
 *                struct_T *pc
 * Return Type  : void
 */
void pacTuningConfig(const double b_gamma[10], const double scale[10],
                     const double maxRate[10], const double lower[10],
                     const double upper[10], const double initial[10],
                     const double options[3], const double history[4],
                     const double pressure[10], const double supervisor[7],
                     g_struct_T *cfg, struct_T *pc)
{
  /*  Experimental signed stiffness compensation; bounds do not guarantee K SPD.
   */
  memcpy(&cfg->gamma[0], &b_gamma[0], 10U * sizeof(double));
  memcpy(&cfg->scale[0], &scale[0], 10U * sizeof(double));
  memcpy(&cfg->maxRate[0], &maxRate[0], 10U * sizeof(double));
  memcpy(&cfg->lower[0], &lower[0], 10U * sizeof(double));
  memcpy(&cfg->upper[0], &upper[0], 10U * sizeof(double));
  memcpy(&cfg->initial[0], &initial[0], 10U * sizeof(double));
  cfg->deadband_mPerS = options[0];
  cfg->residualTolerance_N = options[1];
  cfg->maxUpdateInterval_s = options[2];
  cfg->historyLength_m = history[0];
  cfg->historyBeta = history[1];
  cfg->historyGamma = history[2];
  cfg->historyBound = history[3];
  pc->common_bar = pressure[0];
  pc->sTrim_bar = pressure[1];
  pc->aTrim_bar = pressure[2];
  pc->minimum_bar = pressure[3];
  pc->maximum_bar = pressure[4];
  pc->maxDeltaS_bar = pressure[5];
  pc->maxDeltaA_bar = pressure[6];
  pc->slew_barPerS = pressure[7];
  pc->activation_bar = pressure[8];
  pc->forcePerBar_N = pressure[9];
  pc->pressureStart_s = supervisor[0];
  pc->feedbackStart_s = supervisor[1];
  pc->blendTime_s = supervisor[2];
  pc->invalidGrace_s = supervisor[3];
  pc->maxSampleAge_s = supervisor[4];
  pc->maxAnchorSpeed_mPerS = supervisor[5];
  pc->stopTime_s = supervisor[6];
}

/*
 * File trailer for pacTuningConfig.c
 *
 * [EOF]
 */
