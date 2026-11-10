/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pacTuningConfig.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

#ifndef PACTUNINGCONFIG_H
#define PACTUNINGCONFIG_H

/* Include Files */
#include "pacResidualAdaptiveTunableBlock_internal_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void pacTuningConfig(const double b_gamma[10], const double scale[10],
                     const double maxRate[10], const double lower[10],
                     const double upper[10], const double initial[10],
                     const double options[3], const double history[4],
                     const double pressure[10], const double supervisor[7],
                     g_struct_T *cfg, struct_T *pc);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for pacTuningConfig.h
 *
 * [EOF]
 */
