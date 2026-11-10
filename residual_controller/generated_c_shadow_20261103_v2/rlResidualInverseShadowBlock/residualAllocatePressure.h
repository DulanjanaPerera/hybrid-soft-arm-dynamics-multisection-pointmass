/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualAllocatePressure.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:25
 */

#ifndef RESIDUALALLOCATEPRESSURE_H
#define RESIDUALALLOCATEPRESSURE_H

/* Include Files */
#include "rlResidualInverseShadowBlock_internal_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
boolean_T residualAllocatePressure(
    const double request[2], const double previous[3],
    const double lastRequest[2], const double lastEffective[2], double step,
    const double measured[4], double featureDt, const double h[128],
    const double network_input[3072], const double network_recurrent[49152],
    const double network_bias[384], const double network_output[256],
    const double network_skip[16], const double network_outputBias[2],
    const double network_mean[8], const double network_std[8],
    const b_struct_T *cfg, const double options[8], double pressure[3],
    double b_error[2], unsigned char *iterations);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for residualAllocatePressure.h
 *
 * [EOF]
 */
