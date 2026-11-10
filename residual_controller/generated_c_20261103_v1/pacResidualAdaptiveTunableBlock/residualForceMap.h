/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualForceMap.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

#ifndef RESIDUALFORCEMAP_H
#define RESIDUALFORCEMAP_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void residualForceMap(
    const double p[3], const double measured[6], double dt, const double h[128],
    const double network_input[3072], const double network_recurrent[49152],
    const double network_bias[384], const double network_output[256],
    const double network_skip[16], const double network_outputBias[2],
    const double network_mean[8], const double network_std[8],
    double cfg_activation_bar, double cfg_forcePerBar_N,
    const double gravityShift[2], double total[2], double correction[2]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for residualForceMap.h
 *
 * [EOF]
 */
