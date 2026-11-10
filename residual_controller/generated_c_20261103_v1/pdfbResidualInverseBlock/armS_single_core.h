/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: armS_single_core.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

#ifndef ARMS_SINGLE_CORE_H
#define ARMS_SINGLE_CORE_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void armS_single_core(const double q[2], const double dq[2],
                      const double params_g[3], double M[4], double C[4],
                      double G[2]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for armS_single_core.h
 *
 * [EOF]
 */
