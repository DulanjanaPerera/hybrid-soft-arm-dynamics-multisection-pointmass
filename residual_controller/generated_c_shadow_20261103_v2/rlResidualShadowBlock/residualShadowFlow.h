/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualShadowFlow.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

#ifndef RESIDUALSHADOWFLOW_H
#define RESIDUALSHADOWFLOW_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void residualShadowFlow(const double q[2], const double v[2], const double p[3],
                        const double r[2], double dt, double Q[6]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for residualShadowFlow.h
 *
 * [EOF]
 */
