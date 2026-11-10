/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: forceResidualStep.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

#ifndef FORCERESIDUALSTEP_H
#define FORCERESIDUALSTEP_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void stage(const double q[2], const double v[2], const double force[2],
           double h, double next[2], double newV[2]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for forceResidualStep.h
 *
 * [EOF]
 */
