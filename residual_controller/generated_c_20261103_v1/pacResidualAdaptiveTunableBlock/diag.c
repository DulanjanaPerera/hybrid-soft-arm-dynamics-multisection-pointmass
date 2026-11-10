/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: diag.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "diag.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const double v[3]
 *                double d[9]
 * Return Type  : void
 */
void b_diag(const double v[3], double d[9])
{
  memset(&d[0], 0, 9U * sizeof(double));
  d[0] = v[0];
  d[4] = v[1];
  d[8] = v[2];
}

/*
 * Arguments    : const double v[2]
 *                double d[4]
 * Return Type  : void
 */
void diag(const double v[2], double d[4])
{
  d[1] = 0.0;
  d[2] = 0.0;
  d[0] = v[0];
  d[3] = v[1];
}

/*
 * File trailer for diag.c
 *
 * [EOF]
 */
