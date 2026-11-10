/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mldivide.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

/* Include Files */
#include "mldivide.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : const double A[4]
 *                const double B[2]
 *                double Y[2]
 * Return Type  : void
 */
void mldivide(const double A[4], const double B[2], double Y[2])
{
  double Y_tmp;
  double a21;
  int r1;
  int r2;
  if (fabs(A[1]) > fabs(A[0])) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }
  a21 = A[r2] / A[r1];
  Y_tmp = A[r1 + 2];
  Y[1] = (B[r2] - B[r1] * a21) / (A[r2 + 2] - a21 * Y_tmp);
  Y[0] = (B[r1] - Y[1] * Y_tmp) / A[r1];
}

/*
 * File trailer for mldivide.c
 *
 * [EOF]
 */
