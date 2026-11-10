/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: det.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "det.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : const double x[4]
 * Return Type  : double
 */
double det(const double x[4])
{
  double x_idx_1;
  double x_idx_2;
  double x_idx_3;
  double y;
  int a;
  signed char ipiv_idx_0;
  y = x[0];
  x_idx_1 = x[1];
  x_idx_2 = x[2];
  x_idx_3 = x[3];
  ipiv_idx_0 = 1;
  a = 0;
  if (fabs(x[1]) > fabs(x[0])) {
    a = 1;
  }
  if (x[a] != 0.0) {
    if (a != 0) {
      ipiv_idx_0 = 2;
      y = x[1];
      x_idx_1 = x[0];
      x_idx_2 = x[3];
      x_idx_3 = x[2];
    }
    x_idx_1 /= y;
  }
  if (x_idx_2 != 0.0) {
    x_idx_3 += x_idx_1 * -x_idx_2;
  }
  y *= x_idx_3;
  if (ipiv_idx_0 > 1) {
    y = -y;
  }
  return y;
}

/*
 * File trailer for det.c
 *
 * [EOF]
 */
