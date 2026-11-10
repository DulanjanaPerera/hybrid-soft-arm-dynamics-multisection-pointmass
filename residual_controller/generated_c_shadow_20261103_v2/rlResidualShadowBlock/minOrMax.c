/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: minOrMax.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

/* Include Files */
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : const double x[8]
 * Return Type  : double
 */
double b_maximum(const double x[8])
{
  double ex;
  int b_k;
  int idx;
  if (!rtIsNaN(x[0])) {
    idx = 1;
  } else {
    int k;
    boolean_T exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 9)) {
      if (!rtIsNaN(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }
  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    idx++;
    for (b_k = idx; b_k < 9; b_k++) {
      double d;
      d = x[b_k - 1];
      if (ex < d) {
        ex = d;
      }
    }
  }
  return ex;
}

/*
 * Arguments    : const double x[2]
 * Return Type  : double
 */
double maximum(const double x[2])
{
  double ex;
  if ((x[0] < x[1]) || (rtIsNaN(x[0]) && (!rtIsNaN(x[1])))) {
    ex = x[1];
  } else {
    ex = x[0];
  }
  return ex;
}

/*
 * File trailer for minOrMax.c
 *
 * [EOF]
 */
