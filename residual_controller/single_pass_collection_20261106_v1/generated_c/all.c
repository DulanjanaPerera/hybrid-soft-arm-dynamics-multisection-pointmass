/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: all.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 06-Nov-2026 18:50:26
 */

/* Include Files */
#include "all.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : const boolean_T x[140]
 * Return Type  : boolean_T
 */
boolean_T all(const boolean_T x[140])
{
  int k;
  boolean_T exitg1;
  boolean_T y;
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= 139)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }
  return y;
}

/*
 * Arguments    : const boolean_T x[16]
 * Return Type  : boolean_T
 */
boolean_T b_all(const boolean_T x[16])
{
  int k;
  boolean_T exitg1;
  boolean_T y;
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= 15)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }
  return y;
}

/*
 * File trailer for all.c
 *
 * [EOF]
 */
