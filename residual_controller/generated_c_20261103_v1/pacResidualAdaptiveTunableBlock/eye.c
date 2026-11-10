/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: eye.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:17:11
 */

/* Include Files */
#include "eye.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : double b_I[4]
 * Return Type  : void
 */
void eye(double b_I[4])
{
  b_I[1] = 0.0;
  b_I[2] = 0.0;
  b_I[0] = 1.0;
  b_I[3] = 1.0;
}

/*
 * File trailer for eye.c
 *
 * [EOF]
 */
