/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: forceResidualStep.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

/* Include Files */
#include "forceResidualStep.h"
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    double d;
    y = fabs(u0);
    d = fabs(u1);
    if (rtIsInf(u1)) {
      if (y == 1.0) {
        y = 1.0;
      } else if (y > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d == 0.0) {
      y = 1.0;
    } else if (d == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

/*
 * Arguments    : const double q[2]
 *                const double v[2]
 *                const double force[2]
 *                double h
 *                double next[2]
 *                double newV[2]
 * Return Type  : void
 */
void stage(const double q[2], const double v[2], const double force[2],
           double h, double next[2], double newV[2])
{
  static const double dv[11] = {0.05,
                                -0.0059523809523809521,
                                0.00034722222222222224,
                                -1.2776575276575277E-5,
                                3.2856803690137022E-7,
                                -6.26302709636043E-9,
                                9.1445458654842411E-11,
                                -1.024652859680243E-12,
                                8.54742507061871E-15,
                                -5.0026735363152986E-17,
                                1.743355929321998E-19};
  static const double dv1[11] = {
      0.003968253968253968,    -0.00030864197530864197, 1.2325637325637326E-5,
      -3.2503504725726946E-7,  6.2431444706576982E-9,   -1.0013286157077007E-10,
      2.648557146298413E-12,   -8.04166027354068E-14,   1.5203929502675726E-15,
      -1.7701767897731058E-17, 1.3754287410824441E-19};
  static const double dv2[7] = {0.5,
                                -0.041666666666666664,
                                0.0013888888888888889,
                                -2.48015873015873E-5,
                                2.7557319223985888E-7,
                                -2.08767569878681E-9,
                                1.1470745597729725E-11};
  double delta[2];
  double a11_tmp;
  double a12_tmp;
  double b_a11_tmp;
  double b_rhs1_tmp;
  double d;
  double d1;
  double j12;
  double rhs1_tmp;
  double s;
  int i;
  int outer;
  next[0] = q[0];
  newV[0] = v[0];
  next[1] = q[1];
  newV[1] = v[1];
  d = q[0];
  d1 = q[1];
  s = v[0];
  j12 = v[1];
  a11_tmp = h * h;
  b_a11_tmp = a11_tmp * 1300.0;
  a12_tmp = a11_tmp * 650.0;
  rhs1_tmp = d + h * s;
  b_rhs1_tmp = d1 + h * j12;
  for (outer = 0; outer < 3; outer++) {
    __m128d r;
    double a12;
    double a2;
    double a21;
    double a22;
    double a_tmp;
    double c1;
    double c_idx_0;
    double derivative_idx_1;
    double derivative_idx_2;
    double e1;
    double e2;
    double j11;
    double rhs1;
    double rhs2;
    double w1;
    double w2;
    /* FORCERESIDUALPHYSICS Batched differentiable EXACT nominal M,C,G
     * polynomials. */
    /*  q,v:2xbatch, numeric or unformatted dlarray. P rows M11,M12,M22, */
    /*  C11,C12,C21,C22,G1,G2. Fixed nominal basis, no fitted approximation. */
    /*  Generated rational coefficients; see generate_force_polynomials.py. */
    w1 = 7889.5463510848131 * next[0] + 3944.7731755424065 * next[1];
    w2 = 3944.7731755424065 * next[0] + 7889.5463510848131 * next[1];
    s = next[0] * w1 + next[1] * w2;
    a_tmp = 0.0 * s;
    rhs2 = a_tmp + 1.743355929321998E-19;
    j12 = a_tmp;
    a22 = a_tmp + 1.3754287410824441E-19;
    c_idx_0 = a_tmp;
    for (i = 0; i < 10; i++) {
      j12 = j12 * s + rhs2;
      rhs2 = rhs2 * s + dv[9 - i];
      c_idx_0 = c_idx_0 * s + a22;
      a22 = a22 * s + dv1[9 - i];
    }
    derivative_idx_2 = a_tmp + 1.1470745597729725E-11;
    for (i = 0; i < 6; i++) {
      a_tmp = a_tmp * s + derivative_idx_2;
      derivative_idx_2 = derivative_idx_2 * s + dv2[5 - i];
    }
    double d2;
    double d3;
    double d4;
    s = 2.0 * j12;
    a12 = 2.0 * c_idx_0;
    derivative_idx_2 = 2.0 * a22;
    rhs1 = w1 * w1;
    j12 = derivative_idx_2 * w1;
    e2 = s * w2;
    a21 = e2 * 7889.5463510848131;
    a2 = 0.0030314292099999996 *
         ((a21 + a12 * w2 * rhs1) + j12 * 3944.7731755424065);
    c_idx_0 = w2 * w2;
    e1 = s * w1;
    derivative_idx_1 = e1 * 7889.5463510848131;
    j11 = a12 * w1 * c_idx_0;
    derivative_idx_2 *= w2;
    c1 = 0.0030314292099999996 *
         ((derivative_idx_1 + j11) + derivative_idx_2 * 3944.7731755424065);
    s = rhs2 * 7889.5463510848131;
    d2 = 0.0030314292099999996 * (s + a22 * rhs1);
    d3 = 0.0030314292099999996 * (rhs2 * 3944.7731755424065 + a22 * w1 * w2);
    d4 = 0.0030314292099999996 * (s + a22 * c_idx_0);
    rhs2 = 0.5 * (0.0030314292099999996 *
                      ((derivative_idx_1 + a12 * rt_powd_snf(w1, 3.0)) +
                       j12 * 7889.5463510848131) *
                      newV[0] +
                  a2 * newV[1]);
    derivative_idx_1 =
        0.5 *
        (a2 * newV[0] +
         (2.0 * (0.0030314292099999996 *
                 ((e2 * 3944.7731755424065 + j11) +
                  a22 * (3944.7731755424065 * w2 + 7889.5463510848131 * w1))) -
          c1) *
             newV[1]);
    c_idx_0 =
        0.5 *
        ((2.0 * (0.0030314292099999996 *
                 ((e1 * 3944.7731755424065 + a12 * rhs1 * w2) +
                  a22 * (7889.5463510848131 * w2 + 3944.7731755424065 * w1))) -
          a2) *
             newV[0] +
         c1 * newV[1]);
    j12 = 0.5 * (c1 * newV[0] + 0.0030314292099999996 *
                                    ((a21 + a12 * rt_powd_snf(w2, 3.0)) +
                                     derivative_idx_2 * 7889.5463510848131) *
                                    newV[1]);
    a2 = (d2 + h * (rhs2 + 40.0)) + b_a11_tmp;
    a12 = (d3 + h * derivative_idx_1) + a12_tmp;
    a21 = (d3 + h * c_idx_0) + a12_tmp;
    a22 = (d4 + h * (j12 + 40.0)) + b_a11_tmp;
    s = -0.34160382 * a_tmp;
    rhs1 = ((d2 * rhs1_tmp + d3 * b_rhs1_tmp) +
            h * ((rhs2 + 40.0) * d + d1 * derivative_idx_1)) +
           a11_tmp * (force[0] - s * w1);
    rhs2 = ((d3 * rhs1_tmp + d4 * b_rhs1_tmp) +
            h * (d * c_idx_0 + (j12 + 40.0) * d1)) +
           a11_tmp * (force[1] - s * w2);
    for (i = 0; i < 5; i++) {
      __m128d r1;
      double varargin_1[2];
      s = -(next[0] + next[1]);
      j12 = tanh(2000.0 * (s - 0.035));
      derivative_idx_2 = tanh(2000.0 * (s - -0.035));
      j11 = 500000.0 * ((j12 + 2.0) - derivative_idx_2);
      c1 = j11 + 1.0E+9 * s * (derivative_idx_2 * derivative_idx_2 - j12 * j12);
      j11 *= s;
      c_idx_0 = j11;
      j12 = tanh(2000.0 * (next[0] - 0.035));
      derivative_idx_2 = tanh(2000.0 * (next[0] - -0.035));
      j11 = 500000.0 * ((j12 + 2.0) - derivative_idx_2);
      derivative_idx_1 =
          j11 +
          1.0E+9 * next[0] * (derivative_idx_2 * derivative_idx_2 - j12 * j12);
      j11 *= next[0];
      s = j11;
      j12 = tanh(2000.0 * (next[1] - 0.035));
      derivative_idx_2 = tanh(2000.0 * (next[1] - -0.035));
      j11 = 500000.0 * ((j12 + 2.0) - derivative_idx_2);
      derivative_idx_2 =
          j11 +
          1.0E+9 * next[1] * (derivative_idx_2 * derivative_idx_2 - j12 * j12);
      j11 *= next[1];
      e1 = ((a2 * next[0] + a12 * next[1]) + a11_tmp * (s - c_idx_0)) - rhs1;
      e2 = ((a21 * next[0] + a22 * next[1]) + a11_tmp * (j11 - c_idx_0)) - rhs2;
      j11 = a2 + a11_tmp * (c1 + derivative_idx_1);
      s = a11_tmp * c1;
      j12 = a12 + s;
      c_idx_0 = a21 + s;
      s = a22 + a11_tmp * (c1 + derivative_idx_2);
      derivative_idx_2 = j11 * s - j12 * c_idx_0;
      j12 = (s * e1 - j12 * e2) / derivative_idx_2;
      delta[0] = j12;
      s = (-c_idx_0 * e1 + j11 * e2) / derivative_idx_2;
      delta[1] = s;
      varargin_1[0] = fabs(j12);
      varargin_1[1] = fabs(s);
      r = _mm_loadu_pd(&delta[0]);
      r1 = _mm_loadu_pd(&next[0]);
      _mm_storeu_pd(
          &next[0],
          _mm_sub_pd(r1, _mm_mul_pd(r, _mm_set1_pd(fmin(
                                           1.0, 0.004 / (maximum(varargin_1) +
                                                         1.0E-30))))));
    }
    r = _mm_loadu_pd(&next[0]);
    _mm_storeu_pd(&newV[0], _mm_div_pd(_mm_sub_pd(r, _mm_loadu_pd(&q[0])),
                                       _mm_set1_pd(h)));
  }
}

/*
 * File trailer for forceResidualStep.c
 *
 * [EOF]
 */
