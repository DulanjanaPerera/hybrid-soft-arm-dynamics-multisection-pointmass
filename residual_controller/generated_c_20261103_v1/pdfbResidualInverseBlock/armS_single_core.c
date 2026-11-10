/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: armS_single_core.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

/* Include Files */
#include "armS_single_core.h"
#include "integratedJacobianProduct_compact.h"
#include "pdfbResidualInverseBlock_rtwutil.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <string.h>

/* Function Definitions */
/*
 * One distributed section, translational soft-material kinetic energy only.
 *  q,dq = [l12;l13] and their rates; dm = mi*dXi.
 *
 * Arguments    : const double q[2]
 *                const double dq[2]
 *                const double params_g[3]
 *                double M[4]
 *                double C[4]
 *                double G[2]
 * Return Type  : void
 */
void armS_single_core(const double q[2], const double dq[2],
                      const double params_g[3], double M[4], double C[4],
                      double G[2])
{
  __m128d r;
  __m128d r1;
  double dM[8];
  double b_dv[6];
  double l[3];
  double b_c;
  double b_s_tmp;
  double c;
  double c_c;
  double d;
  double d1;
  double d10;
  double d11;
  double d12;
  double d13;
  double d14;
  double d15;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double d_c;
  double e_c;
  double f_c;
  double g_c;
  double h_c;
  double i_c;
  double j_c;
  double k_c;
  double l_c;
  double m_c;
  double n_c;
  double o_c;
  double p_c;
  double s_tmp;
  int j;
  l[0] = 0.0;
  l[1] = q[0];
  l[2] = q[1];
  /*  integration of; */
  /*  */
  /*    1) (mu) position [3x1] */
  /*    2) (mu_q) Position Jacobian (w.r.t. l2 and l3) [3x2] */
  /*    3) (mu_qq) Hessian of position (w.r.t. l2 and l3) [3x2x2] (two [3x2]s)
   */
  /*        The second derivative of l is in the 3rd dimension. The order is l2
   */
  /*        then l3. */
  /*   */
  /*  Inputs: */
  /*    l   : length change of the PMA [3x1] [constant] (m) */
  /*    L   : Length of the continuum arm [constant] (m) */
  /*    r   : radial offset of the PMA [constant] (m) */
  c = rt_powd_snf(q[0], 4.0);
  b_c = rt_powd_snf(q[0], 3.0);
  c_c = rt_powd_snf(q[1], 3.0);
  d_c = rt_powd_snf(q[1], 4.0);
  e_c = rt_powd_snf(q[0], 6.0);
  f_c = rt_powd_snf(q[0], 5.0);
  g_c = rt_powd_snf(q[1], 5.0);
  h_c = rt_powd_snf(q[1], 6.0);
  i_c = rt_powd_snf(q[0], 8.0);
  j_c = rt_powd_snf(q[0], 7.0);
  k_c = rt_powd_snf(q[1], 7.0);
  l_c = rt_powd_snf(q[1], 8.0);
  m_c = rt_powd_snf(q[0], 10.0);
  n_c = rt_powd_snf(q[0], 9.0);
  o_c = rt_powd_snf(q[1], 9.0);
  p_c = rt_powd_snf(q[1], 10.0);
  /*  Hessian for l2 */
  /*  Hessian for l3 */
  c_integratedJacobianProduct_com(l, M, dM);
  r = _mm_loadu_pd(&M[0]);
  r1 = _mm_set1_pd(0.1);
  _mm_storeu_pd(&M[0], _mm_mul_pd(r1, r));
  r = _mm_loadu_pd(&M[2]);
  _mm_storeu_pd(&M[2], _mm_mul_pd(r1, r));
  r = _mm_loadu_pd(&dM[0]);
  _mm_storeu_pd(&dM[0], _mm_mul_pd(r1, r));
  r = _mm_loadu_pd(&dM[2]);
  _mm_storeu_pd(&dM[2], _mm_mul_pd(r1, r));
  r = _mm_loadu_pd(&dM[4]);
  _mm_storeu_pd(&dM[4], _mm_mul_pd(r1, r));
  r = _mm_loadu_pd(&dM[6]);
  _mm_storeu_pd(&dM[6], _mm_mul_pd(r1, r));
  d = q[0] * q[0];
  s_tmp = 1152.0 * f_c;
  b_s_tmp = 6.9205499999999995 * c;
  d1 = 0.028671245459999994 * b_c;
  d2 = 7.5477053673449971E-5 * d;
  d3 = q[1] * q[1];
  d4 = -1.0714722539482958E-7 * q[0];
  d5 = (6.7904554093973235E-11 - 8.0360419046122186E-8 * d) +
       2.5159017891149992E-5 * c;
  d6 = 0.0038595907349999992 * e_c;
  d7 = 0.35591399999999995 * i_c;
  d8 = 22.0 * m_c;
  d9 = 80.0 * q[0];
  b_dv[0] =
      0.1 *
      (((((((((((-2.08932 * p_c - d9 * 0.17411 * o_c) +
                (0.0344267703 - 47.009699999999995 * d) * l_c) +
               (0.19278991367999995 * q[0] - 104.466 * b_c) * k_c) +
              0.17411 *
                  ((3.0845879999999997 * d - 0.0022054804199999993) -
                   960.0 * c) *
                  h_c) +
             0.17411 *
                 ((-0.0099246618899999986 * q[0] + 5.53644 * b_c) - s_tmp) *
                 g_c) +
            0.17411 *
                (((1.5095410734689995E-5 - 0.021503434094999996 * d) +
                  b_s_tmp) -
                 1050.0 * e_c) *
                d_c) +
           0.17411 *
               (((5.0318035782299983E-5 * q[0] - d1) +
                 6.1691759999999993 * f_c) -
                720.0 * j_c) *
               c_c) +
          0.17411 *
              ((((d2 - 5.357361269741479E-8) - 0.024811654724999993 * c) +
                3.8755079999999995 * e_c) -
               360.0 * i_c) *
              d3) +
         0.17411 *
             ((((d4 + 6.0381642938759981E-5 * b_c) -
                0.013232882519999998 * f_c) +
               1.58184 * j_c) -
              120.0 * n_c) *
             q[1]) +
        0.17411 * (((d5 - d6) + d7) - d8)) /
       1.7921603940369992E-21 / 2.95540245E+9);
  d10 = 80.0 * n_c;
  b_dv[1] =
      0.1 *
      (((((((((((-3.8304199999999997 * p_c - 120.0 * q[0] * 0.17411 * o_c) +
                (0.06196818654 - 62.679599999999994 * d) * l_c) +
               (0.2754141624 * q[0] - 125.35919999999999 * b_c) * k_c) +
              0.17411 *
                  ((3.8755079999999995 * d - 0.0038595907349999992) -
                   1050.0 * c) *
                  h_c) +
             0.17411 *
                 ((-0.013232882519999998 * q[0] + 6.1691759999999993 * b_c) -
                  s_tmp) *
                 g_c) +
            0.17411 *
                (((2.5159017891149992E-5 - 0.024811654724999993 * d) +
                  b_s_tmp) -
                 960.0 * e_c) *
                d_c) +
           0.17411 *
               (((6.0381642938759981E-5 * q[0] - d1) + 5.53644 * f_c) -
                600.0 * j_c) *
               c_c) +
          0.17411 *
              ((((d2 - 8.0360419046122186E-8) - 0.021503434094999996 * c) +
                3.0845879999999997 * e_c) -
               270.0 * i_c) *
              d3) +
         0.17411 *
             ((((d4 + 5.0318035782299983E-5 * b_c) -
                0.0099246618899999986 * f_c) +
               1.1072879999999998 * j_c) -
              d10) *
             q[1]) +
        0.17411 * (((((6.7904554093973235E-11 - 5.357361269741479E-8 * d) +
                      1.5095410734689995E-5 * c) -
                     0.0022054804199999993 * e_c) +
                    0.19773 * i_c) -
                   12.0 * m_c)) /
       1.7921603940369992E-21 / 2.95540245E+9);
  b_s_tmp = 1.5095410734689995E-5 * d * d3;
  d4 = 0.0022054804199999993 * b_c * c_c;
  d11 = 0.59319 * c * d_c;
  d12 = 72.0 * f_c * g_c;
  d13 = 240.0 * j_c * c_c;
  d14 = 40.0 * q[0] * o_c;
  d15 = 8.0 * p_c;
  b_dv[2] =
      0.1 *
      (-0.30156736610581719 *
       (((((((((((((((((((((((((((((((d5 + 2.0127214312919993E-5 * b_c * q[1]) +
                                     b_s_tmp) -
                                    1.0063607156459996E-5 * q[0] * c_c) -
                                   5.0318035782299981E-6 * d_c) -
                                  d6) -
                                 0.0066164412599999988 * f_c * q[1]) -
                                0.0082705515749999983 * c * d3) -
                               d4) +
                              0.0016541103149999997 * d * d_c) +
                             0.0033082206299999994 * q[0] * g_c) +
                            0.0011027402099999997 * h_c) +
                           d7) +
                          0.949104 * j_c * q[1]) +
                         1.6609319999999999 * e_c * d3) +
                        1.4236559999999998 * f_c * c_c) +
                       d11) -
                      0.474552 * b_c * g_c) -
                     0.71182799999999991 * d * h_c) -
                    0.474552 * q[0] * k_c) -
                   0.118638 * l_c) -
                  d8) -
                 d10 * q[1]) -
                180.0 * i_c * d3) -
               d13) -
              210.0 * e_c * d_c) -
             d12) +
            60.0 * c * h_c) +
           120.0 * b_c * k_c) +
          90.0 * d * l_c) +
         d14) +
        d15) /
       1.7921603940369992E-21 / 8.86620735E+9);
  s_tmp = 8.0 * m_c;
  d1 = 40.0 * n_c * q[1];
  d2 = 240.0 * b_c * k_c;
  b_dv[3] =
      0.1 *
      (0.30156736610581719 *
       (((((((((((((((((((((((((((((((((6.7904554093973235E-11 -
                                        8.0360419046122186E-8 * d3) -
                                       5.0318035782299981E-6 * c) -
                                      1.0063607156459996E-5 * b_c * q[1]) +
                                     b_s_tmp) +
                                    2.0127214312919993E-5 * q[0] * c_c) +
                                   2.5159017891149992E-5 * d_c) +
                                  0.0011027402099999997 * e_c) +
                                 0.0033082206299999994 * f_c * q[1]) +
                                0.0016541103149999997 * c * d3) -
                               d4) -
                              0.0082705515749999983 * d * d_c) -
                             0.0066164412599999988 * q[0] * g_c) -
                            0.0038595907349999992 * h_c) -
                           0.118638 * i_c) -
                          0.474552 * j_c * q[1]) -
                         0.71182799999999991 * e_c * d3) -
                        0.474552 * f_c * c_c) +
                       d11) +
                      1.4236559999999998 * b_c * g_c) +
                     1.6609319999999999 * d * h_c) +
                    0.949104 * q[0] * k_c) +
                   0.35591399999999995 * l_c) +
                  s_tmp) +
                 d1) +
                90.0 * i_c * d3) +
               120.0 * j_c * c_c) +
              60.0 * e_c * d_c) -
             d12) -
            210.0 * c * h_c) -
           d2) -
          180.0 * d * l_c) -
         d9 * o_c) -
        22.0 * p_c) /
       1.7921603940369992E-21 / 8.86620735E+9);
  d1 = ((((((((((((((((((((((((((((((((((1.5844395955260423E-10 -
                                         8.3336730862645234E-8 * d) -
                                        8.3336730862645234E-8 * q[0] * q[1]) -
                                       8.3336730862645234E-8 * d3) +
                                      1.7611312523804993E-5 * c) +
                                     3.5222625047609986E-5 * b_c * q[1]) +
                                    5.2833937571414983E-5 * d * d3) +
                                   3.5222625047609986E-5 * q[0] * c_c) +
                                  1.7611312523804993E-5 * d_c) -
                                 0.0020584483919999993 * e_c) -
                                0.0061753451759999989 * f_c * q[1]) -
                               0.012350690351999998 * c * d3) -
                              0.014409138743999997 * b_c * c_c) -
                             0.012350690351999998 * d * d_c) -
                            0.0061753451759999989 * q[0] * g_c) -
                           0.0020584483919999993 * h_c) +
                          0.15378999999999998 * i_c) +
                         0.61515999999999993 * j_c * q[1]) +
                        1.5378999999999998 * e_c * d3) +
                       2.4606399999999997 * f_c * c_c) +
                      2.9220099999999998 * c * d_c) +
                     2.4606399999999997 * b_c * g_c) +
                    1.5378999999999998 * d * h_c) +
                   0.61515999999999993 * q[0] * k_c) +
                  0.15378999999999998 * l_c) -
                 s_tmp) -
                d1) -
               120.0 * i_c * d3) -
              d13) -
             360.0 * e_c * d_c) -
            408.0 * f_c * g_c) -
           360.0 * c * h_c) -
          d2) -
         120.0 * d * l_c) -
        d14) -
       d15;
  b_dv[4] =
      0.1 * (-((2.0 * q[0] + q[1]) * 0.17411 * d1 / 2.3298085122480988E-23) /
             2.068781715E+10);
  b_dv[5] =
      0.1 * (-(0.17411 * (q[0] + 2.0 * q[1]) * d1 / 2.3298085122480988E-23) /
             2.068781715E+10);
  memset(&G[0], 0, sizeof(double) << 1);
  r = _mm_loadu_pd(&b_dv[0]);
  r1 = _mm_loadu_pd(&G[0]);
  _mm_storeu_pd(&G[0], _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(params_g[0]))));
  r = _mm_loadu_pd(&b_dv[2]);
  r1 = _mm_loadu_pd(&G[0]);
  _mm_storeu_pd(&G[0], _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(params_g[1]))));
  r = _mm_loadu_pd(&b_dv[4]);
  r1 = _mm_loadu_pd(&G[0]);
  _mm_storeu_pd(&G[0], _mm_add_pd(r1, _mm_mul_pd(r, _mm_set1_pd(params_g[2]))));
  /*  This function compute the C matrix using the christoffel symbol. The */
  /*  varible joint-space is given as dl. the dl is a vector concatinated of */
  /*  all the joint-space velocities. The input n represent the current section
   */
  /*  which is not used in here. */
  /*   */
  /*  The M is a 3D matrix consisting of partial derivatives of Mass matrix */
  /*  w.r.t. the all the join-space varibles upto section n. The 3rd dimension
   */
  /*  represents the joint-space variable. */
  /*     */
  /*        [[dMass], variable] */
  /*   */
  /*  Inputs: */
  /*    n   : (just for the clarity) current section [constant]  */
  /*    M   : 3D matrix with partial derivatives of M w.r.t. the joints-space */
  /*          varibles upto the section [2n, 2n, h] */
  /*    dl  : the respective joint-space varible velocity vector [2n,1] (m/s) */
  /*   */
  /*  Output: */
  /*    C   : current n section C matrix [2n,2n] */
  d1 = dq[0];
  d2 = dq[1];
  for (j = 0; j < 2; j++) {
    int c_s_tmp;
    c_s_tmp = j << 2;
    /*  + (M(k,j,h) + M(k,h,j) - M(h,j,k)) * dl(h); */
    s_tmp = dM[j + 4] + dM[j + 2];
    /*  + (M(k,j,h) + M(k,h,j) - M(h,j,k)) * dl(h); */
    C[j] = (((dM[j] + dM[j]) - dM[c_s_tmp]) * d1 +
            (s_tmp - dM[c_s_tmp + 2]) * d2) *
           0.5;
    /*  + (M(k,j,h) + M(k,h,j) - M(h,j,k)) * dl(h); */
    b_s_tmp = dM[j + 6];
    /*  + (M(k,j,h) + M(k,h,j) - M(h,j,k)) * dl(h); */
    C[j + 2] = ((s_tmp - dM[c_s_tmp + 1]) * d1 +
                ((b_s_tmp + b_s_tmp) - dM[c_s_tmp + 3]) * d2) *
               0.5;
  }
  /*  approximate Mdot via sum_h dM(:,:,h)*dq(h) */
  /*  this is for the sanity check */
  /*  Mdot = zeros(size(C)); */
  /*  for h=1:length(dl) */
  /*      Mdot = Mdot + M(:,:,h) * dl(h); */
  /*  end */
  /*   */
  /*  S = Mdot - 2*C; */
  /*  fprintf("Skew symmetry error: %g\n", norm(S + S.', 'fro')); */
}

/*
 * File trailer for armS_single_core.c
 *
 * [EOF]
 */
