/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: residualShadowFlow.c
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 20:44:13
 */

/* Include Files */
#include "residualShadowFlow.h"
#include "forceResidualStep.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <string.h>

/* Function Definitions */
/*
 * RESIDUALSHADOWFLOW Columns: augmented, nominal, constant velocity.
 *  Fixed vertical-gravity training model. Held p/r over the prediction
 * interval. Uses the original trained nonlinear SDIRK propagator, not Euler
 * integration.
 *
 * Arguments    : const double q[2]
 *                const double v[2]
 *                const double p[3]
 *                const double r[2]
 *                double dt
 *                double Q[6]
 * Return Type  : void
 */
void residualShadowFlow(const double q[2], const double v[2], const double p[3],
                        const double r[2], double dt, double Q[6])
{
  static const double a[6] = {
      -13.273228961416876, -13.273228961416876, 13.273228961416876, 0.0, 0.0,
      13.273228961416876};
  __m128d b_r;
  __m128d r1;
  __m128d r2;
  double b_force[2];
  double b_q[2];
  double b_v[2];
  double force[2];
  double q1[2];
  double v1[2];
  double d;
  /* FORCERESIDUALSTEP Second-order L-stable SDIRK; nonlinear implicit stages.
   */
  /*  M(q)ddq+(C(q,v)+D)v+Kq+bound(q)=alpha*A*p+residual (N). */
  /*  Two-state batch solve/Newton is differentiable; no learned position jumps.
   */
  b_v[0] = v[0];
  b_q[0] = q[0];
  b_v[1] = v[1];
  b_q[1] = q[1];
  memset(&force[0], 0, sizeof(double) << 1);
  b_r = _mm_loadu_pd(&force[0]);
  _mm_storeu_pd(&force[0], _mm_add_pd(b_r, _mm_mul_pd(_mm_loadu_pd(&a[0]),
                                                      _mm_set1_pd(p[0]))));
  b_r = _mm_loadu_pd(&force[0]);
  _mm_storeu_pd(&force[0], _mm_add_pd(b_r, _mm_mul_pd(_mm_loadu_pd(&a[2]),
                                                      _mm_set1_pd(p[1]))));
  b_r = _mm_loadu_pd(&force[0]);
  _mm_storeu_pd(&force[0], _mm_add_pd(b_r, _mm_mul_pd(_mm_loadu_pd(&a[4]),
                                                      _mm_set1_pd(p[2]))));
  b_r = _mm_loadu_pd(&force[0]);
  _mm_storeu_pd(&b_force[0], _mm_add_pd(b_r, _mm_loadu_pd(&r[0])));
  d = 0.29289321881345254 * (dt / 4.0);
  stage(b_q, b_v, b_force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  r2 = _mm_set1_pd(2.4142135623730945);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, b_force, d, b_q, b_v);
  stage(b_q, b_v, b_force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, b_force, d, b_q, b_v);
  stage(b_q, b_v, b_force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, b_force, d, b_q, b_v);
  stage(b_q, b_v, b_force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, b_force, d, b_q, b_v);
  /* FORCERESIDUALSTEP Second-order L-stable SDIRK; nonlinear implicit stages.
   */
  /*  M(q)ddq+(C(q,v)+D)v+Kq+bound(q)=alpha*A*p+residual (N). */
  /*  Two-state batch solve/Newton is differentiable; no learned position jumps.
   */
  Q[0] = b_q[0];
  b_v[0] = v[0];
  b_q[0] = q[0];
  Q[1] = b_q[1];
  b_v[1] = v[1];
  b_q[1] = q[1];
  stage(b_q, b_v, force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, force, d, b_q, b_v);
  stage(b_q, b_v, force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, force, d, b_q, b_v);
  stage(b_q, b_v, force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, force, d, b_q, b_v);
  stage(b_q, b_v, force, d, q1, v1);
  b_r = _mm_loadu_pd(&q1[0]);
  r1 = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&q1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  b_r = _mm_loadu_pd(&v1[0]);
  r1 = _mm_loadu_pd(&b_v[0]);
  _mm_storeu_pd(&v1[0], _mm_add_pd(r1, _mm_mul_pd(r2, _mm_sub_pd(b_r, r1))));
  stage(q1, v1, force, d, b_q, b_v);
  b_r = _mm_loadu_pd(&b_q[0]);
  _mm_storeu_pd(&Q[2], b_r);
  _mm_storeu_pd(&Q[4],
                _mm_add_pd(_mm_loadu_pd(&q[0]),
                           _mm_mul_pd(_mm_set1_pd(dt), _mm_loadu_pd(&v[0]))));
}

/*
 * File trailer for residualShadowFlow.c
 *
 * [EOF]
 */
