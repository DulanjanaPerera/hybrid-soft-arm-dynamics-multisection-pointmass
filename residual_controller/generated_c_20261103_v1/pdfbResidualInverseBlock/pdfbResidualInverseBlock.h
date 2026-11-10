/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pdfbResidualInverseBlock.h
 *
 * MATLAB Coder version            : 25.1
 * C/C++ source code generated on  : 03-Nov-2026 18:16:46
 */

#ifndef PDFBRESIDUALINVERSEBLOCK_H
#define PDFBRESIDUALINVERSEBLOCK_H

/* Include Files */
#include "pdfbResidualInverseBlock_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void pdfbResidualInverseBlock(
    const double measured[6], const double desiredAcceleration[2],
    boolean_T orientationValid, double t, const double gravity[3],
    boolean_T baseOrientationValid, double measurementTime_s,
    boolean_T enablePressure, boolean_T enableResidual,
    boolean_T enableCommonPressureAdaptation, const double pressureOptions[10],
    const double supervisorOptions[7], const double residualOptions[8],
    double pressure_bar[3], double tauDesired_N[2], double tauAvailable_N[2],
    boolean_T *limited, boolean_T *modelValid, double MassForce_N[2],
    double CoriolisForce_N[2], double GnStiff_N[2], unsigned char *status,
    struct0_T *residualDiagnostic);

void pdfbResidualInverseBlock_init(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for pdfbResidualInverseBlock.h
 *
 * [EOF]
 */
