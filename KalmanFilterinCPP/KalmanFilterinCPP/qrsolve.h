/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qrsolve.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef QRSOLVE_H
#define QRSOLVE_H

/* Include Files */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "butterBandpassOnly_types.h"

/* Function Declarations */
extern void b_qrsolve(const emxArray_creal_T *A, const emxArray_creal_T *B,
                      emxArray_creal_T *Y);
extern void qrsolve(const emxArray_creal_T *A, const emxArray_real_T *B,
                    emxArray_creal_T *Y);
extern int rankFromQR(const emxArray_creal_T *A);

#endif

/*
 * File trailer for qrsolve.h
 *
 * [EOF]
 */
