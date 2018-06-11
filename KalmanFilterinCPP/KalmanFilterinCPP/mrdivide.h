/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mrdivide.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef MRDIVIDE_H
#define MRDIVIDE_H

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
extern void b_mrdivide(const emxArray_creal_T *A, const emxArray_creal_T *B,
  emxArray_creal_T *y);
extern void mrdivide(const emxArray_real_T *A, const emxArray_creal_T *B,
                     emxArray_creal_T *y);

#endif

/*
 * File trailer for mrdivide.h
 *
 * [EOF]
 */
