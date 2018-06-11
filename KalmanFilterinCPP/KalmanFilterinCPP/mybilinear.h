/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mybilinear.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef MYBILINEAR_H
#define MYBILINEAR_H

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
extern void mybilinear(const emxArray_creal_T *z, const emxArray_creal_T *p,
  const emxArray_creal_T *k, const creal_T fs, emxArray_creal_T *zd,
  emxArray_creal_T *pd, emxArray_creal_T *kd, creal_T *dd);

#endif

/*
 * File trailer for mybilinear.h
 *
 * [EOF]
 */
