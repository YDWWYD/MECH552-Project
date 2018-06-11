/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgerc.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef XGERC_H
#define XGERC_H

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
extern void xgerc(int m, int n, const creal_T alpha1, int ix0, const
                  emxArray_creal_T *y, emxArray_creal_T *A, int ia0, int lda);

#endif

/*
 * File trailer for xgerc.h
 *
 * [EOF]
 */
