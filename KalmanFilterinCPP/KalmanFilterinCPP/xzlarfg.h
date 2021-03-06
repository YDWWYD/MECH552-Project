/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzlarfg.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef XZLARFG_H
#define XZLARFG_H

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
extern creal_T b_xzlarfg(creal_T *alpha1, creal_T *x);
extern creal_T xzlarfg(int n, creal_T *alpha1, emxArray_creal_T *x, int ix0);

#endif

/*
 * File trailer for xzlarfg.h
 *
 * [EOF]
 */
