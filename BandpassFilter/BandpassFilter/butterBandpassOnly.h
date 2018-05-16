/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: butterBandpassOnly.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef BUTTERBANDPASSONLY_H
#define BUTTERBANDPASSONLY_H

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
extern void butterBandpassOnly(double n, const double Wn[2], emxArray_real_T
  *num, emxArray_creal_T *den);

#endif

/*
 * File trailer for butterBandpassOnly.h
 *
 * [EOF]
 */
