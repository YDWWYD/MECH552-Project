/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: myzp2ss.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef MYZP2SS_H
#define MYZP2SS_H

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
extern void myzp2ss(emxArray_creal_T *p, double k, emxArray_creal_T *a,
                    emxArray_real_T *b, emxArray_real_T *c, double *d);

#endif

/*
 * File trailer for myzp2ss.h
 *
 * [EOF]
 */
