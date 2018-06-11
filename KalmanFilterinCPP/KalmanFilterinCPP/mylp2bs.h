/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mylp2bs.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef MYLP2BS_H
#define MYLP2BS_H

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
extern void mylp2bs(const emxArray_creal_T *a, const emxArray_real_T *b, const
                    emxArray_real_T *c, double d, double wo, const double
                    bw_data[], emxArray_creal_T *at, emxArray_creal_T *bt,
                    emxArray_creal_T *ct, creal_T *dt);

#endif

/*
 * File trailer for mylp2bs.h
 *
 * [EOF]
 */
