/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzhgeqz.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef XZHGEQZ_H
#define XZHGEQZ_H

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
extern void xzhgeqz(const emxArray_creal_T *A, int ilo, int ihi, int *info,
                    emxArray_creal_T *alpha1, emxArray_creal_T *beta1);

#endif

/*
 * File trailer for xzhgeqz.h
 *
 * [EOF]
 */
