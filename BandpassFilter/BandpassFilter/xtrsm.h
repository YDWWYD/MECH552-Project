/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xtrsm.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef XTRSM_H
#define XTRSM_H

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
extern void b_xtrsm(int m, const emxArray_creal_T *A, int lda, emxArray_creal_T *
                    B);
extern void c_xtrsm(int m, const emxArray_creal_T *A, int lda, emxArray_creal_T *
                    B);
extern void d_xtrsm(int n, const emxArray_creal_T *A, int lda, emxArray_creal_T *
                    B);
extern void e_xtrsm(int n, const emxArray_creal_T *A, int lda, emxArray_creal_T *
                    B);
extern void xtrsm(int m, int n, const emxArray_creal_T *A, int lda,
                  emxArray_creal_T *B, int ldb);

#endif

/*
 * File trailer for xtrsm.h
 *
 * [EOF]
 */
