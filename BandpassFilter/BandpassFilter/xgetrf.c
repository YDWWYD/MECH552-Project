/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgetrf.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xgetrf.h"
#include "xzgetrf.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                emxArray_creal_T *A
 *                int lda
 *                emxArray_int32_T *ipiv
 *                int *info
 * Return Type  : void
 */
void xgetrf(int m, int n, emxArray_creal_T *A, int lda, emxArray_int32_T *ipiv,
            int *info)
{
  int b_info;
  xzgetrf(m, n, A, lda, ipiv, &b_info);
  *info = b_info;
}

/*
 * File trailer for xgetrf.c
 *
 * [EOF]
 */
