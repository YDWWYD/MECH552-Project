/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgerc.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xgerc.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                const creal_T alpha1
 *                int ix0
 *                const emxArray_creal_T *y
 *                emxArray_creal_T *A
 *                int ia0
 *                int lda
 * Return Type  : void
 */
void xgerc(int m, int n, const creal_T alpha1, int ix0, const emxArray_creal_T
           *y, emxArray_creal_T *A, int ia0, int lda)
{
  int jA;
  int jy;
  int j;
  boolean_T b_y;
  double temp_re;
  double temp_im;
  int ix;
  int i15;
  int ijA;
  double A_re;
  double A_im;
  if (!((alpha1.re == 0.0) && (alpha1.im == 0.0))) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 1; j <= n; j++) {
      b_y = ((y->data[jy].re != 0.0) || (y->data[jy].im != 0.0));
      if (b_y) {
        temp_re = y->data[jy].re * alpha1.re + y->data[jy].im * alpha1.im;
        temp_im = y->data[jy].re * alpha1.im - y->data[jy].im * alpha1.re;
        ix = ix0;
        i15 = m + jA;
        for (ijA = jA; ijA + 1 <= i15; ijA++) {
          A_re = A->data[ix - 1].re * temp_re - A->data[ix - 1].im * temp_im;
          A_im = A->data[ix - 1].re * temp_im + A->data[ix - 1].im * temp_re;
          A->data[ijA].re += A_re;
          A->data[ijA].im += A_im;
          ix++;
        }
      }

      jy++;
      jA += lda;
    }
  }
}

/*
 * File trailer for xgerc.c
 *
 * [EOF]
 */
