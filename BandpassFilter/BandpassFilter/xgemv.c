/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgemv.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xgemv.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                const emxArray_creal_T *A
 *                int ia0
 *                int lda
 *                const emxArray_creal_T *x
 *                int ix0
 *                emxArray_creal_T *y
 * Return Type  : void
 */
void xgemv(int m, int n, const emxArray_creal_T *A, int ia0, int lda, const
           emxArray_creal_T *x, int ix0, emxArray_creal_T *y)
{
  int iy;
  int i13;
  int iac;
  int ix;
  double c_re;
  double c_im;
  int i14;
  int ia;
  if (n != 0) {
    for (iy = 1; iy <= n; iy++) {
      y->data[iy - 1].re = 0.0;
      y->data[iy - 1].im = 0.0;
    }

    iy = 0;
    i13 = ia0 + lda * (n - 1);
    iac = ia0;
    while ((lda > 0) && (iac <= i13)) {
      ix = ix0 - 1;
      c_re = 0.0;
      c_im = 0.0;
      i14 = (iac + m) - 1;
      for (ia = iac - 1; ia + 1 <= i14; ia++) {
        c_re += A->data[ia].re * x->data[ix].re + A->data[ia].im * x->data[ix].
          im;
        c_im += A->data[ia].re * x->data[ix].im - A->data[ia].im * x->data[ix].
          re;
        ix++;
      }

      y->data[iy].re += c_re - 0.0 * c_im;
      y->data[iy].im += c_im + 0.0 * c_re;
      iy++;
      iac += lda;
    }
  }
}

/*
 * File trailer for xgemv.c
 *
 * [EOF]
 */
