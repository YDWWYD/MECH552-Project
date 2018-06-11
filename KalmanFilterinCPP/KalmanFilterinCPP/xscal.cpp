/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xscal.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xscal.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T a
 *                emxArray_creal_T *x
 *                int ix0
 * Return Type  : void
 */
void b_xscal(const creal_T a, emxArray_creal_T *x, int ix0)
{
  int k;
  double x_re;
  double x_im;
  for (k = ix0; k < ix0; k++) {
    x_re = x->data[k - 1].re;
    x_im = x->data[k - 1].im;
    x->data[k - 1].re = a.re * x_re - a.im * x_im;
    x->data[k - 1].im = a.re * x_im + a.im * x_re;
  }
}

/*
 * Arguments    : int n
 *                const creal_T a
 *                emxArray_creal_T *x
 *                int ix0
 *                int incx
 * Return Type  : void
 */
void c_xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0, int incx)
{
  int i17;
  int k;
  double x_re;
  double x_im;
  if (!(incx < 1)) {
    i17 = ix0 + incx * (n - 1);
    for (k = ix0; k <= i17; k += incx) {
      x_re = x->data[k - 1].re;
      x_im = x->data[k - 1].im;
      x->data[k - 1].re = a.re * x_re - a.im * x_im;
      x->data[k - 1].im = a.re * x_im + a.im * x_re;
    }
  }
}

/*
 * Arguments    : int n
 *                const creal_T a
 *                emxArray_creal_T *x
 *                int ix0
 * Return Type  : void
 */
void xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0)
{
  int i12;
  int k;
  double x_re;
  double x_im;
  i12 = (ix0 + n) - 1;
  for (k = ix0; k <= i12; k++) {
    x_re = x->data[k - 1].re;
    x_im = x->data[k - 1].im;
    x->data[k - 1].re = a.re * x_re - a.im * x_im;
    x->data[k - 1].im = a.re * x_im + a.im * x_re;
  }
}

/*
 * File trailer for xscal.c
 *
 * [EOF]
 */
