/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qrsolve.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "qrsolve.h"
#include "butterBandpassOnly_emxutil.h"
#include "xunormqr.h"
#include "xgeqp3.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *A
 *                const emxArray_creal_T *B
 *                emxArray_creal_T *Y
 * Return Type  : void
 */
void b_qrsolve(const emxArray_creal_T *A, const emxArray_creal_T *B,
               emxArray_creal_T *Y)
{
  emxArray_creal_T *b_A;
  int i;
  int loop_ub;
  emxArray_creal_T *tau;
  emxArray_int32_T *jpvt;
  int rankA;
  emxArray_creal_T *b_B;
  double Y_re;
  double Y_im;
  double A_re;
  double A_im;
  double brm;
  double bim;
  double s;
  emxInit_creal_T(&b_A, 2);
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity_creal_T(b_A, i);
  loop_ub = A->size[0] * A->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_A->data[i] = A->data[i];
  }

  emxInit_creal_T1(&tau, 1);
  emxInit_int32_T1(&jpvt, 2);
  xgeqp3(b_A, tau, jpvt);
  rankA = rankFromQR(b_A);
  loop_ub = b_A->size[1];
  i = Y->size[0];
  Y->size[0] = loop_ub;
  emxEnsureCapacity_creal_T1(Y, i);
  for (i = 0; i < loop_ub; i++) {
    Y->data[i].re = 0.0;
    Y->data[i].im = 0.0;
  }

  emxInit_creal_T1(&b_B, 1);
  i = b_B->size[0];
  b_B->size[0] = B->size[0];
  emxEnsureCapacity_creal_T1(b_B, i);
  loop_ub = B->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_B->data[i] = B->data[i];
  }

  xunormqr(b_A, b_B, tau);
  i = 0;
  emxFree_creal_T(&tau);
  while (i + 1 <= rankA) {
    Y->data[jpvt->data[i] - 1] = b_B->data[i];
    i++;
  }

  emxFree_creal_T(&b_B);
  for (loop_ub = rankA - 1; loop_ub + 1 > 0; loop_ub--) {
    Y_re = Y->data[jpvt->data[loop_ub] - 1].re;
    Y_im = Y->data[jpvt->data[loop_ub] - 1].im;
    A_re = b_A->data[loop_ub + b_A->size[0] * loop_ub].re;
    A_im = b_A->data[loop_ub + b_A->size[0] * loop_ub].im;
    if (A_im == 0.0) {
      if (Y_im == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_re / A_re;
        Y->data[jpvt->data[loop_ub] - 1].im = 0.0;
      } else if (Y_re == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = 0.0;
        Y->data[jpvt->data[loop_ub] - 1].im = Y_im / A_re;
      } else {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_re / A_re;
        Y->data[jpvt->data[loop_ub] - 1].im = Y_im / A_re;
      }
    } else if (A_re == 0.0) {
      if (Y_re == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_im / A_im;
        Y->data[jpvt->data[loop_ub] - 1].im = 0.0;
      } else if (Y_im == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = 0.0;
        Y->data[jpvt->data[loop_ub] - 1].im = -(Y_re / A_im);
      } else {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_im / A_im;
        Y->data[jpvt->data[loop_ub] - 1].im = -(Y_re / A_im);
      }
    } else {
      brm = fabs(A_re);
      bim = fabs(A_im);
      if (brm > bim) {
        s = A_im / A_re;
        bim = A_re + s * A_im;
        Y->data[jpvt->data[loop_ub] - 1].re = (Y_re + s * Y_im) / bim;
        Y->data[jpvt->data[loop_ub] - 1].im = (Y_im - s * Y_re) / bim;
      } else if (bim == brm) {
        if (A_re > 0.0) {
          s = 0.5;
        } else {
          s = -0.5;
        }

        if (A_im > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        Y->data[jpvt->data[loop_ub] - 1].re = (Y_re * s + Y_im * bim) / brm;
        Y->data[jpvt->data[loop_ub] - 1].im = (Y_im * s - Y_re * bim) / brm;
      } else {
        s = A_re / A_im;
        bim = A_im + s * A_re;
        Y->data[jpvt->data[loop_ub] - 1].re = (s * Y_re + Y_im) / bim;
        Y->data[jpvt->data[loop_ub] - 1].im = (s * Y_im - Y_re) / bim;
      }
    }

    for (i = 0; i + 1 <= loop_ub; i++) {
      Y_re = Y->data[jpvt->data[loop_ub] - 1].re * b_A->data[i + b_A->size[0] *
        loop_ub].re - Y->data[jpvt->data[loop_ub] - 1].im * b_A->data[i +
        b_A->size[0] * loop_ub].im;
      Y_im = Y->data[jpvt->data[loop_ub] - 1].re * b_A->data[i + b_A->size[0] *
        loop_ub].im + Y->data[jpvt->data[loop_ub] - 1].im * b_A->data[i +
        b_A->size[0] * loop_ub].re;
      Y->data[jpvt->data[i] - 1].re -= Y_re;
      Y->data[jpvt->data[i] - 1].im -= Y_im;
    }
  }

  emxFree_int32_T(&jpvt);
  emxFree_creal_T(&b_A);
}

/*
 * Arguments    : const emxArray_creal_T *A
 *                const emxArray_real_T *B
 *                emxArray_creal_T *Y
 * Return Type  : void
 */
void qrsolve(const emxArray_creal_T *A, const emxArray_real_T *B,
             emxArray_creal_T *Y)
{
  emxArray_creal_T *b_A;
  int i;
  int loop_ub;
  emxArray_creal_T *tau;
  emxArray_int32_T *jpvt;
  int rankA;
  emxArray_creal_T *CB;
  double Y_re;
  double Y_im;
  double A_re;
  double A_im;
  double brm;
  double bim;
  double s;
  emxInit_creal_T(&b_A, 2);
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity_creal_T(b_A, i);
  loop_ub = A->size[0] * A->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_A->data[i] = A->data[i];
  }

  emxInit_creal_T1(&tau, 1);
  emxInit_int32_T1(&jpvt, 2);
  xgeqp3(b_A, tau, jpvt);
  rankA = rankFromQR(b_A);
  loop_ub = b_A->size[1];
  i = Y->size[0];
  Y->size[0] = loop_ub;
  emxEnsureCapacity_creal_T1(Y, i);
  for (i = 0; i < loop_ub; i++) {
    Y->data[i].re = 0.0;
    Y->data[i].im = 0.0;
  }

  emxInit_creal_T1(&CB, 1);
  i = CB->size[0];
  CB->size[0] = B->size[0];
  emxEnsureCapacity_creal_T1(CB, i);
  loop_ub = B->size[0];
  for (i = 0; i < loop_ub; i++) {
    CB->data[i].re = B->data[i];
    CB->data[i].im = 0.0;
  }

  xunormqr(b_A, CB, tau);
  i = 0;
  emxFree_creal_T(&tau);
  while (i + 1 <= rankA) {
    Y->data[jpvt->data[i] - 1] = CB->data[i];
    i++;
  }

  emxFree_creal_T(&CB);
  for (loop_ub = rankA - 1; loop_ub + 1 > 0; loop_ub--) {
    Y_re = Y->data[jpvt->data[loop_ub] - 1].re;
    Y_im = Y->data[jpvt->data[loop_ub] - 1].im;
    A_re = b_A->data[loop_ub + b_A->size[0] * loop_ub].re;
    A_im = b_A->data[loop_ub + b_A->size[0] * loop_ub].im;
    if (A_im == 0.0) {
      if (Y_im == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_re / A_re;
        Y->data[jpvt->data[loop_ub] - 1].im = 0.0;
      } else if (Y_re == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = 0.0;
        Y->data[jpvt->data[loop_ub] - 1].im = Y_im / A_re;
      } else {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_re / A_re;
        Y->data[jpvt->data[loop_ub] - 1].im = Y_im / A_re;
      }
    } else if (A_re == 0.0) {
      if (Y_re == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_im / A_im;
        Y->data[jpvt->data[loop_ub] - 1].im = 0.0;
      } else if (Y_im == 0.0) {
        Y->data[jpvt->data[loop_ub] - 1].re = 0.0;
        Y->data[jpvt->data[loop_ub] - 1].im = -(Y_re / A_im);
      } else {
        Y->data[jpvt->data[loop_ub] - 1].re = Y_im / A_im;
        Y->data[jpvt->data[loop_ub] - 1].im = -(Y_re / A_im);
      }
    } else {
      brm = fabs(A_re);
      bim = fabs(A_im);
      if (brm > bim) {
        s = A_im / A_re;
        bim = A_re + s * A_im;
        Y->data[jpvt->data[loop_ub] - 1].re = (Y_re + s * Y_im) / bim;
        Y->data[jpvt->data[loop_ub] - 1].im = (Y_im - s * Y_re) / bim;
      } else if (bim == brm) {
        if (A_re > 0.0) {
          s = 0.5;
        } else {
          s = -0.5;
        }

        if (A_im > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        Y->data[jpvt->data[loop_ub] - 1].re = (Y_re * s + Y_im * bim) / brm;
        Y->data[jpvt->data[loop_ub] - 1].im = (Y_im * s - Y_re * bim) / brm;
      } else {
        s = A_re / A_im;
        bim = A_im + s * A_re;
        Y->data[jpvt->data[loop_ub] - 1].re = (s * Y_re + Y_im) / bim;
        Y->data[jpvt->data[loop_ub] - 1].im = (s * Y_im - Y_re) / bim;
      }
    }

    for (i = 0; i + 1 <= loop_ub; i++) {
      Y_re = Y->data[jpvt->data[loop_ub] - 1].re * b_A->data[i + b_A->size[0] *
        loop_ub].re - Y->data[jpvt->data[loop_ub] - 1].im * b_A->data[i +
        b_A->size[0] * loop_ub].im;
      Y_im = Y->data[jpvt->data[loop_ub] - 1].re * b_A->data[i + b_A->size[0] *
        loop_ub].im + Y->data[jpvt->data[loop_ub] - 1].im * b_A->data[i +
        b_A->size[0] * loop_ub].re;
      Y->data[jpvt->data[i] - 1].re -= Y_re;
      Y->data[jpvt->data[i] - 1].im -= Y_im;
    }
  }

  emxFree_int32_T(&jpvt);
  emxFree_creal_T(&b_A);
}

/*
 * Arguments    : const emxArray_creal_T *A
 * Return Type  : int
 */
int rankFromQR(const emxArray_creal_T *A)
{
  int r;
  int minmn;
  int maxmn;
  double tol;
  r = 0;
  if (A->size[0] < A->size[1]) {
    minmn = A->size[0];
    maxmn = A->size[1];
  } else {
    minmn = A->size[1];
    maxmn = A->size[0];
  }

  if (minmn > 0) {
    tol = (double)maxmn * (fabs(A->data[0].re) + fabs(A->data[0].im)) *
      2.2204460492503131E-16;
    while ((r < minmn) && (!(fabs(A->data[r + A->size[0] * r].re) + fabs(A->
              data[r + A->size[0] * r].im) <= tol))) {
      r++;
    }
  }

  return r;
}

/*
 * File trailer for qrsolve.c
 *
 * [EOF]
 */
