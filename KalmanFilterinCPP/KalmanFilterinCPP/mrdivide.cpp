/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mrdivide.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "mrdivide.h"
#include "butterBandpassOnly_emxutil.h"
#include "xtrsm.h"
#include "xgetrf.h"
#include "qrsolve.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *A
 *                const emxArray_creal_T *B
 *                emxArray_creal_T *y
 * Return Type  : void
 */
void b_mrdivide(const emxArray_creal_T *A, const emxArray_creal_T *B,
                emxArray_creal_T *y)
{
  emxArray_int32_T *ipiv;
  emxArray_creal_T *r5;
  emxArray_creal_T *b_B;
  emxArray_creal_T *b_A;
  unsigned int unnamed_idx_1;
  int info;
  int loop_ub;
  int b_loop_ub;
  int i6;
  double temp_re;
  double temp_im;
  emxInit_int32_T1(&ipiv, 2);
  emxInit_creal_T1(&r5, 1);
  emxInit_creal_T(&b_B, 2);
  emxInit_creal_T1(&b_A, 1);
  if ((A->size[1] == 0) || ((B->size[0] == 0) || (B->size[1] == 0))) {
    unnamed_idx_1 = (unsigned int)B->size[0];
    info = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity_creal_T(y, info);
    loop_ub = (int)unnamed_idx_1;
    for (info = 0; info < loop_ub; info++) {
      y->data[info].re = 0.0;
      y->data[info].im = 0.0;
    }
  } else if (B->size[0] == B->size[1]) {
    info = b_B->size[0] * b_B->size[1];
    b_B->size[0] = B->size[0];
    b_B->size[1] = B->size[1];
    emxEnsureCapacity_creal_T(b_B, info);
    loop_ub = B->size[0] * B->size[1];
    for (info = 0; info < loop_ub; info++) {
      b_B->data[info] = B->data[info];
    }

    xgetrf(B->size[1], B->size[1], b_B, B->size[1], ipiv, &info);
    info = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(y, info);
    loop_ub = A->size[0] * A->size[1];
    for (info = 0; info < loop_ub; info++) {
      y->data[info] = A->data[info];
    }

    d_xtrsm(B->size[1], b_B, B->size[1], y);
    e_xtrsm(B->size[1], b_B, B->size[1], y);
    for (info = B->size[1] - 2; info + 1 > 0; info--) {
      if (ipiv->data[info] != info + 1) {
        temp_re = y->data[y->size[0] * info].re;
        temp_im = y->data[y->size[0] * info].im;
        y->data[y->size[0] * info] = y->data[y->size[0] * (ipiv->data[info] - 1)];
        y->data[y->size[0] * (ipiv->data[info] - 1)].re = temp_re;
        y->data[y->size[0] * (ipiv->data[info] - 1)].im = temp_im;
      }
    }
  } else {
    info = b_B->size[0] * b_B->size[1];
    b_B->size[0] = B->size[1];
    b_B->size[1] = B->size[0];
    emxEnsureCapacity_creal_T(b_B, info);
    loop_ub = B->size[0];
    for (info = 0; info < loop_ub; info++) {
      b_loop_ub = B->size[1];
      for (i6 = 0; i6 < b_loop_ub; i6++) {
        b_B->data[i6 + b_B->size[0] * info].re = B->data[info + B->size[0] * i6]
          .re;
        b_B->data[i6 + b_B->size[0] * info].im = -B->data[info + B->size[0] * i6]
          .im;
      }
    }

    info = b_A->size[0];
    b_A->size[0] = A->size[1];
    emxEnsureCapacity_creal_T1(b_A, info);
    loop_ub = A->size[1];
    for (info = 0; info < loop_ub; info++) {
      b_A->data[info].re = A->data[A->size[0] * info].re;
      b_A->data[info].im = -A->data[A->size[0] * info].im;
    }

    b_qrsolve(b_B, b_A, r5);
    info = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = r5->size[0];
    emxEnsureCapacity_creal_T(y, info);
    loop_ub = r5->size[0];
    for (info = 0; info < loop_ub; info++) {
      y->data[y->size[0] * info].re = r5->data[info].re;
      y->data[y->size[0] * info].im = -r5->data[info].im;
    }
  }

  emxFree_creal_T(&b_A);
  emxFree_creal_T(&b_B);
  emxFree_creal_T(&r5);
  emxFree_int32_T(&ipiv);
}

/*
 * Arguments    : const emxArray_real_T *A
 *                const emxArray_creal_T *B
 *                emxArray_creal_T *y
 * Return Type  : void
 */
void mrdivide(const emxArray_real_T *A, const emxArray_creal_T *B,
              emxArray_creal_T *y)
{
  emxArray_int32_T *ipiv;
  emxArray_creal_T *r3;
  emxArray_creal_T *b_B;
  emxArray_real_T *b_A;
  unsigned int unnamed_idx_1;
  int info;
  int loop_ub;
  int b_loop_ub;
  int i5;
  double temp_re;
  double temp_im;
  emxInit_int32_T1(&ipiv, 2);
  emxInit_creal_T1(&r3, 1);
  emxInit_creal_T(&b_B, 2);
  emxInit_real_T(&b_A, 1);
  if ((A->size[1] == 0) || ((B->size[0] == 0) || (B->size[1] == 0))) {
    unnamed_idx_1 = (unsigned int)B->size[0];
    info = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity_creal_T(y, info);
    loop_ub = (int)unnamed_idx_1;
    for (info = 0; info < loop_ub; info++) {
      y->data[info].re = 0.0;
      y->data[info].im = 0.0;
    }
  } else if (B->size[0] == B->size[1]) {
    info = b_B->size[0] * b_B->size[1];
    b_B->size[0] = B->size[0];
    b_B->size[1] = B->size[1];
    emxEnsureCapacity_creal_T(b_B, info);
    loop_ub = B->size[0] * B->size[1];
    for (info = 0; info < loop_ub; info++) {
      b_B->data[info] = B->data[info];
    }

    xgetrf(B->size[1], B->size[1], b_B, B->size[1], ipiv, &info);
    info = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(y, info);
    loop_ub = A->size[0] * A->size[1];
    for (info = 0; info < loop_ub; info++) {
      y->data[info].re = A->data[info];
      y->data[info].im = 0.0;
    }

    d_xtrsm(B->size[1], b_B, B->size[1], y);
    e_xtrsm(B->size[1], b_B, B->size[1], y);
    for (info = B->size[1] - 2; info + 1 > 0; info--) {
      if (ipiv->data[info] != info + 1) {
        temp_re = y->data[y->size[0] * info].re;
        temp_im = y->data[y->size[0] * info].im;
        y->data[y->size[0] * info] = y->data[y->size[0] * (ipiv->data[info] - 1)];
        y->data[y->size[0] * (ipiv->data[info] - 1)].re = temp_re;
        y->data[y->size[0] * (ipiv->data[info] - 1)].im = temp_im;
      }
    }
  } else {
    info = b_B->size[0] * b_B->size[1];
    b_B->size[0] = B->size[1];
    b_B->size[1] = B->size[0];
    emxEnsureCapacity_creal_T(b_B, info);
    loop_ub = B->size[0];
    for (info = 0; info < loop_ub; info++) {
      b_loop_ub = B->size[1];
      for (i5 = 0; i5 < b_loop_ub; i5++) {
        b_B->data[i5 + b_B->size[0] * info].re = B->data[info + B->size[0] * i5]
          .re;
        b_B->data[i5 + b_B->size[0] * info].im = -B->data[info + B->size[0] * i5]
          .im;
      }
    }

    info = b_A->size[0];
    b_A->size[0] = A->size[1];
    emxEnsureCapacity_real_T1(b_A, info);
    loop_ub = A->size[1];
    for (info = 0; info < loop_ub; info++) {
      b_A->data[info] = A->data[A->size[0] * info];
    }

    qrsolve(b_B, b_A, r3);
    info = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = r3->size[0];
    emxEnsureCapacity_creal_T(y, info);
    loop_ub = r3->size[0];
    for (info = 0; info < loop_ub; info++) {
      y->data[y->size[0] * info].re = r3->data[info].re;
      y->data[y->size[0] * info].im = -r3->data[info].im;
    }
  }

  emxFree_real_T(&b_A);
  emxFree_creal_T(&b_B);
  emxFree_creal_T(&r3);
  emxFree_int32_T(&ipiv);
}

/*
 * File trailer for mrdivide.c
 *
 * [EOF]
 */
