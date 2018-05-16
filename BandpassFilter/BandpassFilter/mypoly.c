/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mypoly.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "mypoly.h"
#include "butterBandpassOnly_emxutil.h"
#include "isequal.h"
#include "sortButter.h"
#include "schur.h"
#include "eig.h"
#include "xzgeev.h"
#include "anyNonFinite.h"

/* Function Definitions */

/*
 * POLY Convert roots to polynomial.
 *    POLY(A), when A is an N by N matrix, is a row vector with
 *    N+1 elements which are the coefficients of the
 *    characteristic polynomial, det(lambda*eye(size(A)) - A).
 *
 *    POLY(V), when V is a vector, is a vector whose elements are
 *    the coefficients of the polynomial whose roots are the
 *    elements of V. For vectors, ROOTS and POLY are inverse
 *    functions of each other, up to ordering, scaling, and
 *    roundoff error.
 *
 *    Examples:
 *
 *    roots(poly(1:20)) generates Wilkinson's famous example.
 *
 *    Class support for inputs A,V:
 *       float: double, single
 *
 *    See also ROOTS, CONV, RESIDUE, POLYVAL.
 * Arguments    : const emxArray_creal_T *x
 *                emxArray_creal_T *c
 * Return Type  : void
 */
void mypoly(const emxArray_creal_T *x, emxArray_creal_T *c)
{
  emxArray_creal_T *e;
  emxArray_creal_T *alpha1;
  emxArray_creal_T *beta1;
  int i7;
  emxArray_creal_T *b_e;
  int info;
  boolean_T p;
  int j;
  boolean_T exitg2;
  emxArray_creal_T *b_x;
  int i;
  int exitg1;
  double x_re;
  double x_im;
  boolean_T c_x;
  emxArray_creal_T *b_c;
  emxArray_real_T *r6;
  double e_re;
  double e_im;
  emxArray_int32_T *r7;
  emxArray_int32_T *r8;

  /*    Copyright 1984-2014 The MathWorks, Inc. */
  /*  if m==n square matrix */
  emxInit_creal_T(&e, 2);
  emxInit_creal_T1(&alpha1, 1);
  emxInit_creal_T1(&beta1, 1);
  if (x->size[0] == x->size[1]) {
    /*  Characteristic polynomial (square x) */
    emxInit_creal_T1(&b_e, 1);
    if ((x->size[0] == 0) || (x->size[1] == 0)) {
      i7 = b_e->size[0];
      b_e->size[0] = x->size[0];
      emxEnsureCapacity_creal_T1(b_e, i7);
      info = x->size[0];
      for (i7 = 0; i7 < info; i7++) {
        b_e->data[i7].re = 0.0;
        b_e->data[i7].im = 0.0;
      }
    } /*else if (anyNonFinite(x)) {
      if ((x->size[0] == 1) && (x->size[1] == 1)) {
        i7 = b_e->size[0];
        b_e->size[0] = 1;
        emxEnsureCapacity_creal_T1(b_e, i7);
        b_e->data[0].re = rtNaN;
        b_e->data[0].im = 0.0;
      } else {
        i7 = b_e->size[0];
        b_e->size[0] = x->size[0];
        emxEnsureCapacity_creal_T1(b_e, i7);
        info = x->size[0];
        for (i7 = 0; i7 < info; i7++) {
          b_e->data[i7].re = rtNaN;
          b_e->data[i7].im = 0.0;
        }
      }
    }*/ else if ((x->size[0] == 1) && (x->size[1] == 1)) {
      i7 = b_e->size[0];
      b_e->size[0] = 1;
      emxEnsureCapacity_creal_T1(b_e, i7);
      b_e->data[0] = x->data[0];
    } else {
      p = (x->size[0] == x->size[1]);
      if (p) {
        j = 0;
        exitg2 = false;
        while ((!exitg2) && (j <= x->size[1] - 1)) {
          i = 0;
          do {
            exitg1 = 0;
            if (i <= j) {
              x_re = x->data[j + x->size[0] * i].re;
              x_im = -x->data[j + x->size[0] * i].im;
              c_x = ((x->data[i + x->size[0] * j].re == x_re) && (x->data[i +
                      x->size[0] * j].im == x_im));
              if (!c_x) {
                p = false;
                exitg1 = 1;
              } else {
                i++;
              }
            } else {
              j++;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      }

      if (p) {
        emxInit_creal_T(&b_x, 2);
        i7 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = x->size[0];
        b_x->size[1] = x->size[1];
        emxEnsureCapacity_creal_T(b_x, i7);
        info = x->size[0] * x->size[1];
        for (i7 = 0; i7 < info; i7++) {
          b_x->data[i7] = x->data[i7];
        }

        schur(b_x, e);
        i7 = b_e->size[0];
        b_e->size[0] = e->size[0];
        emxEnsureCapacity_creal_T1(b_e, i7);
        info = 0;
        emxFree_creal_T(&b_x);
        while (info + 1 <= e->size[0]) {
          b_e->data[info] = e->data[info + e->size[0] * info];
          info++;
        }
      } else {
        xzgeev(x, &info, alpha1, beta1);
        makeD(alpha1, beta1, b_e);
      }
    }

    i7 = e->size[0] * e->size[1];
    e->size[0] = b_e->size[0];
    e->size[1] = 1;
    emxEnsureCapacity_creal_T(e, i7);
    info = b_e->size[0];
    for (i7 = 0; i7 < info; i7++) {
      e->data[i7] = b_e->data[i7];
    }

    emxFree_creal_T(&b_e);

    /* elseif (m==1) || (n==1) */
  } else if ((x->size[0] == 1) || (x->size[1] == 1)) {
    i7 = e->size[0] * e->size[1];
    e->size[0] = x->size[0];
    e->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(e, i7);
    info = x->size[0] * x->size[1];
    for (i7 = 0; i7 < info; i7++) {
      e->data[i7] = x->data[i7];
    }
  } else {
    /* error(message('MATLAB:poly:InputSize')) */
    i7 = e->size[0] * e->size[1];
    e->size[0] = 0;
    e->size[1] = 0;
    emxEnsureCapacity_creal_T(e, i7);
  }

  /*  Strip out infinities */
  /*  e = e( isfinite(e) ); */
  /*  Expand recursion formula */
  if ((e->size[0] == 0) || (e->size[1] == 0)) {
    i = 0;
  } else {
    info = e->size[0];
    i = e->size[1];
    if (info > i) {
      i = info;
    }
  }

  i7 = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = 1 + i;
  emxEnsureCapacity_creal_T(c, i7);
  c->data[0].re = 1.0;
  c->data[0].im = 0.0;
  for (i7 = 0; i7 < i; i7++) {
    c->data[c->size[0] * (i7 + 1)].re = 0.0;
    c->data[c->size[0] * (i7 + 1)].im = 0.0;
  }

  /* c = [1 zeros(1,n,class(x))]; */
  /* c = [1 zeros(1,n)]; */
  j = 0;
  emxInit_creal_T(&b_c, 2);
  while (j <= i - 1) {
    info = (int)((1.0 + (double)j) + 1.0) - 2;
    x_re = e->data[j].re;
    x_im = e->data[j].im;
    i7 = b_c->size[0] * b_c->size[1];
    b_c->size[0] = 1;
    b_c->size[1] = (int)((1.0 + (double)j) + 1.0) - 1;
    emxEnsureCapacity_creal_T(b_c, i7);
    for (i7 = 0; i7 <= info; i7++) {
      e_re = x_re * c->data[i7].re - x_im * c->data[i7].im;
      e_im = x_re * c->data[i7].im + x_im * c->data[i7].re;
      b_c->data[b_c->size[0] * i7].re = c->data[1 + i7].re - e_re;
      b_c->data[b_c->size[0] * i7].im = c->data[1 + i7].im - e_im;
    }

    info = b_c->size[1];
    for (i7 = 0; i7 < info; i7++) {
      c->data[1 + i7] = b_c->data[b_c->size[0] * i7];
    }

    j++;
  }

  emxFree_creal_T(&b_c);
  emxInit_real_T1(&r6, 2);

  /*  The result should be real if the roots are complex conjugates. */
  i7 = r6->size[0] * r6->size[1];
  r6->size[0] = e->size[0];
  r6->size[1] = e->size[1];
  emxEnsureCapacity_real_T(r6, i7);
  info = e->size[0] * e->size[1];
  for (i7 = 0; i7 < info; i7++) {
    r6->data[i7] = e->data[i7].im;
  }

  j = r6->size[0] * r6->size[1] - 1;
  info = 0;
  for (i = 0; i <= j; i++) {
    if (e->data[i].im > 0.0) {
      info++;
    }
  }

  emxInit_int32_T(&r7, 1);
  i7 = r7->size[0];
  r7->size[0] = info;
  emxEnsureCapacity_int32_T(r7, i7);
  info = 0;
  for (i = 0; i <= j; i++) {
    if (e->data[i].im > 0.0) {
      r7->data[info] = i + 1;
      info++;
    }
  }

  i7 = r6->size[0] * r6->size[1];
  r6->size[0] = e->size[0];
  r6->size[1] = e->size[1];
  emxEnsureCapacity_real_T(r6, i7);
  info = e->size[0] * e->size[1];
  for (i7 = 0; i7 < info; i7++) {
    r6->data[i7] = e->data[i7].im;
  }

  j = r6->size[0] * r6->size[1] - 1;
  info = 0;
  emxFree_real_T(&r6);
  for (i = 0; i <= j; i++) {
    if (e->data[i].im < 0.0) {
      info++;
    }
  }

  emxInit_int32_T(&r8, 1);
  i7 = r8->size[0];
  r8->size[0] = info;
  emxEnsureCapacity_int32_T(r8, i7);
  info = 0;
  for (i = 0; i <= j; i++) {
    if (e->data[i].im < 0.0) {
      r8->data[info] = i + 1;
      info++;
    }
  }

  i7 = alpha1->size[0];
  alpha1->size[0] = r7->size[0];
  emxEnsureCapacity_creal_T1(alpha1, i7);
  info = r7->size[0];
  for (i7 = 0; i7 < info; i7++) {
    alpha1->data[i7] = e->data[r7->data[i7] - 1];
  }

  emxFree_int32_T(&r7);
  sortButter(alpha1);
  i7 = beta1->size[0];
  beta1->size[0] = r8->size[0];
  emxEnsureCapacity_creal_T1(beta1, i7);
  info = r8->size[0];
  for (i7 = 0; i7 < info; i7++) {
    beta1->data[i7].re = e->data[r8->data[i7] - 1].re;
    beta1->data[i7].im = -e->data[r8->data[i7] - 1].im;
  }

  emxFree_int32_T(&r8);
  emxFree_creal_T(&e);
  sortButter(beta1);
  if (isequal(alpha1, beta1)) {
    i7 = c->size[0] * c->size[1];
    c->size[0] = 1;
    emxEnsureCapacity_creal_T(c, i7);
    info = c->size[0];
    j = c->size[1];
    info *= j;
    for (i7 = 0; i7 < info; i7++) {
      x_re = c->data[i7].re;
      c->data[i7].re = x_re;
      c->data[i7].im = 0.0;
    }
  }

  emxFree_creal_T(&beta1);
  emxFree_creal_T(&alpha1);
}

/*
 * File trailer for mypoly.c
 *
 * [EOF]
 */
