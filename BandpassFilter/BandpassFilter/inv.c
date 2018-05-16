/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: inv.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "inv.h"
#include "butterBandpassOnly_emxutil.h"
#include "xtrsm.h"
#include "colon.h"
#include "xzgetrf.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *x
 *                emxArray_creal_T *y
 * Return Type  : void
 */
void inv(const emxArray_creal_T *x, emxArray_creal_T *y)
{
  int n;
  int pipk;
  int c;
  emxArray_creal_T *b_x;
  emxArray_int32_T *p;
  emxArray_int32_T *ipiv;
  int k;
  boolean_T b_y;
  int i;
  double y_re;
  double y_im;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    pipk = y->size[0] * y->size[1];
    y->size[0] = x->size[0];
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(y, pipk);
    c = x->size[0] * x->size[1];
    for (pipk = 0; pipk < c; pipk++) {
      y->data[pipk] = x->data[pipk];
    }
  } else {
    n = x->size[0];
    pipk = y->size[0] * y->size[1];
    y->size[0] = x->size[0];
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(y, pipk);
    c = x->size[0] * x->size[1];
    for (pipk = 0; pipk < c; pipk++) {
      y->data[pipk].re = 0.0;
      y->data[pipk].im = 0.0;
    }

    emxInit_creal_T(&b_x, 2);
    pipk = b_x->size[0] * b_x->size[1];
    b_x->size[0] = x->size[0];
    b_x->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(b_x, pipk);
    c = x->size[0] * x->size[1];
    for (pipk = 0; pipk < c; pipk++) {
      b_x->data[pipk] = x->data[pipk];
    }

    emxInit_int32_T1(&p, 2);
    emxInit_int32_T1(&ipiv, 2);
    xzgetrf(x->size[0], x->size[0], b_x, x->size[0], ipiv, &pipk);
    eml_signed_integer_colon(x->size[0], p);
    for (k = 0; k < ipiv->size[1]; k++) {
      if (ipiv->data[k] > 1 + k) {
        pipk = p->data[ipiv->data[k] - 1];
        p->data[ipiv->data[k] - 1] = p->data[k];
        p->data[k] = pipk;
      }
    }

    emxFree_int32_T(&ipiv);
    for (k = 0; k + 1 <= n; k++) {
      c = p->data[k] - 1;
      y->data[k + y->size[0] * (p->data[k] - 1)].re = 1.0;
      y->data[k + y->size[0] * (p->data[k] - 1)].im = 0.0;
      for (pipk = k; pipk + 1 <= n; pipk++) {
        b_y = ((y->data[pipk + y->size[0] * c].re != 0.0) || (y->data[pipk +
                y->size[0] * c].im != 0.0));
        if (b_y) {
          for (i = pipk + 1; i + 1 <= n; i++) {
            y_re = y->data[pipk + y->size[0] * c].re * b_x->data[i + b_x->size[0]
              * pipk].re - y->data[pipk + y->size[0] * c].im * b_x->data[i +
              b_x->size[0] * pipk].im;
            y_im = y->data[pipk + y->size[0] * c].re * b_x->data[i + b_x->size[0]
              * pipk].im + y->data[pipk + y->size[0] * c].im * b_x->data[i +
              b_x->size[0] * pipk].re;
            y->data[i + y->size[0] * c].re -= y_re;
            y->data[i + y->size[0] * c].im -= y_im;
          }
        }
      }
    }

    emxFree_int32_T(&p);
    xtrsm(x->size[0], x->size[0], b_x, x->size[0], y, x->size[0]);
    emxFree_creal_T(&b_x);
  }
}

/*
 * File trailer for inv.c
 *
 * [EOF]
 */
