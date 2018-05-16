/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: sort1.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "sortButter.h"
#include "relop.h"
#include "butterBandpassOnly_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_creal_T *x
 * Return Type  : void
 */
void sortButter(emxArray_creal_T *x)
{
  int dim;
  int i19;
  emxArray_creal_T *vwork;
  int i2;
  int vstride;
  int k;
  int j;
  emxArray_int32_T *iidx;
  emxArray_int32_T *iwork;
  emxArray_creal_T *xwork;
  int n;
  unsigned int unnamed_idx_0;
  creal_T b_vwork;
  creal_T c_vwork;
  boolean_T p;
  int b_j;
  int pEnd;
  int b_p;
  int q;
  int qEnd;
  int kEnd;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  if (dim <= 1) {
    i19 = x->size[0];
  } else {
    i19 = 1;
  }

  emxInit_creal_T1(&vwork, 1);
  i2 = vwork->size[0];
  vwork->size[0] = i19;
  emxEnsureCapacity_creal_T1(vwork, i2);
  vstride = 1;
  k = 1;
  while (k <= dim - 1) {
    vstride *= x->size[0];
    k = 2;
  }

  j = 0;
  emxInit_int32_T(&iidx, 1);
  emxInit_int32_T(&iwork, 1);
  emxInit_creal_T1(&xwork, 1);
  while (j + 1 <= vstride) {
    for (k = 0; k + 1 <= i19; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    n = vwork->size[0];
    unnamed_idx_0 = (unsigned int)vwork->size[0];
    i2 = iidx->size[0];
    iidx->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity_int32_T(iidx, i2);
    dim = (int)unnamed_idx_0;
    for (i2 = 0; i2 < dim; i2++) {
      iidx->data[i2] = 0;
    }

    if (vwork->size[0] != 0) {
      i2 = iwork->size[0];
      iwork->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity_int32_T(iwork, i2);
      for (k = 1; k <= n - 1; k += 2) {
        b_vwork = vwork->data[k - 1];
        c_vwork = vwork->data[k];
        if (relop(b_vwork, c_vwork) || (rtIsNaN(vwork->data[k].re) || rtIsNaN
             (vwork->data[k].im))) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          iidx->data[k - 1] = k;
          iidx->data[k] = k + 1;
        } else {
          iidx->data[k - 1] = k + 1;
          iidx->data[k] = k;
        }
      }

      if ((vwork->size[0] & 1) != 0) {
        iidx->data[vwork->size[0] - 1] = vwork->size[0];
      }

      dim = 2;
      while (dim < n) {
        i2 = dim << 1;
        b_j = 1;
        for (pEnd = 1 + dim; pEnd < n + 1; pEnd = qEnd + dim) {
          b_p = b_j;
          q = pEnd - 1;
          qEnd = b_j + i2;
          if (qEnd > n + 1) {
            qEnd = n + 1;
          }

          k = 0;
          kEnd = qEnd - b_j;
          while (k + 1 <= kEnd) {
            b_vwork = vwork->data[iidx->data[b_p - 1] - 1];
            c_vwork = vwork->data[iidx->data[q] - 1];
            if (relop(b_vwork, c_vwork) || (rtIsNaN(vwork->data[iidx->data[q] -
                  1].re) || rtIsNaN(vwork->data[iidx->data[q] - 1].im))) {
              p = true;
            } else {
              p = false;
            }

            if (p) {
              iwork->data[k] = iidx->data[b_p - 1];
              b_p++;
              if (b_p == pEnd) {
                while (q + 1 < qEnd) {
                  k++;
                  iwork->data[k] = iidx->data[q];
                  q++;
                }
              }
            } else {
              iwork->data[k] = iidx->data[q];
              q++;
              if (q + 1 == qEnd) {
                while (b_p < pEnd) {
                  k++;
                  iwork->data[k] = iidx->data[b_p - 1];
                  b_p++;
                }
              }
            }

            k++;
          }

          for (k = 0; k + 1 <= kEnd; k++) {
            iidx->data[(b_j + k) - 1] = iwork->data[k];
          }

          b_j = qEnd;
        }

        dim = i2;
      }

      unnamed_idx_0 = (unsigned int)vwork->size[0];
      i2 = xwork->size[0];
      xwork->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity_creal_T1(xwork, i2);
      for (k = 0; k + 1 <= n; k++) {
        xwork->data[k] = vwork->data[k];
      }

      for (k = 0; k + 1 <= n; k++) {
        vwork->data[k] = xwork->data[iidx->data[k] - 1];
      }
    }

    for (k = 0; k + 1 <= i19; k++) {
      x->data[j + k * vstride] = vwork->data[k];
    }

    j++;
  }

  emxFree_creal_T(&xwork);
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&iidx);
  emxFree_creal_T(&vwork);
}

/*
 * File trailer for sort1.c
 *
 * [EOF]
 */
