/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xunormqr.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xunormqr.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *Q
 *                emxArray_creal_T *C
 *                const emxArray_creal_T *tau
 * Return Type  : void
 */
void xunormqr(const emxArray_creal_T *Q, emxArray_creal_T *C, const
              emxArray_creal_T *tau)
{
  int m;
  int u0;
  int mn;
  double tauj_re;
  double tauj_im;
  double wj_re;
  double wj_im;
  int i;
  double b_wj_re;
  m = Q->size[0];
  u0 = Q->size[0];
  mn = Q->size[1];
  if (u0 < mn) {
    mn = u0;
  }

  for (u0 = 0; u0 + 1 <= mn; u0++) {
    tauj_re = tau->data[u0].re;
    tauj_im = -tau->data[u0].im;
    if ((tauj_re != 0.0) || (tauj_im != 0.0)) {
      wj_re = C->data[u0].re;
      wj_im = C->data[u0].im;
      for (i = u0 + 1; i + 1 <= m; i++) {
        wj_re += Q->data[i + Q->size[0] * u0].re * C->data[i].re + Q->data[i +
          Q->size[0] * u0].im * C->data[i].im;
        wj_im += Q->data[i + Q->size[0] * u0].re * C->data[i].im - Q->data[i +
          Q->size[0] * u0].im * C->data[i].re;
      }

      b_wj_re = wj_re;
      wj_re = tauj_re * wj_re - tauj_im * wj_im;
      wj_im = tauj_re * wj_im + tauj_im * b_wj_re;
      if ((wj_re != 0.0) || (wj_im != 0.0)) {
        C->data[u0].re -= wj_re;
        C->data[u0].im -= wj_im;
        for (i = u0 + 1; i + 1 <= m; i++) {
          tauj_re = Q->data[i + Q->size[0] * u0].re * wj_re - Q->data[i +
            Q->size[0] * u0].im * wj_im;
          tauj_im = Q->data[i + Q->size[0] * u0].re * wj_im + Q->data[i +
            Q->size[0] * u0].im * wj_re;
          C->data[i].re -= tauj_re;
          C->data[i].im -= tauj_im;
        }
      }
    }
  }
}

/*
 * File trailer for xunormqr.c
 *
 * [EOF]
 */
