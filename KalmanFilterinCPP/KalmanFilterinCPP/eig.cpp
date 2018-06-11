/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: eig.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "eig.h"
#include "butterBandpassOnly_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *alpha1
 *                const emxArray_creal_T *beta1
 *                emxArray_creal_T *D
 * Return Type  : void
 */
void makeD(const emxArray_creal_T *alpha1, const emxArray_creal_T *beta1,
           emxArray_creal_T *D)
{
  int i9;
  int loop_ub;
  double alpha1_re;
  double alpha1_im;
  double beta1_re;
  double beta1_im;
  double brm;
  double bim;
  double s;
  i9 = D->size[0];
  D->size[0] = alpha1->size[0];
  emxEnsureCapacity_creal_T1(D, i9);
  loop_ub = alpha1->size[0];
  for (i9 = 0; i9 < loop_ub; i9++) {
    alpha1_re = alpha1->data[i9].re;
    alpha1_im = alpha1->data[i9].im;
    beta1_re = beta1->data[i9].re;
    beta1_im = beta1->data[i9].im;
    if (beta1_im == 0.0) {
      if (alpha1_im == 0.0) {
        D->data[i9].re = alpha1_re / beta1_re;
        D->data[i9].im = 0.0;
      } else if (alpha1_re == 0.0) {
        D->data[i9].re = 0.0;
        D->data[i9].im = alpha1_im / beta1_re;
      } else {
        D->data[i9].re = alpha1_re / beta1_re;
        D->data[i9].im = alpha1_im / beta1_re;
      }
    } else if (beta1_re == 0.0) {
      if (alpha1_re == 0.0) {
        D->data[i9].re = alpha1_im / beta1_im;
        D->data[i9].im = 0.0;
      } else if (alpha1_im == 0.0) {
        D->data[i9].re = 0.0;
        D->data[i9].im = -(alpha1_re / beta1_im);
      } else {
        D->data[i9].re = alpha1_im / beta1_im;
        D->data[i9].im = -(alpha1_re / beta1_im);
      }
    } else {
      brm = fabs(beta1_re);
      bim = fabs(beta1_im);
      if (brm > bim) {
        s = beta1_im / beta1_re;
        bim = beta1_re + s * beta1_im;
        D->data[i9].re = (alpha1_re + s * alpha1_im) / bim;
        D->data[i9].im = (alpha1_im - s * alpha1_re) / bim;
      } else if (bim == brm) {
        if (beta1_re > 0.0) {
          s = 0.5;
        } else {
          s = -0.5;
        }

        if (beta1_im > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        D->data[i9].re = (alpha1_re * s + alpha1_im * bim) / brm;
        D->data[i9].im = (alpha1_im * s - alpha1_re * bim) / brm;
      } else {
        s = beta1_re / beta1_im;
        bim = beta1_im + s * beta1_re;
        D->data[i9].re = (s * alpha1_re + alpha1_im) / bim;
        D->data[i9].im = (s * alpha1_im - alpha1_re) / bim;
      }
    }
  }
}

/*
 * File trailer for eig.c
 *
 * [EOF]
 */
