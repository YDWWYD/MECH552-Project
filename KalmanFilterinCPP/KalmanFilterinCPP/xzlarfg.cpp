/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzlarfg.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xzlarfg.h"
#include "xscal.h"
#include "recip.h"
#include "xdlapy3.h"
#include "xnrm2.h"
#include "relop.h"
#include "butterBandpassOnly_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : creal_T *alpha1
 *                creal_T *x
 * Return Type  : creal_T
 */
creal_T b_xzlarfg(creal_T *alpha1, creal_T *x)
{
  creal_T tau;
  double xnorm;
  double beta1;
  int knt;
  double ai;
  creal_T b_alpha1;
  double x_re;
  double x_im;
  int k;
  tau.re = 0.0;
  tau.im = 0.0;
  xnorm = rt_hypotd_snf(x->re, x->im);
  if ((xnorm != 0.0) || (alpha1->im != 0.0)) {
    beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
    if (alpha1->re >= 0.0) {
      beta1 = -beta1;
    }

    if (fabs(beta1) < 1.0020841800044864E-292) {
      knt = 0;
      do {
        knt++;
        x->re *= 9.9792015476736E+291;
        x->im *= 9.9792015476736E+291;
        beta1 *= 9.9792015476736E+291;
        alpha1->re *= 9.9792015476736E+291;
        alpha1->im *= 9.9792015476736E+291;
      } while (!(fabs(beta1) >= 1.0020841800044864E-292));

      beta1 = xdlapy3(alpha1->re, alpha1->im, rt_hypotd_snf(x->re, x->im));
      if (alpha1->re >= 0.0) {
        beta1 = -beta1;
      }

      xnorm = beta1 - alpha1->re;
      ai = 0.0 - alpha1->im;
      if (ai == 0.0) {
        tau.re = xnorm / beta1;
        tau.im = 0.0;
      } else if (xnorm == 0.0) {
        tau.re = 0.0;
        tau.im = ai / beta1;
      } else {
        tau.re = xnorm / beta1;
        tau.im = ai / beta1;
      }

      b_alpha1.re = alpha1->re - beta1;
      b_alpha1.im = alpha1->im;
      *alpha1 = recip(b_alpha1);
      xnorm = alpha1->re;
      ai = alpha1->im;
      x_re = x->re;
      x_im = x->im;
      x->re = xnorm * x_re - ai * x_im;
      x->im = xnorm * x_im + ai * x_re;
      for (k = 1; k <= knt; k++) {
        beta1 *= 1.0020841800044864E-292;
      }

      alpha1->re = beta1;
      alpha1->im = 0.0;
    } else {
      xnorm = beta1 - alpha1->re;
      ai = 0.0 - alpha1->im;
      if (ai == 0.0) {
        tau.re = xnorm / beta1;
        tau.im = 0.0;
      } else if (xnorm == 0.0) {
        tau.re = 0.0;
        tau.im = ai / beta1;
      } else {
        tau.re = xnorm / beta1;
        tau.im = ai / beta1;
      }

      b_alpha1.re = alpha1->re - beta1;
      b_alpha1.im = alpha1->im;
      *alpha1 = recip(b_alpha1);
      xnorm = alpha1->re;
      ai = alpha1->im;
      x_re = x->re;
      x_im = x->im;
      x->re = xnorm * x_re - ai * x_im;
      x->im = xnorm * x_im + ai * x_re;
      alpha1->re = beta1;
      alpha1->im = 0.0;
    }
  }

  return tau;
}

/*
 * Arguments    : int n
 *                creal_T *alpha1
 *                emxArray_creal_T *x
 *                int ix0
 * Return Type  : creal_T
 */
creal_T xzlarfg(int n, creal_T *alpha1, emxArray_creal_T *x, int ix0)
{
  creal_T tau;
  double xnorm;
  double beta1;
  int knt;
  double ai;
  int i11;
  int k;
  creal_T b_alpha1;
  tau.re = 0.0;
  tau.im = 0.0;
  if (!(n <= 0)) {
    xnorm = xnrm2(n - 1, x, ix0);
    if ((xnorm != 0.0) || (alpha1->im != 0.0)) {
      beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
      if (alpha1->re >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = 0;
        i11 = (ix0 + n) - 2;
        do {
          knt++;
          for (k = ix0; k <= i11; k++) {
            xnorm = x->data[k - 1].re;
            ai = x->data[k - 1].im;
            x->data[k - 1].re = 9.9792015476736E+291 * xnorm - 0.0 * ai;
            x->data[k - 1].im = 9.9792015476736E+291 * ai + 0.0 * xnorm;
          }

          beta1 *= 9.9792015476736E+291;
          alpha1->re *= 9.9792015476736E+291;
          alpha1->im *= 9.9792015476736E+291;
        } while (!(fabs(beta1) >= 1.0020841800044864E-292));

        beta1 = xdlapy3(alpha1->re, alpha1->im, xnrm2(n - 1, x, ix0));
        if (alpha1->re >= 0.0) {
          beta1 = -beta1;
        }

        xnorm = beta1 - alpha1->re;
        ai = 0.0 - alpha1->im;
        if (ai == 0.0) {
          tau.re = xnorm / beta1;
          tau.im = 0.0;
        } else if (xnorm == 0.0) {
          tau.re = 0.0;
          tau.im = ai / beta1;
        } else {
          tau.re = xnorm / beta1;
          tau.im = ai / beta1;
        }

        b_alpha1.re = alpha1->re - beta1;
        b_alpha1.im = alpha1->im;
        *alpha1 = recip(b_alpha1);
        xscal(n - 1, *alpha1, x, ix0);
        for (k = 1; k <= knt; k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        alpha1->re = beta1;
        alpha1->im = 0.0;
      } else {
        xnorm = beta1 - alpha1->re;
        ai = 0.0 - alpha1->im;
        if (ai == 0.0) {
          tau.re = xnorm / beta1;
          tau.im = 0.0;
        } else if (xnorm == 0.0) {
          tau.re = 0.0;
          tau.im = ai / beta1;
        } else {
          tau.re = xnorm / beta1;
          tau.im = ai / beta1;
        }

        b_alpha1.re = alpha1->re - beta1;
        b_alpha1.im = alpha1->im;
        *alpha1 = recip(b_alpha1);
        xscal(n - 1, *alpha1, x, ix0);
        alpha1->re = beta1;
        alpha1->im = 0.0;
      }
    }
  }

  return tau;
}

/*
 * File trailer for xzlarfg.c
 *
 * [EOF]
 */
