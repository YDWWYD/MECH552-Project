/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzgetrf.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xzgetrf.h"
#include "colon.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                emxArray_creal_T *A
 *                int lda
 *                emxArray_int32_T *ipiv
 *                int *info
 * Return Type  : void
 */
void xzgetrf(int m, int n, emxArray_creal_T *A, int lda, emxArray_int32_T *ipiv,
             int *info)
{
  int iy;
  int b_info;
  int u0;
  int j;
  int mmj;
  boolean_T b_A;
  int c;
  int ix;
  double smax;
  int jA;
  int i10;
  double s;
  int jy;
  int b_j;
  double A_re;
  double A_im;
  double b_A_re;
  double b_A_im;
  double brm;
  int ijA;
  if (m < n) {
    iy = m;
  } else {
    iy = n;
  }

  eml_signed_integer_colon(iy, ipiv);
  b_info = 0;
  if ((m < 1) || (n < 1)) {
  } else {
    u0 = m - 1;
    if (!(u0 < n)) {
      u0 = n;
    }

    for (j = 0; j + 1 <= u0; j++) {
      mmj = m - j;
      c = j * (lda + 1);
      if (mmj < 1) {
        iy = -1;
      } else {
        iy = 0;
        if (mmj > 1) {
          ix = c;
          smax = fabs(A->data[c].re) + fabs(A->data[c].im);
          for (jA = 1; jA + 1 <= mmj; jA++) {
            ix++;
            s = fabs(A->data[ix].re) + fabs(A->data[ix].im);
            if (s > smax) {
              iy = jA;
              smax = s;
            }
          }
        }
      }

      b_A = ((A->data[c + iy].re != 0.0) || (A->data[c + iy].im != 0.0));
      if (b_A) {
        if (iy != 0) {
          ipiv->data[j] = (j + iy) + 1;
          ix = j;
          iy += j;
          for (jA = 1; jA <= n; jA++) {
            smax = A->data[ix].re;
            s = A->data[ix].im;
            A->data[ix] = A->data[iy];
            A->data[iy].re = smax;
            A->data[iy].im = s;
            ix += lda;
            iy += lda;
          }
        }

        i10 = c + mmj;
        for (iy = c + 1; iy + 1 <= i10; iy++) {
          A_re = A->data[iy].re;
          A_im = A->data[iy].im;
          b_A_re = A->data[c].re;
          b_A_im = A->data[c].im;
          if (b_A_im == 0.0) {
            if (A_im == 0.0) {
              A->data[iy].re = A_re / b_A_re;
              A->data[iy].im = 0.0;
            } else if (A_re == 0.0) {
              A->data[iy].re = 0.0;
              A->data[iy].im = A_im / b_A_re;
            } else {
              A->data[iy].re = A_re / b_A_re;
              A->data[iy].im = A_im / b_A_re;
            }
          } else if (b_A_re == 0.0) {
            if (A_re == 0.0) {
              A->data[iy].re = A_im / b_A_im;
              A->data[iy].im = 0.0;
            } else if (A_im == 0.0) {
              A->data[iy].re = 0.0;
              A->data[iy].im = -(A_re / b_A_im);
            } else {
              A->data[iy].re = A_im / b_A_im;
              A->data[iy].im = -(A_re / b_A_im);
            }
          } else {
            brm = fabs(b_A_re);
            smax = fabs(b_A_im);
            if (brm > smax) {
              s = b_A_im / b_A_re;
              smax = b_A_re + s * b_A_im;
              A->data[iy].re = (A_re + s * A_im) / smax;
              A->data[iy].im = (A_im - s * A_re) / smax;
            } else if (smax == brm) {
              if (b_A_re > 0.0) {
                s = 0.5;
              } else {
                s = -0.5;
              }

              if (b_A_im > 0.0) {
                smax = 0.5;
              } else {
                smax = -0.5;
              }

              A->data[iy].re = (A_re * s + A_im * smax) / brm;
              A->data[iy].im = (A_im * s - A_re * smax) / brm;
            } else {
              s = b_A_re / b_A_im;
              smax = b_A_im + s * b_A_re;
              A->data[iy].re = (s * A_re + A_im) / smax;
              A->data[iy].im = (s * A_im - A_re) / smax;
            }
          }
        }
      } else {
        b_info = j + 1;
      }

      iy = (n - j) - 1;
      jA = c + lda;
      jy = c + lda;
      for (b_j = 1; b_j <= iy; b_j++) {
        b_A = ((A->data[jy].re != 0.0) || (A->data[jy].im != 0.0));
        if (b_A) {
          smax = -A->data[jy].re - A->data[jy].im * 0.0;
          s = A->data[jy].re * 0.0 + -A->data[jy].im;
          ix = c + 1;
          i10 = mmj + jA;
          for (ijA = 1 + jA; ijA + 1 <= i10; ijA++) {
            A_re = A->data[ix].re * smax - A->data[ix].im * s;
            A_im = A->data[ix].re * s + A->data[ix].im * smax;
            A->data[ijA].re += A_re;
            A->data[ijA].im += A_im;
            ix++;
          }
        }

        jy += lda;
        jA += lda;
      }
    }

    if ((b_info == 0) && (m <= n)) {
      b_A = ((A->data[(m + A->size[0] * (m - 1)) - 1].re != 0.0) || (A->data[(m
               + A->size[0] * (m - 1)) - 1].im != 0.0));
      if (!b_A) {
        b_info = m;
      }
    }
  }

  *info = b_info;
}

/*
 * File trailer for xzgetrf.c
 *
 * [EOF]
 */
