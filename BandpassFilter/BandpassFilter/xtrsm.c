/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xtrsm.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xtrsm.h"
#include "recip.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                const emxArray_creal_T *A
 *                int lda
 *                emxArray_creal_T *B
 * Return Type  : void
 */
void b_xtrsm(int m, const emxArray_creal_T *A, int lda, emxArray_creal_T *B)
{
  int k;
  int kAcol;
  boolean_T b_B;
  int i;
  double B_re;
  double B_im;
  for (k = 0; k + 1 <= m; k++) {
    kAcol = lda * k;
    b_B = ((B->data[k].re != 0.0) || (B->data[k].im != 0.0));
    if (b_B) {
      for (i = k + 1; i + 1 <= m; i++) {
        B_re = B->data[k].re * A->data[i + kAcol].re - B->data[k].im * A->data[i
          + kAcol].im;
        B_im = B->data[k].re * A->data[i + kAcol].im + B->data[k].im * A->data[i
          + kAcol].re;
        B->data[i].re -= B_re;
        B->data[i].im -= B_im;
      }
    }
  }
}

/*
 * Arguments    : int m
 *                const emxArray_creal_T *A
 *                int lda
 *                emxArray_creal_T *B
 * Return Type  : void
 */
void c_xtrsm(int m, const emxArray_creal_T *A, int lda, emxArray_creal_T *B)
{
  int k;
  int kAcol;
  boolean_T b_B;
  double B_re;
  double B_im;
  double A_re;
  double A_im;
  double brm;
  double bim;
  double s;
  int i;
  for (k = m - 1; k + 1 > 0; k--) {
    kAcol = lda * k;
    b_B = ((B->data[k].re != 0.0) || (B->data[k].im != 0.0));
    if (b_B) {
      B_re = B->data[k].re;
      B_im = B->data[k].im;
      A_re = A->data[k + kAcol].re;
      A_im = A->data[k + kAcol].im;
      if (A_im == 0.0) {
        if (B_im == 0.0) {
          B->data[k].re = B_re / A_re;
          B->data[k].im = 0.0;
        } else if (B_re == 0.0) {
          B->data[k].re = 0.0;
          B->data[k].im = B_im / A_re;
        } else {
          B->data[k].re = B_re / A_re;
          B->data[k].im = B_im / A_re;
        }
      } else if (A_re == 0.0) {
        if (B_re == 0.0) {
          B->data[k].re = B_im / A_im;
          B->data[k].im = 0.0;
        } else if (B_im == 0.0) {
          B->data[k].re = 0.0;
          B->data[k].im = -(B_re / A_im);
        } else {
          B->data[k].re = B_im / A_im;
          B->data[k].im = -(B_re / A_im);
        }
      } else {
        brm = fabs(A_re);
        bim = fabs(A_im);
        if (brm > bim) {
          s = A_im / A_re;
          bim = A_re + s * A_im;
          B->data[k].re = (B_re + s * B_im) / bim;
          B->data[k].im = (B_im - s * B_re) / bim;
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

          B->data[k].re = (B_re * s + B_im * bim) / brm;
          B->data[k].im = (B_im * s - B_re * bim) / brm;
        } else {
          s = A_re / A_im;
          bim = A_im + s * A_re;
          B->data[k].re = (s * B_re + B_im) / bim;
          B->data[k].im = (s * B_im - B_re) / bim;
        }
      }

      for (i = 0; i + 1 <= k; i++) {
        B_re = B->data[k].re * A->data[i + kAcol].re - B->data[k].im * A->data[i
          + kAcol].im;
        B_im = B->data[k].re * A->data[i + kAcol].im + B->data[k].im * A->data[i
          + kAcol].re;
        B->data[i].re -= B_re;
        B->data[i].im -= B_im;
      }
    }
  }
}

/*
 * Arguments    : int n
 *                const emxArray_creal_T *A
 *                int lda
 *                emxArray_creal_T *B
 * Return Type  : void
 */
void d_xtrsm(int n, const emxArray_creal_T *A, int lda, emxArray_creal_T *B)
{
  int j;
  int jAcol;
  int k;
  creal_T b_A;
  boolean_T c_A;
  double B_re;
  double B_im;
  if ((n == 0) || (B->size[1] == 0)) {
  } else {
    for (j = 0; j + 1 <= n; j++) {
      jAcol = lda * j;
      for (k = 1; k <= j; k++) {
        c_A = ((A->data[(k + jAcol) - 1].re != 0.0) || (A->data[(k + jAcol) - 1]
                .im != 0.0));
        if (c_A) {
          B_re = A->data[(k + jAcol) - 1].re * B->data[k - 1].re - A->data[(k +
            jAcol) - 1].im * B->data[k - 1].im;
          B_im = A->data[(k + jAcol) - 1].re * B->data[k - 1].im + A->data[(k +
            jAcol) - 1].im * B->data[k - 1].re;
          B->data[j].re -= B_re;
          B->data[j].im -= B_im;
        }
      }

      b_A = A->data[j + jAcol];
      b_A = recip(b_A);
      B_re = B->data[j].re;
      B_im = B->data[j].im;
      B->data[j].re = b_A.re * B_re - b_A.im * B_im;
      B->data[j].im = b_A.re * B_im + b_A.im * B_re;
    }
  }
}

/*
 * Arguments    : int n
 *                const emxArray_creal_T *A
 *                int lda
 *                emxArray_creal_T *B
 * Return Type  : void
 */
void e_xtrsm(int n, const emxArray_creal_T *A, int lda, emxArray_creal_T *B)
{
  int j;
  int jAcol;
  int k;
  boolean_T b_A;
  double A_re;
  double A_im;
  if ((n == 0) || (B->size[1] == 0)) {
  } else {
    for (j = n - 1; j + 1 > 0; j--) {
      jAcol = lda * j - 1;
      for (k = j + 2; k <= n; k++) {
        b_A = ((A->data[k + jAcol].re != 0.0) || (A->data[k + jAcol].im != 0.0));
        if (b_A) {
          A_re = A->data[k + jAcol].re * B->data[k - 1].re - A->data[k + jAcol].
            im * B->data[k - 1].im;
          A_im = A->data[k + jAcol].re * B->data[k - 1].im + A->data[k + jAcol].
            im * B->data[k - 1].re;
          B->data[j].re -= A_re;
          B->data[j].im -= A_im;
        }
      }
    }
  }
}

/*
 * Arguments    : int m
 *                int n
 *                const emxArray_creal_T *A
 *                int lda
 *                emxArray_creal_T *B
 *                int ldb
 * Return Type  : void
 */
void xtrsm(int m, int n, const emxArray_creal_T *A, int lda, emxArray_creal_T *B,
           int ldb)
{
  int j;
  int jBcol;
  int k;
  int kAcol;
  boolean_T b_B;
  double B_re;
  double B_im;
  double A_re;
  double A_im;
  double brm;
  double bim;
  double s;
  int i;
  if ((n == 0) || ((B->size[0] == 0) || (B->size[1] == 0))) {
  } else {
    for (j = 1; j <= n; j++) {
      jBcol = ldb * (j - 1) - 1;
      for (k = m; k > 0; k--) {
        kAcol = lda * (k - 1) - 1;
        b_B = ((B->data[k + jBcol].re != 0.0) || (B->data[k + jBcol].im != 0.0));
        if (b_B) {
          B_re = B->data[k + jBcol].re;
          B_im = B->data[k + jBcol].im;
          A_re = A->data[k + kAcol].re;
          A_im = A->data[k + kAcol].im;
          if (A_im == 0.0) {
            if (B_im == 0.0) {
              B->data[k + jBcol].re = B_re / A_re;
              B->data[k + jBcol].im = 0.0;
            } else if (B_re == 0.0) {
              B->data[k + jBcol].re = 0.0;
              B->data[k + jBcol].im = B_im / A_re;
            } else {
              B->data[k + jBcol].re = B_re / A_re;
              B->data[k + jBcol].im = B_im / A_re;
            }
          } else if (A_re == 0.0) {
            if (B_re == 0.0) {
              B->data[k + jBcol].re = B_im / A_im;
              B->data[k + jBcol].im = 0.0;
            } else if (B_im == 0.0) {
              B->data[k + jBcol].re = 0.0;
              B->data[k + jBcol].im = -(B_re / A_im);
            } else {
              B->data[k + jBcol].re = B_im / A_im;
              B->data[k + jBcol].im = -(B_re / A_im);
            }
          } else {
            brm = fabs(A_re);
            bim = fabs(A_im);
            if (brm > bim) {
              s = A_im / A_re;
              bim = A_re + s * A_im;
              B->data[k + jBcol].re = (B_re + s * B_im) / bim;
              B->data[k + jBcol].im = (B_im - s * B_re) / bim;
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

              B->data[k + jBcol].re = (B_re * s + B_im * bim) / brm;
              B->data[k + jBcol].im = (B_im * s - B_re * bim) / brm;
            } else {
              s = A_re / A_im;
              bim = A_im + s * A_re;
              B->data[k + jBcol].re = (s * B_re + B_im) / bim;
              B->data[k + jBcol].im = (s * B_im - B_re) / bim;
            }
          }

          for (i = 1; i < k; i++) {
            B_re = B->data[k + jBcol].re * A->data[i + kAcol].re - B->data[k +
              jBcol].im * A->data[i + kAcol].im;
            B_im = B->data[k + jBcol].re * A->data[i + kAcol].im + B->data[k +
              jBcol].im * A->data[i + kAcol].re;
            B->data[i + jBcol].re -= B_re;
            B->data[i + jBcol].im -= B_im;
          }
        }
      }
    }
  }
}

/*
 * File trailer for xtrsm.c
 *
 * [EOF]
 */
