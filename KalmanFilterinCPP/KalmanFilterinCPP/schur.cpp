/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: schur.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "schur.h"
#include "butterBandpassOnly_emxutil.h"
#include "xgerc.h"
#include "xgemv.h"
#include "xzlarfg.h"
#include "xzhseqr.h"
#include "anyNonFinite.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_creal_T *A
 *                emxArray_creal_T *V
 * Return Type  : void
 */
void schur(emxArray_creal_T *A, emxArray_creal_T *V)
{
  int n;
  int i8;
  int ntau;
  unsigned int uv0[2];
  emxArray_creal_T *tau;
  emxArray_creal_T *work;
  int istart;
  int jend;
  int i;
  int j;
  int im1n;
  int in;
  creal_T alpha1;
  int jy;
  boolean_T b_tau;
  int lastv;
  int lastc;
  boolean_T exitg1;
  double c_re;
  double c_im;
  int ix;
  int exitg2;
  double temp_re;
  creal_T c;
  double temp_im;
  double work_re;
  double work_im;
  if (anyNonFinite(A)) {
    for (i8 = 0; i8 < 2; i8++) {
      uv0[i8] = (unsigned int)A->size[i8];
    }

    i8 = V->size[0] * V->size[1];
    V->size[0] = (int)uv0[0];
    V->size[1] = (int)uv0[1];
    emxEnsureCapacity_creal_T(V, i8);
    ntau = (int)uv0[0] * (int)uv0[1];
    for (i8 = 0; i8 < ntau; i8++) {
      V->data[i8].re = rtNaN;
      V->data[i8].im = 0.0;
    }

    ntau = V->size[0];
    if ((V->size[0] == 0) || (V->size[1] == 0) || (1 >= V->size[0])) {
    } else {
      istart = 2;
      if (V->size[0] - 2 < V->size[1] - 1) {
        jend = V->size[0] - 1;
      } else {
        jend = V->size[1];
      }

      for (j = 1; j <= jend; j++) {
        for (i = istart; i <= ntau; i++) {
          V->data[(i + V->size[0] * (j - 1)) - 1].re = 0.0;
          V->data[(i + V->size[0] * (j - 1)) - 1].im = 0.0;
        }

        istart++;
      }
    }
  } else {
    n = A->size[0];
    if (A->size[0] < 1) {
      ntau = 0;
    } else {
      ntau = A->size[0] - 1;
    }

    emxInit_creal_T1(&tau, 1);
    emxInit_creal_T1(&work, 1);
    i8 = tau->size[0];
    tau->size[0] = ntau;
    emxEnsureCapacity_creal_T1(tau, i8);
    ntau = A->size[0];
    i8 = work->size[0];
    work->size[0] = ntau;
    emxEnsureCapacity_creal_T1(work, i8);
    for (i8 = 0; i8 < ntau; i8++) {
      work->data[i8].re = 0.0;
      work->data[i8].im = 0.0;
    }

    for (i = 0; i + 1 < n; i++) {
      im1n = i * n + 2;
      in = (i + 1) * n;
      alpha1 = A->data[(i + A->size[0] * i) + 1];
      ntau = i + 3;
      if (!(ntau < n)) {
        ntau = n;
      }

      tau->data[i] = xzlarfg((n - i) - 1, &alpha1, A, ntau + i * n);
      A->data[(i + A->size[0] * i) + 1].re = 1.0;
      A->data[(i + A->size[0] * i) + 1].im = 0.0;
      ntau = (n - i) - 3;
      jy = (i + im1n) - 1;
      b_tau = ((tau->data[i].re != 0.0) || (tau->data[i].im != 0.0));
      if (b_tau) {
        lastv = ntau + 2;
        ntau += jy;
        exitg1 = false;
        while ((!exitg1) && (lastv > 0)) {
          b_tau = ((A->data[ntau + 1].re == 0.0) && (A->data[ntau + 1].im == 0.0));
          if (b_tau) {
            lastv--;
            ntau--;
          } else {
            exitg1 = true;
          }
        }

        lastc = n;
        exitg1 = false;
        while ((!exitg1) && (lastc > 0)) {
          ntau = in + lastc;
          j = ntau;
          do {
            exitg2 = 0;
            if ((n > 0) && (j <= ntau + (lastv - 1) * n)) {
              b_tau = ((A->data[j - 1].re != 0.0) || (A->data[j - 1].im != 0.0));
              if (b_tau) {
                exitg2 = 1;
              } else {
                j += n;
              }
            } else {
              lastc--;
              exitg2 = 2;
            }
          } while (exitg2 == 0);

          if (exitg2 == 1) {
            exitg1 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        if (lastc != 0) {
          for (ntau = 1; ntau <= lastc; ntau++) {
            work->data[ntau - 1].re = 0.0;
            work->data[ntau - 1].im = 0.0;
          }

          ix = jy;
          i8 = (in + n * (lastv - 1)) + 1;
          istart = in + 1;
          while ((n > 0) && (istart <= i8)) {
            c_re = A->data[ix].re - 0.0 * A->data[ix].im;
            c_im = A->data[ix].im + 0.0 * A->data[ix].re;
            ntau = 0;
            jend = (istart + lastc) - 1;
            for (j = istart; j <= jend; j++) {
              temp_re = A->data[j - 1].re * c_re - A->data[j - 1].im * c_im;
              temp_im = A->data[j - 1].re * c_im + A->data[j - 1].im * c_re;
              work->data[ntau].re += temp_re;
              work->data[ntau].im += temp_im;
              ntau++;
            }

            ix++;
            istart += n;
          }
        }

        c_re = -tau->data[i].re;
        c_im = -tau->data[i].im;
        if (!((c_re == 0.0) && (c_im == 0.0))) {
          ntau = in;
          for (j = 1; j <= lastv; j++) {
            b_tau = ((A->data[jy].re != 0.0) || (A->data[jy].im != 0.0));
            if (b_tau) {
              temp_re = A->data[jy].re * c_re + A->data[jy].im * c_im;
              temp_im = A->data[jy].re * c_im - A->data[jy].im * c_re;
              ix = 0;
              i8 = lastc + ntau;
              for (istart = ntau; istart + 1 <= i8; istart++) {
                work_re = work->data[ix].re * temp_re - work->data[ix].im *
                  temp_im;
                work_im = work->data[ix].re * temp_im + work->data[ix].im *
                  temp_re;
                A->data[istart].re += work_re;
                A->data[istart].im += work_im;
                ix++;
              }
            }

            jy++;
            ntau += n;
          }
        }
      }

      ntau = (n - i) - 3;
      jend = i + im1n;
      istart = (i + in) + 2;
      c_re = tau->data[i].re;
      c_im = -tau->data[i].im;
      if ((c_re != 0.0) || (c_im != 0.0)) {
        lastv = ntau + 2;
        ntau += jend;
        exitg1 = false;
        while ((!exitg1) && (lastv > 0)) {
          b_tau = ((A->data[ntau].re == 0.0) && (A->data[ntau].im == 0.0));
          if (b_tau) {
            lastv--;
            ntau--;
          } else {
            exitg1 = true;
          }
        }

        lastc = (n - i) - 1;
        exitg1 = false;
        while ((!exitg1) && (lastc > 0)) {
          ntau = istart + (lastc - 1) * n;
          j = ntau;
          do {
            exitg2 = 0;
            if (j <= (ntau + lastv) - 1) {
              b_tau = ((A->data[j - 1].re != 0.0) || (A->data[j - 1].im != 0.0));
              if (b_tau) {
                exitg2 = 1;
              } else {
                j++;
              }
            } else {
              lastc--;
              exitg2 = 2;
            }
          } while (exitg2 == 0);

          if (exitg2 == 1) {
            exitg1 = true;
          }
        }
      } else {
        lastv = 0;
        lastc = 0;
      }

      if (lastv > 0) {
        xgemv(lastv, lastc, A, istart, n, A, jend, work);
        c.re = -c_re;
        c.im = -c_im;
        xgerc(lastv, lastc, c, jend, work, A, istart, n);
      }

      A->data[(i + A->size[0] * i) + 1] = alpha1;
    }

    emxFree_creal_T(&work);
    emxFree_creal_T(&tau);
    i8 = V->size[0] * V->size[1];
    V->size[0] = A->size[0];
    V->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(V, i8);
    ntau = A->size[0] * A->size[1];
    for (i8 = 0; i8 < ntau; i8++) {
      V->data[i8] = A->data[i8];
    }

    eml_zlahqr(V);
    ntau = V->size[0];
    if ((V->size[0] == 0) || (V->size[1] == 0) || (3 >= V->size[0])) {
    } else {
      istart = 4;
      if (V->size[0] - 4 < V->size[1] - 1) {
        jend = V->size[0] - 3;
      } else {
        jend = V->size[1];
      }

      for (j = 1; j <= jend; j++) {
        for (i = istart; i <= ntau; i++) {
          V->data[(i + V->size[0] * (j - 1)) - 1].re = 0.0;
          V->data[(i + V->size[0] * (j - 1)) - 1].im = 0.0;
        }

        istart++;
      }
    }
  }
}

/*
 * File trailer for schur.c
 *
 * [EOF]
 */
