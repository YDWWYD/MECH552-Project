/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgeqp3.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xgeqp3.h"
#include "xnrm2.h"
#include "relop.h"
#include "xgerc.h"
#include "xgemv.h"
#include "xzlarfg.h"
#include "xscal.h"
#include "recip.h"
#include "xdlapy3.h"
#include "ixamax.h"
#include "butterBandpassOnly_emxutil.h"
#include "colon.h"
#include "butterBandpassOnly_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_creal_T *A
 *                emxArray_creal_T *tau
 *                emxArray_int32_T *jpvt
 * Return Type  : void
 */
void xgeqp3(emxArray_creal_T *A, emxArray_creal_T *tau, emxArray_int32_T *jpvt)
{
  int m;
  int n;
  int ix;
  int mn;
  int iy;
  emxArray_creal_T *work;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int k;
  int i;
  int i_i;
  int nmi;
  int mmi;
  int pvt;
  creal_T atmp;
  double temp2;
  double temp1;
  creal_T temp;
  double beta1;
  double ai;
  int lastv;
  boolean_T exitg1;
  creal_T b_temp;
  boolean_T b_A;
  int exitg2;
  m = A->size[0];
  n = A->size[1];
  ix = A->size[0];
  mn = A->size[1];
  if (ix < mn) {
    mn = ix;
  }

  iy = tau->size[0];
  tau->size[0] = mn;
  emxEnsureCapacity_creal_T1(tau, iy);
  eml_signed_integer_colon(A->size[1], jpvt);
  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    emxInit_creal_T1(&work, 1);
    ix = A->size[1];
    iy = work->size[0];
    work->size[0] = ix;
    emxEnsureCapacity_creal_T1(work, iy);
    for (iy = 0; iy < ix; iy++) {
      work->data[iy].re = 0.0;
      work->data[iy].im = 0.0;
    }

    emxInit_real_T(&vn1, 1);
    emxInit_real_T(&vn2, 1);
    ix = A->size[1];
    iy = vn1->size[0];
    vn1->size[0] = ix;
    emxEnsureCapacity_real_T1(vn1, iy);
    iy = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity_real_T1(vn2, iy);
    k = 1;
    for (ix = 0; ix + 1 <= n; ix++) {
      vn1->data[ix] = xnrm2(m, A, k);
      vn2->data[ix] = vn1->data[ix];
      k += m;
    }

    for (i = 0; i + 1 <= mn; i++) {
      i_i = i + i * m;
      nmi = n - i;
      mmi = m - i;
      pvt = (i + ixamax(nmi, vn1, i + 1)) - 1;
      if (pvt + 1 != i + 1) {
        ix = m * pvt;
        iy = m * i;
        for (k = 1; k <= m; k++) {
          temp = A->data[ix];
          A->data[ix] = A->data[iy];
          A->data[iy] = temp;
          ix++;
          iy++;
        }

        ix = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i];
        jpvt->data[i] = ix;
        vn1->data[pvt] = vn1->data[i];
        vn2->data[pvt] = vn2->data[i];
      }

      if (i + 1 < m) {
        atmp = A->data[i_i];
        tau->data[i] = xzlarfg(mmi, &atmp, A, i_i + 2);
        A->data[i_i] = atmp;
      } else {
        atmp = A->data[i_i];
        temp2 = 0.0;
        temp1 = 0.0;
        if (A->data[i_i].im != 0.0) {
          beta1 = xdlapy3(A->data[i_i].re, A->data[i_i].im, 0.0);
          if (A->data[i_i].re >= 0.0) {
            beta1 = -beta1;
          }

          if (fabs(beta1) < 1.0020841800044864E-292) {
            iy = 0;
            do {
              iy++;
              for (k = i_i; k + 1 <= i_i; k++) {
                temp1 = A->data[k].re;
                temp2 = A->data[k].im;
                A->data[k].re = 9.9792015476736E+291 * temp1 - 0.0 * temp2;
                A->data[k].im = 9.9792015476736E+291 * temp2 + 0.0 * temp1;
              }

              beta1 *= 9.9792015476736E+291;
              atmp.re *= 9.9792015476736E+291;
              atmp.im *= 9.9792015476736E+291;
            } while (!(fabs(beta1) >= 1.0020841800044864E-292));

            beta1 = xdlapy3(atmp.re, atmp.im, 0.0);
            if (atmp.re >= 0.0) {
              beta1 = -beta1;
            }

            temp1 = beta1 - atmp.re;
            if (0.0 - atmp.im == 0.0) {
              temp2 = temp1 / beta1;
              temp1 = 0.0;
            } else if (temp1 == 0.0) {
              temp2 = 0.0;
              temp1 = (0.0 - atmp.im) / beta1;
            } else {
              temp2 = temp1 / beta1;
              temp1 = (0.0 - atmp.im) / beta1;
            }

            temp.re = atmp.re - beta1;
            temp.im = atmp.im;
            b_xscal(recip(temp), A, i_i + 1);
            for (k = 1; k <= iy; k++) {
              beta1 *= 1.0020841800044864E-292;
            }

            atmp.re = beta1;
            atmp.im = 0.0;
          } else {
            temp1 = beta1 - A->data[i_i].re;
            ai = 0.0 - A->data[i_i].im;
            if (ai == 0.0) {
              temp2 = temp1 / beta1;
              temp1 = 0.0;
            } else if (temp1 == 0.0) {
              temp2 = 0.0;
              temp1 = ai / beta1;
            } else {
              temp2 = temp1 / beta1;
              temp1 = ai / beta1;
            }

            temp.re = A->data[i_i].re - beta1;
            temp.im = A->data[i_i].im;
            b_xscal(recip(temp), A, i_i + 1);
            atmp.re = beta1;
            atmp.im = 0.0;
          }
        }

        tau->data[i].re = temp2;
        tau->data[i].im = temp1;
        A->data[i_i] = atmp;
      }

      if (i + 1 < n) {
        atmp = A->data[i_i];
        A->data[i_i].re = 1.0;
        A->data[i_i].im = 0.0;
        k = (i + (i + 1) * m) + 1;
        temp.re = tau->data[i].re;
        temp.im = -tau->data[i].im;
        if ((temp.re != 0.0) || (temp.im != 0.0)) {
          lastv = mmi;
          ix = (i_i + mmi) - 1;
          exitg1 = false;
          while ((!exitg1) && (lastv > 0)) {
            b_A = ((A->data[ix].re == 0.0) && (A->data[ix].im == 0.0));
            if (b_A) {
              lastv--;
              ix--;
            } else {
              exitg1 = true;
            }
          }

          ix = nmi - 1;
          exitg1 = false;
          while ((!exitg1) && (ix > 0)) {
            iy = k + (ix - 1) * m;
            pvt = iy;
            do {
              exitg2 = 0;
              if (pvt <= (iy + lastv) - 1) {
                b_A = ((A->data[pvt - 1].re != 0.0) || (A->data[pvt - 1].im !=
                        0.0));
                if (b_A) {
                  exitg2 = 1;
                } else {
                  pvt++;
                }
              } else {
                ix--;
                exitg2 = 2;
              }
            } while (exitg2 == 0);

            if (exitg2 == 1) {
              exitg1 = true;
            }
          }
        } else {
          lastv = 0;
          ix = 0;
        }

        if (lastv > 0) {
          xgemv(lastv, ix, A, k, m, A, i_i + 1, work);
          b_temp.re = -temp.re;
          b_temp.im = -temp.im;
          xgerc(lastv, ix, b_temp, i_i + 1, work, A, k, m);
        }

        A->data[i_i] = atmp;
      }

      for (ix = i + 1; ix + 1 <= n; ix++) {
        if (vn1->data[ix] != 0.0) {
          temp1 = rt_hypotd_snf(A->data[i + A->size[0] * ix].re, A->data[i +
                                A->size[0] * ix].im) / vn1->data[ix];
          temp1 = 1.0 - temp1 * temp1;
          if (temp1 < 0.0) {
            temp1 = 0.0;
          }

          temp2 = vn1->data[ix] / vn2->data[ix];
          temp2 = temp1 * (temp2 * temp2);
          if (temp2 <= 1.4901161193847656E-8) {
            if (i + 1 < m) {
              vn1->data[ix] = xnrm2(mmi - 1, A, (i + m * ix) + 2);
              vn2->data[ix] = vn1->data[ix];
            } else {
              vn1->data[ix] = 0.0;
              vn2->data[ix] = 0.0;
            }
          } else {
            vn1->data[ix] *= sqrt(temp1);
          }
        }
      }
    }

    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
    emxFree_creal_T(&work);
  }
}

/*
 * File trailer for xgeqp3.c
 *
 * [EOF]
 */
