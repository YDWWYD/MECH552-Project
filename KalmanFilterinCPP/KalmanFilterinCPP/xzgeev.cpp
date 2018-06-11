/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzgeev.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xzgeev.h"
#include "butterBandpassOnly_emxutil.h"
#include "xzlartg.h"
#include "xzlascl.h"
#include "xzhgeqz.h"
#include "relop.h"
#include "butterBandpassOnly_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *A
 *                int *info
 *                emxArray_creal_T *alpha1
 *                emxArray_creal_T *beta1
 * Return Type  : void
 */
void xzgeev(const emxArray_creal_T *A, int *info, emxArray_creal_T *alpha1,
            emxArray_creal_T *beta1)
{
  emxArray_creal_T *At;
  int jcol;
  int ii;
  int nzcount;
  double anrm;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  int ilo;
  double ctoc;
  int ihi;
  boolean_T notdone;
  int exitg3;
  double cfrom1;
  int i;
  double cto1;
  int j;
  double mul;
  creal_T b_At;
  creal_T c_At;
  double c;
  creal_T atmp;
  boolean_T exitg4;
  int exitg2;
  boolean_T d_At;
  double stemp_re;
  emxInit_creal_T(&At, 2);
  jcol = At->size[0] * At->size[1];
  At->size[0] = A->size[0];
  At->size[1] = A->size[1];
  emxEnsureCapacity_creal_T(At, jcol);
  ii = A->size[0] * A->size[1];
  for (jcol = 0; jcol < ii; jcol++) {
    At->data[jcol] = A->data[jcol];
  }

  nzcount = 0;
  jcol = alpha1->size[0];
  alpha1->size[0] = At->size[0];
  emxEnsureCapacity_creal_T1(alpha1, jcol);
  ii = At->size[0];
  for (jcol = 0; jcol < ii; jcol++) {
    alpha1->data[jcol].re = 0.0;
    alpha1->data[jcol].im = 0.0;
  }

  jcol = beta1->size[0];
  beta1->size[0] = At->size[0];
  emxEnsureCapacity_creal_T1(beta1, jcol);
  ii = At->size[0];
  for (jcol = 0; jcol < ii; jcol++) {
    beta1->data[jcol].re = 0.0;
    beta1->data[jcol].im = 0.0;
  }

  if (!((At->size[0] == 0) || (At->size[1] == 0))) {
    anrm = 0.0;
    jcol = 0;
    exitg1 = false;
    while ((!exitg1) && (jcol <= At->size[0] * At->size[1] - 1)) {
      absxk = rt_hypotd_snf(At->data[jcol].re, At->data[jcol].im);
      if (rtIsNaN(absxk)) {
        anrm = rtNaN;
        exitg1 = true;
      } else {
        if (absxk > anrm) {
          anrm = absxk;
        }

        jcol++;
      }
    }

    if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
      jcol = alpha1->size[0];
      alpha1->size[0] = At->size[0];
      emxEnsureCapacity_creal_T1(alpha1, jcol);
      ii = At->size[0];
      for (jcol = 0; jcol < ii; jcol++) {
        alpha1->data[jcol].re = rtNaN;
        alpha1->data[jcol].im = 0.0;
      }

      jcol = beta1->size[0];
      beta1->size[0] = At->size[0];
      emxEnsureCapacity_creal_T1(beta1, jcol);
      ii = At->size[0];
      for (jcol = 0; jcol < ii; jcol++) {
        beta1->data[jcol].re = rtNaN;
        beta1->data[jcol].im = 0.0;
      }
    } else {
      ilascl = false;
      anrmto = anrm;
      if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
        anrmto = 6.7178761075670888E-139;
        ilascl = true;
      } else {
        if (anrm > 1.4885657073574029E+138) {
          anrmto = 1.4885657073574029E+138;
          ilascl = true;
        }
      }

      if (ilascl) {
        absxk = anrm;
        ctoc = anrmto;
        notdone = true;
        while (notdone) {
          cfrom1 = absxk * 2.0041683600089728E-292;
          cto1 = ctoc / 4.9896007738368E+291;
          if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
            mul = 2.0041683600089728E-292;
            absxk = cfrom1;
          } else if (cto1 > absxk) {
            mul = 4.9896007738368E+291;
            ctoc = cto1;
          } else {
            mul = ctoc / absxk;
            notdone = false;
          }

          jcol = At->size[0] * At->size[1];
          emxEnsureCapacity_creal_T(At, jcol);
          jcol = At->size[0];
          ii = At->size[1];
          ii *= jcol;
          for (jcol = 0; jcol < ii; jcol++) {
            At->data[jcol].re *= mul;
            At->data[jcol].im *= mul;
          }
        }
      }

      ilo = 0;
      ihi = At->size[0];
      if (At->size[0] <= 1) {
        ihi = 1;
      } else {
        do {
          exitg3 = 0;
          i = 0;
          j = 0;
          notdone = false;
          ii = ihi;
          exitg1 = false;
          while ((!exitg1) && (ii > 0)) {
            nzcount = 0;
            i = ii;
            j = ihi;
            jcol = 1;
            exitg4 = false;
            while ((!exitg4) && (jcol <= ihi)) {
              d_At = ((At->data[(ii + At->size[0] * (jcol - 1)) - 1].re != 0.0) ||
                      (At->data[(ii + At->size[0] * (jcol - 1)) - 1].im != 0.0));
              if (d_At || (ii == jcol)) {
                if (nzcount == 0) {
                  j = jcol;
                  nzcount = 1;
                  jcol++;
                } else {
                  nzcount = 2;
                  exitg4 = true;
                }
              } else {
                jcol++;
              }
            }

            if (nzcount < 2) {
              notdone = true;
              exitg1 = true;
            } else {
              ii--;
            }
          }

          if (!notdone) {
            exitg3 = 2;
          } else {
            nzcount = At->size[0];
            if (i != ihi) {
              for (jcol = 0; jcol + 1 <= nzcount; jcol++) {
                atmp = At->data[(i + At->size[0] * jcol) - 1];
                At->data[(i + At->size[0] * jcol) - 1] = At->data[(ihi +
                  At->size[0] * jcol) - 1];
                At->data[(ihi + At->size[0] * jcol) - 1] = atmp;
              }
            }

            if (j != ihi) {
              for (jcol = 0; jcol + 1 <= ihi; jcol++) {
                atmp = At->data[jcol + At->size[0] * (j - 1)];
                At->data[jcol + At->size[0] * (j - 1)] = At->data[jcol +
                  At->size[0] * (ihi - 1)];
                At->data[jcol + At->size[0] * (ihi - 1)] = atmp;
              }
            }

            ihi--;
            if (ihi == 1) {
              exitg3 = 1;
            }
          }
        } while (exitg3 == 0);

        if (exitg3 == 1) {
        } else {
          do {
            exitg2 = 0;
            i = 0;
            j = 0;
            notdone = false;
            jcol = ilo + 1;
            exitg1 = false;
            while ((!exitg1) && (jcol <= ihi)) {
              nzcount = 0;
              i = ihi;
              j = jcol;
              ii = ilo + 1;
              exitg4 = false;
              while ((!exitg4) && (ii <= ihi)) {
                d_At = ((At->data[(ii + At->size[0] * (jcol - 1)) - 1].re != 0.0)
                        || (At->data[(ii + At->size[0] * (jcol - 1)) - 1].im !=
                            0.0));
                if (d_At || (ii == jcol)) {
                  if (nzcount == 0) {
                    i = ii;
                    nzcount = 1;
                    ii++;
                  } else {
                    nzcount = 2;
                    exitg4 = true;
                  }
                } else {
                  ii++;
                }
              }

              if (nzcount < 2) {
                notdone = true;
                exitg1 = true;
              } else {
                jcol++;
              }
            }

            if (!notdone) {
              exitg2 = 1;
            } else {
              nzcount = At->size[0];
              if (i != ilo + 1) {
                for (jcol = ilo; jcol + 1 <= nzcount; jcol++) {
                  atmp = At->data[(i + At->size[0] * jcol) - 1];
                  At->data[(i + At->size[0] * jcol) - 1] = At->data[ilo +
                    At->size[0] * jcol];
                  At->data[ilo + At->size[0] * jcol] = atmp;
                }
              }

              if (j != ilo + 1) {
                for (jcol = 0; jcol + 1 <= ihi; jcol++) {
                  atmp = At->data[jcol + At->size[0] * (j - 1)];
                  At->data[jcol + At->size[0] * (j - 1)] = At->data[jcol +
                    At->size[0] * ilo];
                  At->data[jcol + At->size[0] * ilo] = atmp;
                }
              }

              ilo++;
              if (ilo + 1 == ihi) {
                exitg2 = 1;
              }
            }
          } while (exitg2 == 0);
        }
      }

      nzcount = At->size[0];
      if ((!(At->size[0] <= 1)) && (!(ihi < ilo + 3))) {
        for (jcol = ilo; jcol + 1 < ihi - 1; jcol++) {
          for (ii = ihi - 1; ii + 1 > jcol + 2; ii--) {
            b_At = At->data[(ii + At->size[0] * jcol) - 1];
            c_At = At->data[ii + At->size[0] * jcol];
            xzlartg(b_At, c_At, &c, &atmp, &At->data[(ii + At->size[0] * jcol) -
                    1]);
            At->data[ii + At->size[0] * jcol].re = 0.0;
            At->data[ii + At->size[0] * jcol].im = 0.0;
            for (j = jcol + 1; j + 1 <= nzcount; j++) {
              absxk = atmp.re * At->data[ii + At->size[0] * j].re - atmp.im *
                At->data[ii + At->size[0] * j].im;
              ctoc = atmp.re * At->data[ii + At->size[0] * j].im + atmp.im *
                At->data[ii + At->size[0] * j].re;
              stemp_re = c * At->data[(ii + At->size[0] * j) - 1].re + absxk;
              absxk = c * At->data[(ii + At->size[0] * j) - 1].im + ctoc;
              ctoc = At->data[(ii + At->size[0] * j) - 1].re;
              cfrom1 = At->data[(ii + At->size[0] * j) - 1].im;
              cto1 = At->data[(ii + At->size[0] * j) - 1].im;
              mul = At->data[(ii + At->size[0] * j) - 1].re;
              At->data[ii + At->size[0] * j].re = c * At->data[ii + At->size[0] *
                j].re - (atmp.re * ctoc + atmp.im * cfrom1);
              At->data[ii + At->size[0] * j].im = c * At->data[ii + At->size[0] *
                j].im - (atmp.re * cto1 - atmp.im * mul);
              At->data[(ii + At->size[0] * j) - 1].re = stemp_re;
              At->data[(ii + At->size[0] * j) - 1].im = absxk;
            }

            atmp.re = -atmp.re;
            atmp.im = -atmp.im;
            for (i = 0; i + 1 <= ihi; i++) {
              absxk = atmp.re * At->data[i + At->size[0] * (ii - 1)].re -
                atmp.im * At->data[i + At->size[0] * (ii - 1)].im;
              ctoc = atmp.re * At->data[i + At->size[0] * (ii - 1)].im + atmp.im
                * At->data[i + At->size[0] * (ii - 1)].re;
              stemp_re = c * At->data[i + At->size[0] * ii].re + absxk;
              absxk = c * At->data[i + At->size[0] * ii].im + ctoc;
              ctoc = At->data[i + At->size[0] * ii].re;
              cfrom1 = At->data[i + At->size[0] * ii].im;
              cto1 = At->data[i + At->size[0] * ii].im;
              mul = At->data[i + At->size[0] * ii].re;
              At->data[i + At->size[0] * (ii - 1)].re = c * At->data[i +
                At->size[0] * (ii - 1)].re - (atmp.re * ctoc + atmp.im * cfrom1);
              At->data[i + At->size[0] * (ii - 1)].im = c * At->data[i +
                At->size[0] * (ii - 1)].im - (atmp.re * cto1 - atmp.im * mul);
              At->data[i + At->size[0] * ii].re = stemp_re;
              At->data[i + At->size[0] * ii].im = absxk;
            }
          }
        }
      }

      xzhgeqz(At, ilo + 1, ihi, &nzcount, alpha1, beta1);
      if ((nzcount == 0) && ilascl) {
        xzlascl(anrmto, anrm, alpha1);
      }
    }
  }

  emxFree_creal_T(&At);
  *info = nzcount;
}

/*
 * File trailer for xzgeev.c
 *
 * [EOF]
 */
