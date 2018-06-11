/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzlascl.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "xzlascl.h"
#include "butterBandpassOnly_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : double cfrom
 *                double cto
 *                emxArray_creal_T *A
 * Return Type  : void
 */
void xzlascl(double cfrom, double cto, emxArray_creal_T *A)
{
  double cfromc;
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double mul;
  int i18;
  int loop_ub;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }

    i18 = A->size[0];
    emxEnsureCapacity_creal_T1(A, i18);
    loop_ub = A->size[0];
    for (i18 = 0; i18 < loop_ub; i18++) {
      A->data[i18].re *= mul;
      A->data[i18].im *= mul;
    }
  }
}

/*
 * File trailer for xzlascl.c
 *
 * [EOF]
 */
