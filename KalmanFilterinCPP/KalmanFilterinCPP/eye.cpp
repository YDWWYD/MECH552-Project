/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: eye.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "eye.h"
#include "butterBandpassOnly_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const double varargin_1[2]
 *                emxArray_real_T *I
 * Return Type  : void
 */
void b_eye(const double varargin_1[2], emxArray_real_T *I)
{
  int u0;
  int d;
  int loop_ub;
  u0 = (int)varargin_1[0];
  d = (int)varargin_1[1];
  if (u0 < d) {
    d = u0;
  }

  u0 = I->size[0] * I->size[1];
  I->size[0] = (int)varargin_1[0];
  I->size[1] = (int)varargin_1[1];
  emxEnsureCapacity_real_T(I, u0);
  loop_ub = (int)varargin_1[0] * (int)varargin_1[1];
  for (u0 = 0; u0 < loop_ub; u0++) {
    I->data[u0] = 0.0;
  }

  if (d > 0) {
    for (u0 = 0; u0 + 1 <= d; u0++) {
      I->data[u0 + I->size[0] * u0] = 1.0;
    }
  }
}

/*
 * Arguments    : double varargin_1
 *                emxArray_real_T *I
 * Return Type  : void
 */
void eye(double varargin_1, emxArray_real_T *I)
{
  int k;
  int loop_ub;
  k = I->size[0] * I->size[1];
  I->size[0] = (int)varargin_1;
  I->size[1] = (int)varargin_1;
  emxEnsureCapacity_real_T(I, k);
  loop_ub = (int)varargin_1 * (int)varargin_1;
  for (k = 0; k < loop_ub; k++) {
    I->data[k] = 0.0;
  }

  if ((int)varargin_1 > 0) {
    for (k = 0; k + 1 <= (int)varargin_1; k++) {
      I->data[k + I->size[0] * k] = 1.0;
    }
  }
}

/*
 * File trailer for eye.c
 *
 * [EOF]
 */
