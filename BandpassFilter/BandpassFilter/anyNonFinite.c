/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: anyNonFinite.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "anyNonFinite.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *x
 * Return Type  : boolean_T
 */
boolean_T anyNonFinite(const emxArray_creal_T *x)
{
  boolean_T p;
  int nx;
  int k;
  nx = x->size[0] * x->size[1];
  p = true;
  for (k = 0; k + 1 <= nx; k++) {
    if (p && ((!(rtIsInf(x->data[k].re) || rtIsInf(x->data[k].im))) &&
              (!(rtIsNaN(x->data[k].re) || rtIsNaN(x->data[k].im))))) {
      p = true;
    } else {
      p = false;
    }
  }

  return !p;
}

/*
 * File trailer for anyNonFinite.c
 *
 * [EOF]
 */
