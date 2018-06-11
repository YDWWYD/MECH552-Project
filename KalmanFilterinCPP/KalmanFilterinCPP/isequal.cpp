/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: isequal.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "isequal.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_creal_T *varargin_1
 *                const emxArray_creal_T *varargin_2
 * Return Type  : boolean_T
 */
boolean_T isequal(const emxArray_creal_T *varargin_1, const emxArray_creal_T
                  *varargin_2)
{
  boolean_T p;
  boolean_T b_p;
  int k;
  boolean_T exitg1;
  boolean_T b_varargin_1;
  p = false;
  b_p = false;
  if (varargin_1->size[0] == varargin_2->size[0]) {
    b_p = true;
  }

  if (b_p && (!(varargin_1->size[0] == 0)) && (!(varargin_2->size[0] == 0))) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= varargin_2->size[0] - 1)) {
      b_varargin_1 = ((varargin_1->data[k].re == varargin_2->data[k].re) &&
                      (varargin_1->data[k].im == varargin_2->data[k].im));
      if (!b_varargin_1) {
        b_p = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (b_p) {
    p = true;
  }

  return p;
}

/*
 * File trailer for isequal.c
 *
 * [EOF]
 */
