/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: sort1.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 11-May-2018 15:51:52
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "buttordBandpassOnly.h"
#include "sortButtord.h"

/* Function Definitions */

/*
 * Arguments    : double x[2]
 * Return Type  : void
 */
void sortButtord(double x[2])
{
  boolean_T p;
  double tmp;
  if ((x[0] <= x[1]) || rtIsNaN(x[1])) {
    p = true;
  } else {
    p = false;
  }

  if (!p) {
    tmp = x[0];
    x[0] = x[1];
    x[1] = tmp;
  }
}

/*
 * File trailer for sort1.c
 *
 * [EOF]
 */
