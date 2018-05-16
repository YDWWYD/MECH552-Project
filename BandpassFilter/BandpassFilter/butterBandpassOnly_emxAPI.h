/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: butterBandpassOnly_emxAPI.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef BUTTERBANDPASSONLY_EMXAPI_H
#define BUTTERBANDPASSONLY_EMXAPI_H

/* Include Files */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "butterBandpassOnly_types.h"

/* Function Declarations */
extern emxArray_creal_T *emxCreateND_creal_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_creal_T *emxCreateWrapperND_creal_T(creal_T *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_creal_T *emxCreateWrapper_creal_T(creal_T *data, int rows, int
  cols);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_creal_T *emxCreate_creal_T(int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_creal_T(emxArray_creal_T *emxArray);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxInitArray_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for butterBandpassOnly_emxAPI.h
 *
 * [EOF]
 */
