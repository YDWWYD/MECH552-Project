/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: butterBandpassOnly_emxutil.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

#ifndef BUTTERBANDPASSONLY_EMXUTIL_H
#define BUTTERBANDPASSONLY_EMXUTIL_H

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
extern void emxEnsureCapacity_creal_T(emxArray_creal_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_creal_T1(emxArray_creal_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_int32_T1(emxArray_int32_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_int8_T(emxArray_int8_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_int8_T1(emxArray_int8_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T1(emxArray_real_T *emxArray, int oldNumel);
extern void emxFreeMatrix_cell_wrap_2(cell_wrap_2 pMatrix[2]);
extern void emxFree_creal_T(emxArray_creal_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_int8_T(emxArray_int8_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInitMatrix_cell_wrap_2(cell_wrap_2 pMatrix[2]);
extern void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);
extern void emxInit_creal_T1(emxArray_creal_T **pEmxArray, int numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
extern void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions);
extern void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions);
extern void emxInit_int8_T1(emxArray_int8_T **pEmxArray, int numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for butterBandpassOnly_emxutil.h
 *
 * [EOF]
 */
