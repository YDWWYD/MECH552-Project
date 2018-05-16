/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: myzp2ss.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "myzp2ss.h"
#include "butterBandpassOnly_emxutil.h"
#include "relop.h"
#include "butterBandpassOnly_rtwutil.h"

/* Function Definitions */

/*
 * ZP2SS  Zero-pole to state-space conversion.
 *    [A,B,C,D] = ZP2SS(Z,P,K)  calculates a state-space representation:
 *        .
 *        x = Ax + Bu
 *        y = Cx + Du
 *
 *    for a system given a set of pole locations in column vector P,
 *    a matrix Z with the zero locations in as many columns as there are
 *    outputs, and the gains for each numerator transfer function in
 *    vector K.  The A,B,C,D matrices are returned in block diagonal
 *    form.
 *
 *    The poles and zeros must correspond to a proper system. If the poles
 *    or zeros are complex, they must appear in complex conjugate pairs,
 *    i.e., the corresponding transfer function must be real.
 *
 *    See also SS2ZP, ZP2TF, TF2ZP, TF2SS, SS2TF.
 * Arguments    : emxArray_creal_T *p
 *                double k
 *                emxArray_creal_T *a
 *                emxArray_real_T *b
 *                emxArray_real_T *c
 *                double *d
 * Return Type  : void
 */
void myzp2ss(emxArray_creal_T *p, double k, emxArray_creal_T *a, emxArray_real_T
             *b, emxArray_real_T *c, double *d)
{
  int r1;
  int i1;
  int np;
  double b_d;
  unsigned int i;
  emxArray_creal_T *result;
  emxArray_int8_T *reshapes_f2;
  emxArray_real_T *varargin_2;
  emxArray_real_T *b1;
  emxArray_real_T *r2;
  creal_T b_c[3];
  creal_T x[2];
  double wn;
  double den[3];
  int b_k;
  int b_r2;
  double b_b1[2];
  double t[4];
  double B[4];
  double a22;
  double Y[4];
  boolean_T empty_non_axis_sizes;
  int b_result;
  int i2;
  double b_Y[4];

  /*    Thanks to G.F. Franklin */
  /*    Copyright 1984-2008 The MathWorks, Inc. */
  /* [z,p,k,isSIMO] = parse_input(z,p,k); */
  /* isSIMO = 0; */
  r1 = p->size[0];
  i1 = p->size[0];
  p->size[0] = r1;
  emxEnsureCapacity_creal_T1(p, i1);

  /*  if isSIMO */
  /*      % If it's multi-output, we can't use the nice algorithm */
  /*      % that follows, so use the numerically unreliable method */
  /*      % of going through polynomial form, and then return. */
  /*      [num,den] = zp2tf(z,p,k); % Suppress compile-time diagnostics */
  /*      [a,b,c,d] = tf2ss(num,den); */
  /*      return */
  /*  end */
  /*  Strip infinities and throw away. */
  /* p = p(isfinite(p)); */
  /* z = z(isfinite(z)); */
  /*  Group into complex pairs */
  np = p->size[0];

  /*  try */
  /*      % z and p should have real elements and exact complex conjugate pair. */
  /*      z = cplxpair(z,0); */
  /*      p = cplxpair(p,0); */
  /*  catch */
  /*      % If fail, revert to use the old default tolerance. */
  /*      % The use of tolerance in checking for real entries and conjugate pairs */
  /*      % may result in misinterpretation for edge cases. Please review the */
  /*      % process of how z and p are generated. */
  /*      z = cplxpair(z,1e6*nz*norm(z)*eps + eps); */
  /*      p = cplxpair(p,1e6*np*norm(p)*eps + eps); */
  /*  end */
  /*  Initialize state-space matrices for running series */
  i1 = a->size[0] * a->size[1];
  a->size[0] = 0;
  a->size[1] = 0;
  emxEnsureCapacity_creal_T(a, i1);
  i1 = b->size[0];
  b->size[0] = 0;
  emxEnsureCapacity_real_T1(b, i1);
  i1 = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = 0;
  emxEnsureCapacity_real_T(c, i1);
  b_d = 1.0;

  /*  If odd number of poles AND zeros, convert the pole and zero */
  /*  at the end into state-space. */
  /*    H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2)) */
  r1 = (int)rt_remd_snf(p->size[0], 2.0);
  if (r1 != 0) {
    i1 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = 1;
    emxEnsureCapacity_creal_T(a, i1);
    a->data[0] = p->data[p->size[0] - 1];
    i1 = b->size[0];
    b->size[0] = 1;
    emxEnsureCapacity_real_T1(b, i1);
    b->data[0] = 1.0;
    i1 = c->size[0] * c->size[1];
    c->size[0] = 1;
    c->size[1] = 1;
    emxEnsureCapacity_real_T(c, i1);
    c->data[0] = 1.0;
    b_d = 0.0;
    np = p->size[0] - 1;

    /* end  */
    /*  If odd number of zeros only, convert the zero at the */
    /*  end, along with a pole-pair into state-space. */
    /*    H(s) = (s+num(2))/(s^2+den(2)s+den(3))  */
  }

  /*  Now we have an even number of poles and zeros, although not  */
  /*  necessarily the same number - there may be more poles. */
  /*    H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3)) */
  /*  Loop through rest of pairs, connecting in series to build the model. */
  i = 1U;

  /*  Take care of any left over unmatched pole pairs. */
  /*    H(s) = 1/(s^2+den(2)s+den(3)) */
  emxInit_creal_T(&result, 2);
  emxInit_int8_T1(&reshapes_f2, 2);
  emxInit_real_T1(&varargin_2, 2);
  emxInit_real_T1(&b1, 2);
  emxInit_real_T1(&r2, 2);
  while ((double)i < np) {
    for (i1 = 0; i1 < 2; i1++) {
      x[i1] = p->data[(i1 + (int)i) - 1];
    }

    b_c[0].re = 1.0;
    b_c[0].im = 0.0;
    for (r1 = 0; r1 < 2; r1++) {
      wn = b_c[r1].re;
      b_c[r1 + 1].re = -x[r1].re * b_c[r1].re - -x[r1].im * b_c[r1].im;
      b_c[r1 + 1].im = -x[r1].re * b_c[r1].im + -x[r1].im * wn;
      b_k = r1;
      while (b_k + 1 > 1) {
        b_c[1].re -= x[r1].re * b_c[0].re - x[r1].im * b_c[0].im;
        b_c[1].im -= x[r1].re * b_c[0].im + x[r1].im * b_c[0].re;
        b_k = 0;
      }
    }

    for (i1 = 0; i1 < 3; i1++) {
      den[i1] = b_c[i1].re;
    }

    for (i1 = 0; i1 < 2; i1++) {
      x[i1] = p->data[(i1 + (int)i) - 1];
    }

    for (b_k = 0; b_k < 2; b_k++) {
      b_b1[b_k] = rt_hypotd_snf(x[b_k].re, x[b_k].im);
    }

    wn = sqrt(b_b1[0] * b_b1[1]);
    if (wn == 0.0) {
      wn = 1.0;
    }

    b_b1[0] = 1.0;
    b_b1[1] = 1.0 / wn;
    for (i1 = 0; i1 < 4; i1++) {
      t[i1] = 0.0;
    }

    /*  Balancing transformation */
    B[0] = -den[1];
    B[2] = -den[2];
    for (r1 = 0; r1 < 2; r1++) {
      t[r1 + (r1 << 1)] = b_b1[r1];
      B[1 + (r1 << 1)] = 1.0 - (double)r1;
    }

    if (t[1] > t[0]) {
      r1 = 1;
      b_r2 = 0;
    } else {
      r1 = 0;
      b_r2 = 1;
    }

    wn = t[b_r2] / t[r1];
    a22 = t[2 + b_r2] - wn * t[2 + r1];
    for (b_k = 0; b_k < 2; b_k++) {
      Y[1 + (b_k << 1)] = (B[b_r2 + (b_k << 1)] - B[r1 + (b_k << 1)] * wn) / a22;
      Y[b_k << 1] = (B[r1 + (b_k << 1)] - Y[1 + (b_k << 1)] * t[2 + r1]) / t[r1];
    }

    if (t[1] > t[0]) {
      r1 = 1;
      b_r2 = 0;
    } else {
      r1 = 0;
      b_r2 = 1;
    }

    wn = t[b_r2] / t[r1];
    b_b1[1] = ((1.0 - (double)b_r2) - (1.0 - (double)r1) * wn) / (t[2 + b_r2] -
      wn * t[2 + r1]);
    b_b1[0] = ((1.0 - (double)r1) - b_b1[1] * t[2 + r1]) / t[r1];

    /*  [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1); */
    /*  Next lines perform series connection  */
    if (!((a->size[0] == 0) || (a->size[1] == 0))) {
      b_k = a->size[0];
    } else {
      r1 = a->size[0];
      if (!(r1 == 0)) {
        b_k = a->size[0];
      } else {
        b_k = a->size[0];
        if (!(b_k > 0)) {
          b_k = 0;
        }

        r1 = a->size[0];
        if (r1 > b_k) {
          b_k = a->size[0];
        }
      }
    }

    empty_non_axis_sizes = (b_k == 0);
    if (empty_non_axis_sizes || (!((a->size[0] == 0) || (a->size[1] == 0)))) {
      b_result = a->size[1];
    } else {
      b_result = 0;
    }

    if (empty_non_axis_sizes) {
      r1 = 2;
    } else {
      r1 = a->size[0];
      if (!(r1 == 0)) {
        r1 = 2;
      } else {
        r1 = 0;
      }
    }

    i1 = reshapes_f2->size[0] * reshapes_f2->size[1];
    reshapes_f2->size[0] = b_k;
    reshapes_f2->size[1] = r1;
    emxEnsureCapacity_int8_T1(reshapes_f2, i1);
    b_r2 = b_k * r1;
    for (i1 = 0; i1 < b_r2; i1++) {
      reshapes_f2->data[i1] = 0;
    }

    i1 = result->size[0] * result->size[1];
    result->size[0] = b_k;
    result->size[1] = b_result + reshapes_f2->size[1];
    emxEnsureCapacity_creal_T(result, i1);
    for (i1 = 0; i1 < b_result; i1++) {
      for (i2 = 0; i2 < b_k; i2++) {
        result->data[i2 + result->size[0] * i1] = a->data[i2 + b_k * i1];
      }
    }

    b_r2 = reshapes_f2->size[1];
    for (i1 = 0; i1 < b_r2; i1++) {
      r1 = reshapes_f2->size[0];
      for (i2 = 0; i2 < r1; i2++) {
        result->data[i2 + result->size[0] * (i1 + b_result)].re =
          reshapes_f2->data[i2 + reshapes_f2->size[0] * i1];
        result->data[i2 + result->size[0] * (i1 + b_result)].im = 0.0;
      }
    }

    i1 = b1->size[0] * b1->size[1];
    b1->size[0] = 2;
    b1->size[1] = c->size[1];
    emxEnsureCapacity_real_T(b1, i1);
    for (i1 = 0; i1 < 2; i1++) {
      b_r2 = c->size[1];
      for (i2 = 0; i2 < b_r2; i2++) {
        b1->data[i1 + b1->size[0] * i2] = b_b1[i1] * c->data[c->size[0] * i2];
      }
    }

    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < 2; i2++) {
        b_Y[i1 + (i2 << 1)] = 0.0;
        for (r1 = 0; r1 < 2; r1++) {
          b_Y[i1 + (i2 << 1)] += Y[i1 + (r1 << 1)] * t[r1 + (i2 << 1)];
        }
      }
    }

    i1 = varargin_2->size[0] * varargin_2->size[1];
    varargin_2->size[0] = 2;
    varargin_2->size[1] = b1->size[1] + 2;
    emxEnsureCapacity_real_T(varargin_2, i1);
    b_r2 = b1->size[1];
    for (i1 = 0; i1 < b_r2; i1++) {
      for (i2 = 0; i2 < 2; i2++) {
        varargin_2->data[i2 + varargin_2->size[0] * i1] = b1->data[i2 + b1->
          size[0] * i1];
      }
    }

    for (i1 = 0; i1 < 2; i1++) {
      for (i2 = 0; i2 < 2; i2++) {
        varargin_2->data[i2 + varargin_2->size[0] * (i1 + b1->size[1])] = b_Y[i2
          + (i1 << 1)];
      }
    }

    if (!((result->size[0] == 0) || (result->size[1] == 0))) {
      b_result = result->size[1];
    } else {
      b_result = varargin_2->size[1];
    }

    if (!((result->size[0] == 0) || (result->size[1] == 0))) {
      r1 = result->size[0];
    } else {
      r1 = 0;
    }

    i1 = a->size[0] * a->size[1];
    a->size[0] = r1 + 2;
    a->size[1] = b_result;
    emxEnsureCapacity_creal_T(a, i1);
    for (i1 = 0; i1 < b_result; i1++) {
      for (i2 = 0; i2 < r1; i2++) {
        a->data[i2 + a->size[0] * i1] = result->data[i2 + r1 * i1];
      }
    }

    for (i1 = 0; i1 < b_result; i1++) {
      for (i2 = 0; i2 < 2; i2++) {
        a->data[(i2 + r1) + a->size[0] * i1].re = varargin_2->data[i2 + (i1 << 1)];
        a->data[(i2 + r1) + a->size[0] * i1].im = 0.0;
      }
    }

    r1 = b->size[0];
    i1 = b->size[0];
    b->size[0] = r1 + 2;
    emxEnsureCapacity_real_T1(b, i1);
    for (i1 = 0; i1 < 2; i1++) {
      b->data[r1 + i1] = b_b1[i1] * b_d;
    }

    for (i1 = 0; i1 < 2; i1++) {
      b_b1[i1] = 0.0;
      for (i2 = 0; i2 < 2; i2++) {
        b_b1[i1] += (double)i2 * t[i2 + (i1 << 1)];
      }
    }

    i1 = r2->size[0] * r2->size[1];
    r2->size[0] = 1;
    r2->size[1] = c->size[1] + 2;
    emxEnsureCapacity_real_T(r2, i1);
    b_r2 = c->size[1];
    for (i1 = 0; i1 < b_r2; i1++) {
      r2->data[r2->size[0] * i1] = 0.0 * c->data[c->size[0] * i1];
    }

    for (i1 = 0; i1 < 2; i1++) {
      r2->data[r2->size[0] * (i1 + c->size[1])] = b_b1[i1];
    }

    i1 = c->size[0] * c->size[1];
    c->size[0] = 1;
    c->size[1] = r2->size[1];
    emxEnsureCapacity_real_T(c, i1);
    b_r2 = r2->size[1];
    for (i1 = 0; i1 < b_r2; i1++) {
      c->data[c->size[0] * i1] = r2->data[r2->size[0] * i1];
    }

    b_d = 0.0;
    i += 2U;
  }

  emxFree_real_T(&r2);
  emxFree_real_T(&b1);
  emxFree_real_T(&varargin_2);
  emxFree_int8_T(&reshapes_f2);
  emxFree_creal_T(&result);

  /*  Apply gain k: */
  i1 = c->size[0] * c->size[1];
  c->size[0] = 1;
  emxEnsureCapacity_real_T(c, i1);
  r1 = c->size[0];
  b_r2 = c->size[1];
  b_r2 *= r1;
  for (i1 = 0; i1 < b_r2; i1++) {
    c->data[i1] *= k;
  }

  b_d *= k;

  /* ---------------------------------------------------------------------------- */
  /*  function [z,p,k,isSIMO] = parse_input(z,p,k) */
  /*  %PARSE_INPUT   Make sure input args are valid. */
  /*   */
  /*  % Initially assume it is a SISO system */
  /*  isSIMO = 0; */
  /*   */
  /*  % Check that p is a vector */
  /*  if ~any(size(p)<2), */
  /*     ctrlMsgUtils.error('Controllib:general:pNotVector') */
  /*  end */
  /*  % Columnize p */
  /*  p = p(:); */
  /*   */
  /*  % Check that k is a vector */
  /*  if ~any(size(k)<2), */
  /*     ctrlMsgUtils.error('Controllib:general:kNotVector') */
  /*  end */
  /*  % Columnize k */
  /*  k = k(:); */
  /*   */
  /*  % Check size of z */
  /*  if any(size(z)<2), */
  /*     % z is a vector or an empty, columnize it */
  /*     z = z(:); */
  /*  else */
  /*     % z is a matrix */
  /*     isSIMO = 1; */
  /*  end */
  /*   */
  /*  % Check for properness */
  /*  if size(z,1) > length(p), */
  /*     % improper */
  /*     ctrlMsgUtils.error('Controllib:general:improperSystem') */
  /*  end */
  /*   */
  /*  % Check for the appropriate length of k */
  /*  if length(k) ~= size(z,2) && (~isempty(z)) */
  /*     ctrlMsgUtils.error('Controllib:general:zkLengthMismatch') */
  /*  end */
  *d = b_d;
}

/*
 * File trailer for myzp2ss.c
 *
 * [EOF]
 */
