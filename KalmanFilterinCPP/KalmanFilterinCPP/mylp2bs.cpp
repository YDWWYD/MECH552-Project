/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mylp2bs.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "mylp2bs.h"
#include "butterBandpassOnly_emxutil.h"
#include "mrdivide.h"
#include "xtrsm.h"
#include "xgetrf.h"
#include "qrsolve.h"
#include "eye.h"
#include "inv.h"

/* Function Definitions */

/*
 * LP2BS Lowpass to bandstop analog filter transformation.
 *    [NUMT,DENT] = LP2BS(NUM,DEN,Wo,Bw) transforms the lowpass filter
 *    prototype NUM(s)/DEN(s) with unity cutoff frequency to a
 *    bandstop filter with center frequency Wo and bandwidth Bw.
 *    [AT,BT,CT,DT] = LP2BS(A,B,C,D,Wo,Bw) does the same when the
 *    filter is described in state-space form.
 *
 *    % Example:
 *    %   Design a Elliptic analog lowpass filter with 5dB of ripple in the
 *    %   passband and a stopband 90 decibels down. Transform this filter to
 *    %   a bandstop filter with center frequency 24Hz and Bandwidth of 10Hz.
 *
 *    Wo= 24; Bw=10                   % Define Center frequency and Bandwidth
 *    [z,p,k]=ellipap(6,5,90);        % Lowpass filter prototype
 *    [b,a]=zp2tf(z,p,k);             % Specify filter in polynomial form
 *    [num,den]=lp2bs(b,a,Wo,Bw);     % Convert LPF to BSF
 *    freqs(num,den)                  % Frequency response of analog filter
 *
 *    See also BILINEAR, IMPINVAR, LP2BP, LP2LP and LP2HP
 * Arguments    : const emxArray_creal_T *a
 *                const emxArray_real_T *b
 *                const emxArray_real_T *c
 *                double d
 *                double wo
 *                const double bw_data[]
 *                emxArray_creal_T *at
 *                emxArray_creal_T *bt
 *                emxArray_creal_T *ct
 *                creal_T *dt
 * Return Type  : void
 */
void mylp2bs(const emxArray_creal_T *a, const emxArray_real_T *b, const
             emxArray_real_T *c, double d, double wo, const double bw_data[],
             emxArray_creal_T *at, emxArray_creal_T *bt, emxArray_creal_T *ct,
             creal_T *dt)
{
  emxArray_creal_T *b_b;
  double temp_re;
  double y;
  int i3;
  int c_b;
  int d_b;
  emxArray_real_T *e_b;
  emxArray_creal_T *result;
  boolean_T empty_non_axis_sizes;
  int b_result;
  int i4;
  emxArray_real_T *c_result;
  int c_idx_0;
  cell_wrap_2 reshapes[2];
  emxArray_creal_T *Y;
  emxArray_int32_T *ipiv;
  unsigned int unnamed_idx_0;
  emxArray_creal_T *b_a;
  double temp_im;
  double a_re;
  emxInit_creal_T(&b_b, 2);

  /*    Author(s): J.N. Little and G.F. Franklin, 8-4-87 */
  /*    Copyright 1988-2013 The MathWorks, Inc. */
  /*  if nargin == 4		% Transfer function case */
  /*      % handle column vector inputs: convert to rows */
  /*      if size(a,2) == 1 */
  /*          a = a(:).'; */
  /*      end */
  /*      if size(b,2) == 1 */
  /*          b = b(:).'; */
  /*      end */
  /*       */
  /*      % Cast to enforce precision rules */
  /*      wo = signal.internal.sigcasttofloat(c,'double','lp2bs','Wo',... */
  /*        'allownumeric'); */
  /*      bw = signal.internal.sigcasttofloat(d,'double','lp2bs','Bw',... */
  /*        'allownumeric'); */
  /*      if any([signal.internal.sigcheckfloattype(b,'single','lp2bs','DEN')... */
  /*          signal.internal.sigcheckfloattype(a,'single','lp2bs','NUM')]) */
  /*        b = single(b); */
  /*        a = single(a); */
  /*      end */
  /*      % Transform to state-space */
  /*      [a,b,c,d] = tf2ss(a,b); */
  /*   */
  /*  % Cast to enforce precision rules */
  /*  elseif (nargin == 6) */
  /*    if any([signal.internal.sigcheckfloattype(a,'single','lp2bs','A')... */
  /*        signal.internal.sigcheckfloattype(b,'single','lp2bs','B')... */
  /*        signal.internal.sigcheckfloattype(c,'single','lp2bs','C')... */
  /*        signal.internal.sigcheckfloattype(d,'single','lp2bs','D')]) */
  /*      a = single(a); */
  /*      b = single(b);     */
  /*      c = single(c); */
  /*      d = single(d); */
  /*    end */
  /*    wo = signal.internal.sigcasttofloat(wo,'double','lp2bs','Wo',... */
  /*      'allownumeric'); */
  /*    bw = signal.internal.sigcasttofloat(bw,'double','lp2bs','Bw',... */
  /*      'allownumeric'); */
  /*  else */
  /*    error(message('signal:lp2bs:MustHaveNInputs')); */
  /*  end */
  /*   */
  /*  error(abcdchk(a,b,c,d)); */
  /*  Transform lowpass to bandstop */
  temp_re = 1.0 / bw_data[0] * wo;
  y = 1.0 / temp_re * wo;
  inv(a, b_b);
  i3 = b_b->size[0] * b_b->size[1];
  emxEnsureCapacity_creal_T(b_b, i3);
  c_b = b_b->size[0];
  d_b = b_b->size[1];
  d_b *= c_b;
  for (i3 = 0; i3 < d_b; i3++) {
    b_b->data[i3].re *= y;
    b_b->data[i3].im *= y;
  }

  emxInit_real_T1(&e_b, 2);
  eye(c->size[1], e_b);
  i3 = e_b->size[0] * e_b->size[1];
  emxEnsureCapacity_real_T(e_b, i3);
  c_b = e_b->size[0];
  d_b = e_b->size[1];
  d_b *= c_b;
  for (i3 = 0; i3 < d_b; i3++) {
    e_b->data[i3] *= wo;
  }

  emxInit_creal_T(&result, 2);
  if (!((b_b->size[0] == 0) || (b_b->size[1] == 0))) {
    c_b = b_b->size[0];
  } else if (!((e_b->size[0] == 0) || (e_b->size[1] == 0))) {
    c_b = e_b->size[0];
  } else {
    c_b = b_b->size[0];
    if (!(c_b > 0)) {
      c_b = 0;
    }

    if (e_b->size[0] > c_b) {
      c_b = e_b->size[0];
    }
  }

  empty_non_axis_sizes = (c_b == 0);
  if (empty_non_axis_sizes || (!((b_b->size[0] == 0) || (b_b->size[1] == 0)))) {
    b_result = b_b->size[1];
  } else {
    b_result = 0;
  }

  if (empty_non_axis_sizes || (!((e_b->size[0] == 0) || (e_b->size[1] == 0)))) {
    d_b = e_b->size[1];
  } else {
    d_b = 0;
  }

  i3 = result->size[0] * result->size[1];
  result->size[0] = c_b;
  result->size[1] = b_result + d_b;
  emxEnsureCapacity_creal_T(result, i3);
  for (i3 = 0; i3 < b_result; i3++) {
    for (i4 = 0; i4 < c_b; i4++) {
      result->data[i4 + result->size[0] * i3] = b_b->data[i4 + c_b * i3];
    }
  }

  for (i3 = 0; i3 < d_b; i3++) {
    for (i4 = 0; i4 < c_b; i4++) {
      result->data[i4 + result->size[0] * (i3 + b_result)].re = e_b->data[i4 +
        c_b * i3];
      result->data[i4 + result->size[0] * (i3 + b_result)].im = 0.0;
    }
  }

  eye(c->size[1], e_b);
  i3 = e_b->size[0] * e_b->size[1];
  emxEnsureCapacity_real_T(e_b, i3);
  c_b = e_b->size[0];
  d_b = e_b->size[1];
  d_b *= c_b;
  for (i3 = 0; i3 < d_b; i3++) {
    e_b->data[i3] *= -wo;
  }

  emxInit_real_T1(&c_result, 2);
  if (!((e_b->size[0] == 0) || (e_b->size[1] == 0))) {
    c_idx_0 = e_b->size[0];
  } else {
    c_b = c->size[1];
    d_b = c->size[1];
    if (!((c_b == 0) || (d_b == 0))) {
      c_idx_0 = c->size[1];
    } else {
      c_idx_0 = e_b->size[0];
      if (!(c_idx_0 > 0)) {
        c_idx_0 = 0;
      }

      c_b = c->size[1];
      if (c_b > c_idx_0) {
        c_idx_0 = c->size[1];
      }
    }
  }

  empty_non_axis_sizes = (c_idx_0 == 0);
  if (empty_non_axis_sizes || (!((e_b->size[0] == 0) || (e_b->size[1] == 0)))) {
    b_result = e_b->size[1];
  } else {
    b_result = 0;
  }

  if (empty_non_axis_sizes) {
    d_b = c->size[1];
  } else {
    c_b = c->size[1];
    d_b = c->size[1];
    if (!((c_b == 0) || (d_b == 0))) {
      d_b = c->size[1];
    } else {
      d_b = 0;
    }
  }

  emxInitMatrix_cell_wrap_2(reshapes);
  i3 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
  reshapes[1].f1->size[0] = c_idx_0;
  reshapes[1].f1->size[1] = d_b;
  emxEnsureCapacity_real_T(reshapes[1].f1, i3);
  d_b *= c_idx_0;
  for (i3 = 0; i3 < d_b; i3++) {
    reshapes[1].f1->data[i3] = 0.0;
  }

  i3 = c_result->size[0] * c_result->size[1];
  c_result->size[0] = c_idx_0;
  c_result->size[1] = b_result + reshapes[1].f1->size[1];
  emxEnsureCapacity_real_T(c_result, i3);
  for (i3 = 0; i3 < b_result; i3++) {
    for (i4 = 0; i4 < c_idx_0; i4++) {
      c_result->data[i4 + c_result->size[0] * i3] = e_b->data[i4 + c_idx_0 * i3];
    }
  }

  emxFree_real_T(&e_b);
  d_b = reshapes[1].f1->size[1];
  for (i3 = 0; i3 < d_b; i3++) {
    c_b = reshapes[1].f1->size[0];
    for (i4 = 0; i4 < c_b; i4++) {
      c_result->data[i4 + c_result->size[0] * (i3 + b_result)] = reshapes[1].
        f1->data[i4 + reshapes[1].f1->size[0] * i3];
    }
  }

  emxFreeMatrix_cell_wrap_2(reshapes);
  if (!((result->size[0] == 0) || (result->size[1] == 0))) {
    c_b = result->size[1];
  } else if (!((c_result->size[0] == 0) || (c_result->size[1] == 0))) {
    c_b = c_result->size[1];
  } else {
    c_b = result->size[1];
    if (!(c_b > 0)) {
      c_b = 0;
    }

    if (c_result->size[1] > c_b) {
      c_b = c_result->size[1];
    }
  }

  empty_non_axis_sizes = (c_b == 0);
  if (empty_non_axis_sizes || (!((result->size[0] == 0) || (result->size[1] == 0))))
  {
    b_result = result->size[0];
  } else {
    b_result = 0;
  }

  if (empty_non_axis_sizes || (!((c_result->size[0] == 0) || (c_result->size[1] ==
         0)))) {
    d_b = c_result->size[0];
  } else {
    d_b = 0;
  }

  i3 = at->size[0] * at->size[1];
  at->size[0] = b_result + d_b;
  at->size[1] = c_b;
  emxEnsureCapacity_creal_T(at, i3);
  for (i3 = 0; i3 < c_b; i3++) {
    for (i4 = 0; i4 < b_result; i4++) {
      at->data[i4 + at->size[0] * i3] = result->data[i4 + b_result * i3];
    }
  }

  emxFree_creal_T(&result);
  for (i3 = 0; i3 < c_b; i3++) {
    for (i4 = 0; i4 < d_b; i4++) {
      at->data[(i4 + b_result) + at->size[0] * i3].re = c_result->data[i4 + d_b *
        i3];
      at->data[(i4 + b_result) + at->size[0] * i3].im = 0.0;
    }
  }

  emxFree_real_T(&c_result);
  emxInit_creal_T1(&Y, 1);
  y = 1.0 / temp_re * wo;
  emxInit_int32_T1(&ipiv, 2);
  if ((a->size[0] == 0) || (a->size[1] == 0) || (b->size[0] == 0)) {
    unnamed_idx_0 = (unsigned int)a->size[1];
    i3 = Y->size[0];
    Y->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity_creal_T1(Y, i3);
    d_b = (int)unnamed_idx_0;
    for (i3 = 0; i3 < d_b; i3++) {
      Y->data[i3].re = 0.0;
      Y->data[i3].im = 0.0;
    }
  } else if (a->size[0] == a->size[1]) {
    i3 = b_b->size[0] * b_b->size[1];
    b_b->size[0] = a->size[0];
    b_b->size[1] = a->size[1];
    emxEnsureCapacity_creal_T(b_b, i3);
    d_b = a->size[0] * a->size[1];
    for (i3 = 0; i3 < d_b; i3++) {
      b_b->data[i3] = a->data[i3];
    }

    xgetrf(a->size[1], a->size[1], b_b, a->size[1], ipiv, &c_b);
    i3 = Y->size[0];
    Y->size[0] = b->size[0];
    emxEnsureCapacity_creal_T1(Y, i3);
    d_b = b->size[0];
    for (i3 = 0; i3 < d_b; i3++) {
      Y->data[i3].re = b->data[i3];
      Y->data[i3].im = 0.0;
    }

    for (c_b = 0; c_b + 1 < a->size[1]; c_b++) {
      if (ipiv->data[c_b] != c_b + 1) {
        temp_re = Y->data[c_b].re;
        temp_im = Y->data[c_b].im;
        Y->data[c_b] = Y->data[ipiv->data[c_b] - 1];
        Y->data[ipiv->data[c_b] - 1].re = temp_re;
        Y->data[ipiv->data[c_b] - 1].im = temp_im;
      }
    }

    b_xtrsm(a->size[1], b_b, a->size[1], Y);
    c_xtrsm(a->size[1], b_b, a->size[1], Y);
  } else {
    qrsolve(a, b, Y);
  }

  emxFree_int32_T(&ipiv);
  emxFree_creal_T(&b_b);
  i3 = bt->size[0];
  bt->size[0] = Y->size[0] + c->size[1];
  emxEnsureCapacity_creal_T1(bt, i3);
  d_b = Y->size[0];
  for (i3 = 0; i3 < d_b; i3++) {
    bt->data[i3].re = -(y * Y->data[i3].re);
    bt->data[i3].im = -(y * Y->data[i3].im);
  }

  d_b = c->size[1];
  for (i3 = 0; i3 < d_b; i3++) {
    bt->data[i3 + Y->size[0]].re = -0.0;
    bt->data[i3 + Y->size[0]].im = -0.0;
  }

  emxFree_creal_T(&Y);
  emxInit_creal_T(&b_a, 2);
  mrdivide(c, a, b_a);
  c_b = c->size[1];
  i3 = ct->size[0] * ct->size[1];
  ct->size[0] = 1;
  ct->size[1] = b_a->size[1] + c_b;
  emxEnsureCapacity_creal_T(ct, i3);
  d_b = b_a->size[1];
  for (i3 = 0; i3 < d_b; i3++) {
    ct->data[ct->size[0] * i3] = b_a->data[b_a->size[0] * i3];
  }

  for (i3 = 0; i3 < c_b; i3++) {
    ct->data[ct->size[0] * (i3 + b_a->size[1])].re = 0.0;
    ct->data[ct->size[0] * (i3 + b_a->size[1])].im = 0.0;
  }

  mrdivide(c, a, b_a);
  temp_im = 0.0;
  y = 0.0;
  for (i3 = 0; i3 < b_a->size[1]; i3++) {
    temp_re = b->data[i3];
    a_re = b_a->data[b_a->size[0] * i3].re * temp_re - b_a->data[b_a->size[0] *
      i3].im * 0.0;
    temp_re = b_a->data[b_a->size[0] * i3].re * 0.0 + b_a->data[b_a->size[0] *
      i3].im * temp_re;
    temp_im += a_re;
    y += temp_re;
  }

  emxFree_creal_T(&b_a);
  dt->re = d - temp_im;
  dt->im = 0.0 - y;

  /*  if nargin == 4		% Transfer function case */
  /*  Transform back to transfer function */
  /*      zinf = ltipack.getTolerance('infzero',true); */
  /*      [z,k] = ltipack.sszero(at,bt,ct,dt,[],zinf); */
  /*      num = k * poly(z); */
  /*      den = poly(at); */
  /*      at = num; */
  /*      bt = den; */
  /*  end */
}

/*
 * File trailer for mylp2bs.c
 *
 * [EOF]
 */
