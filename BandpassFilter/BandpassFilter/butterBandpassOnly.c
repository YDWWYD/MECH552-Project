/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: butterBandpassOnly.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "butterBandpassOnly_emxutil.h"
#include "exp.h"
#include "isequal.h"
#include "sortButter.h"
#include "relop.h"
#include "mypoly.h"
#include "mybilinear.h"
#include "mylp2bs.h"
#include "myzp2ss.h"
#include "butterBandpassOnly_rtwutil.h"

/* Function Definitions */

//#include "stdio.h"
//
//int main(void)
//{
//	double Wn[] = { 0.3,0.6 };
//	emxArray_real_T *num;
//	emxArray_creal_T *den;
//	emxInit_real_T1(&num, 2);
//	emxInit_creal_T1(&den, 2);
//	butterBandpassOnly(2, Wn, num, den);
//
//	for (int i = 0; i < num->size[1]; i++)
//	{
//		printf("num[%d] is %.10f\n", i, num->data[i]);
//	}
//
//	for (int i = 0; i < den->size[1]; i++)
//	{
//		printf("Re(den[%d]) is %.8f		Im(den[%d]) is %f\n", i, den->data[i].re, i, den->data[i].im);
//	}
//
//	system("pause");
//	return 0;
//}
/*
 * Cast to enforce precision rules
 * Wn = double(Wn);
 * Arguments    : double n
 *                const double Wn[2]
 *                emxArray_real_T *num
 *                emxArray_creal_T *den
 * Return Type  : void
 */
void butterBandpassOnly(double n, const double Wn[2], emxArray_real_T *num,
  emxArray_creal_T *den)
{
  int k;
  emxArray_real_T *y;
  double u[2];
  double b_Wn;
  int i0;
  emxArray_creal_T *p;
  double B;
  emxArray_creal_T *b_p;
  double im;
  emxArray_creal_T *c_p;
  emxArray_creal_T *x;
  creal_T b_y;
  emxArray_real_T *r;
  emxArray_creal_T *kern;
  emxArray_creal_T *a;
  emxArray_creal_T *b_a;
  double u_data[1];
  creal_T d;
  emxArray_creal_T *e;
  emxArray_int8_T *b_r;
  int j;
  emxArray_creal_T *d_p;
  int end;
  emxArray_int32_T *r0;
  double kern_re;
  double p_re;
  emxArray_int32_T *r1;
  double b_d;
  double c_d;

  /*  step 1: get analog, pre-warped frequencies */
  for (k = 0; k < 2; k++) {
    u[k] = 4.0 * tan(3.1415926535897931 * Wn[k] / 2.0);
  }

  emxInit_real_T1(&y, 2);

  /*  step 2: convert to low-pass prototype estimate */
  b_Wn = sqrt(u[0] * u[1]);

  /*  center frequency */
  /*  Cast to enforce precision rules */
  /* validateattributes(n,{'numeric'},{'scalar','integer','positive'},'butter','N'); */
  /* n = double(n); */
  /*  step 3: Get N-th order Butterworth analog lowpass prototype */
  /* BUTTAP Butterworth analog lowpass filter prototype. */
  /*    [Z,P,K] = BUTTAP(N) returns the zeros, poles, and gain */
  /*    for an N-th order normalized prototype Butterworth analog */
  /*    lowpass filter.  The resulting filter has N poles around */
  /*    the unit circle in the left half plane, and no zeros. */
  /*  */
  /*    % Example: */
  /*    %   Design a 9th order Butterworth analog lowpass filter and display */
  /*    %   its frequency response. */
  /*  */
  /*    [z,p,k]=buttap(9);          % Butterworth filter prototype */
  /*    [num,den]=zp2tf(z,p,k);     % Convert to transfer function form */
  /*    freqs(num,den)              % Frequency response of analog filter           */
  /*  */
  /*    See also BUTTER, CHEB1AP, CHEB2AP, ELLIPAP. */
  /*    Author(s): J.N. Little and J.O. Smith, 1-14-87 */
  /*    	   L. Shure, 1-13-88, revised */
  /*    Copyright 1988-2002 The MathWorks, Inc. */
  /* validateattributes(n,{'numeric'},{'scalar','integer','positive'},'buttap','N'); */
  /*  Cast to enforce precision rules */
  /*  Poles are on the unit circle in the left-half plane. */
  if (rtIsNaN(n - 1.0)) {
    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity_real_T(y, i0);
    y->data[0] = rtNaN;
  } else if (n - 1.0 < 1.0) {
    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 0;
    emxEnsureCapacity_real_T(y, i0);
  } else if (rtIsInf(n - 1.0) && (1.0 == n - 1.0)) {
    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity_real_T(y, i0);
    y->data[0] = rtNaN;
  } else {
    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)floor(((n - 1.0) - 1.0) / 2.0) + 1;
    emxEnsureCapacity_real_T(y, i0);
    k = (int)floor(((n - 1.0) - 1.0) / 2.0);
    for (i0 = 0; i0 <= k; i0++) {
      y->data[y->size[0] * i0] = 1.0 + 2.0 * (double)i0;
    }
  }

  emxInit_creal_T(&p, 2);
  B = 2.0 * n;
  i0 = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = y->size[1];
  emxEnsureCapacity_creal_T(p, i0);
  k = y->size[0] * y->size[1];
  for (i0 = 0; i0 < k; i0++) {
    im = 3.1415926535897931 * y->data[i0];
    p->data[i0].re = (im / B + 1.5707963267948966) * 0.0;
    p->data[i0].im = im / B + 1.5707963267948966;
  }

  emxInit_creal_T(&b_p, 2);
  b_exp(p);
  i0 = b_p->size[0] * b_p->size[1];
  b_p->size[0] = 2;
  b_p->size[1] = p->size[1];
  emxEnsureCapacity_creal_T(b_p, i0);
  k = p->size[1];
  for (i0 = 0; i0 < k; i0++) {
    b_p->data[b_p->size[0] * i0] = p->data[p->size[0] * i0];
  }

  k = p->size[1];
  for (i0 = 0; i0 < k; i0++) {
    b_p->data[1 + b_p->size[0] * i0].re = p->data[p->size[0] * i0].re;
    b_p->data[1 + b_p->size[0] * i0].im = -p->data[p->size[0] * i0].im;
  }

  emxInit_creal_T1(&c_p, 1);
  i0 = c_p->size[0];
  c_p->size[0] = b_p->size[1] << 1;
  emxEnsureCapacity_creal_T1(c_p, i0);
  k = b_p->size[1] << 1;
  for (i0 = 0; i0 < k; i0++) {
    c_p->data[i0] = b_p->data[i0];
  }

  if (rt_remd_snf(n, 2.0) == 1.0) {
    /*  n is odd */
    i0 = c_p->size[0];
    c_p->size[0] = (b_p->size[1] << 1) + 1;
    emxEnsureCapacity_creal_T1(c_p, i0);
    k = b_p->size[1] << 1;
    for (i0 = 0; i0 < k; i0++) {
      c_p->data[i0] = b_p->data[i0];
    }

    c_p->data[b_p->size[1] << 1].re = -1.0;
    c_p->data[b_p->size[1] << 1].im = 0.0;
  }

  emxFree_creal_T(&b_p);
  emxInit_creal_T1(&x, 1);
  i0 = x->size[0];
  x->size[0] = c_p->size[0];
  emxEnsureCapacity_creal_T1(x, i0);
  k = c_p->size[0];
  for (i0 = 0; i0 < k; i0++) {
    x->data[i0].re = -c_p->data[i0].re;
    x->data[i0].im = -c_p->data[i0].im;
  }

  if (x->size[0] == 0) {
    b_y.re = 1.0;
    b_y.im = 0.0;
  } else {
    b_y = x->data[0];
    for (k = 2; k <= x->size[0]; k++) {
      B = b_y.re;
      b_y.re = b_y.re * x->data[k - 1].re - b_y.im * x->data[k - 1].im;
      b_y.im = B * x->data[k - 1].im + b_y.im * x->data[k - 1].re;
    }
  }

  emxInit_real_T(&r, 1);
  emxInit_creal_T(&kern, 2);
  emxInit_creal_T(&a, 2);
  emxInit_creal_T(&b_a, 2);

  /*  z = []; */
  /*  p = exp(1i*(pi*(1:2:n-1)/(2*n) + pi/2)); */
  /*  p = [p; conj(p)]; */
  /*  p = p(:); */
  /*  if rem(n,2)==1   % n is odd */
  /*      p = [p; -1]; */
  /*  end */
  /*  k = real(prod(-p)); */
  /*  Transform to state-space */
  myzp2ss(c_p, b_y.re, a, r, y, &B);

  /*  step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn */
  u_data[0] = u[1] - u[0];
  mylp2bs(a, r, y, B, b_Wn, u_data, b_a, c_p, p, &b_y);

  /*  step 5: Use Bilinear transformation to find discrete equivalent: */
  mybilinear(b_a, c_p, p, b_y, a, x, kern, &d);

  /*  nargout <= 2 */
  mypoly(a, den);

  /*  buttnum: num = buttnum(btype,n,Wn,Bw,analog,den); */
  b_Wn = 2.0 * rt_atan2d_snf(b_Wn, 4.0);
  i0 = r->size[0];
  r->size[0] = (int)n + (int)n;
  emxEnsureCapacity_real_T1(r, i0);
  k = (int)n;
  emxFree_creal_T(&b_a);
  emxFree_creal_T(&a);
  for (i0 = 0; i0 < k; i0++) {
    r->data[i0] = 1.0;
  }

  k = (int)n;
  for (i0 = 0; i0 < k; i0++) {
    r->data[i0 + (int)n] = -1.0;
  }

  /* POLY Convert roots to polynomial. */
  /*    POLY(A), when A is an N by N matrix, is a row vector with */
  /*    N+1 elements which are the coefficients of the */
  /*    characteristic polynomial, det(lambda*eye(size(A)) - A). */
  /*  */
  /*    POLY(V), when V is a vector, is a vector whose elements are */
  /*    the coefficients of the polynomial whose roots are the */
  /*    elements of V. For vectors, ROOTS and POLY are inverse */
  /*    functions of each other, up to ordering, scaling, and */
  /*    roundoff error. */
  /*  */
  /*    Examples: */
  /*  */
  /*    roots(poly(1:20)) generates Wilkinson's famous example. */
  /*  */
  /*    Class support for inputs A,V: */
  /*       float: double, single */
  /*  */
  /*    See also ROOTS, CONV, RESIDUE, POLYVAL. */
  /*    Copyright 1984-2014 The MathWorks, Inc. */
  /*  if m==n square matrix */
  emxInit_creal_T1(&e, 1);
  if (r->size[0] == 1) {
    /*  Characteristic polynomial (square x) */
    emxInit_int8_T(&b_r, 1);
    i0 = b_r->size[0];
    b_r->size[0] = r->size[0];
    emxEnsureCapacity_int8_T(b_r, i0);
    k = r->size[0];
    for (i0 = 0; i0 < k; i0++) {
      b_r->data[i0] = (signed char)r->data[i0];
    }

    i0 = e->size[0];
    e->size[0] = 1;
    emxEnsureCapacity_creal_T1(e, i0);
    e->data[0].re = b_r->data[0];
    e->data[0].im = 0.0;
    emxFree_int8_T(&b_r);

    /* elseif (m==1) || (n==1) */
  } else {
    i0 = e->size[0];
    e->size[0] = r->size[0];
    emxEnsureCapacity_creal_T1(e, i0);
    k = r->size[0];
    for (i0 = 0; i0 < k; i0++) {
      e->data[i0].re = r->data[i0];
      e->data[i0].im = 0.0;
    }
  }

  emxFree_real_T(&r);

  /*  Strip out infinities */
  /*  e = e( isfinite(e) ); */
  /*  Expand recursion formula */
  k = e->size[0];
  i0 = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = 1 + k;
  emxEnsureCapacity_creal_T(p, i0);
  p->data[0].re = 1.0;
  p->data[0].im = 0.0;
  for (i0 = 0; i0 < k; i0++) {
    p->data[p->size[0] * (i0 + 1)].re = 0.0;
    p->data[p->size[0] * (i0 + 1)].im = 0.0;
  }

  /* c = [1 zeros(1,n,class(x))]; */
  /* c = [1 zeros(1,n)]; */
  j = 0;
  emxInit_creal_T(&d_p, 2);
  while (j <= e->size[0] - 1) {
    k = (int)((1.0 + (double)j) + 1.0) - 2;
    B = e->data[j].re;
    im = e->data[j].im;
    i0 = d_p->size[0] * d_p->size[1];
    d_p->size[0] = 1;
    d_p->size[1] = (int)((1.0 + (double)j) + 1.0) - 1;
    emxEnsureCapacity_creal_T(d_p, i0);
    for (i0 = 0; i0 <= k; i0++) {
      kern_re = B * p->data[i0].re - im * p->data[i0].im;
      p_re = B * p->data[i0].im + im * p->data[i0].re;
      d_p->data[d_p->size[0] * i0].re = p->data[1 + i0].re - kern_re;
      d_p->data[d_p->size[0] * i0].im = p->data[1 + i0].im - p_re;
    }

    k = d_p->size[1];
    for (i0 = 0; i0 < k; i0++) {
      p->data[1 + i0] = d_p->data[d_p->size[0] * i0];
    }

    j++;
  }

  emxFree_creal_T(&d_p);

  /*  The result should be real if the roots are complex conjugates. */
  end = e->size[0] - 1;
  k = 0;
  for (j = 0; j <= end; j++) {
    if (e->data[j].im > 0.0) {
      k++;
    }
  }

  emxInit_int32_T(&r0, 1);
  i0 = r0->size[0];
  r0->size[0] = k;
  emxEnsureCapacity_int32_T(r0, i0);
  k = 0;
  for (j = 0; j <= end; j++) {
    if (e->data[j].im > 0.0) {
      r0->data[k] = j + 1;
      k++;
    }
  }

  end = e->size[0] - 1;
  k = 0;
  for (j = 0; j <= end; j++) {
    if (e->data[j].im < 0.0) {
      k++;
    }
  }

  emxInit_int32_T(&r1, 1);
  i0 = r1->size[0];
  r1->size[0] = k;
  emxEnsureCapacity_int32_T(r1, i0);
  k = 0;
  for (j = 0; j <= end; j++) {
    if (e->data[j].im < 0.0) {
      r1->data[k] = j + 1;
      k++;
    }
  }

  i0 = x->size[0];
  x->size[0] = r0->size[0];
  emxEnsureCapacity_creal_T1(x, i0);
  k = r0->size[0];
  for (i0 = 0; i0 < k; i0++) {
    x->data[i0] = e->data[r0->data[i0] - 1];
  }

  emxFree_int32_T(&r0);
  sortButter(x);
  i0 = c_p->size[0];
  c_p->size[0] = r1->size[0];
  emxEnsureCapacity_creal_T1(c_p, i0);
  k = r1->size[0];
  for (i0 = 0; i0 < k; i0++) {
    c_p->data[i0].re = e->data[r1->data[i0] - 1].re;
    c_p->data[i0].im = -e->data[r1->data[i0] - 1].im;
  }

  emxFree_int32_T(&r1);
  emxFree_creal_T(&e);
  sortButter(c_p);
  if (isequal(x, c_p)) {
    i0 = p->size[0] * p->size[1];
    p->size[0] = 1;
    emxEnsureCapacity_creal_T(p, i0);
    k = p->size[0];
    j = p->size[1];
    k *= j;
    for (i0 = 0; i0 < k; i0++) {
      p_re = p->data[i0].re;
      p->data[i0].re = p_re;
      p->data[i0].im = 0.0;
    }
  }

  emxFree_creal_T(&x);
  emxFree_creal_T(&c_p);

  /*  now normalize so |H(w)| == 1: */
  b_y.re = b_Wn * 0.0;
  b_y.im = -b_Wn;
  k = p->size[1];
  i0 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)((double)k - 1.0) + 1;
  emxEnsureCapacity_real_T(y, i0);
  k = (int)((double)k - 1.0);
  for (i0 = 0; i0 <= k; i0++) {
    y->data[y->size[0] * i0] = i0;
  }

  i0 = kern->size[0] * kern->size[1];
  kern->size[0] = 1;
  kern->size[1] = y->size[1];
  emxEnsureCapacity_creal_T(kern, i0);
  k = y->size[0] * y->size[1];
  for (i0 = 0; i0 < k; i0++) {
    kern->data[i0].re = y->data[i0] * b_y.re;
    kern->data[i0].im = y->data[i0] * b_y.im;
  }

  emxFree_real_T(&y);
  b_exp(kern);
  if ((kern->size[1] == 1) || (den->size[1] == 1)) {
    B = 0.0;
    im = 0.0;
    for (i0 = 0; i0 < kern->size[1]; i0++) {
      kern_re = kern->data[kern->size[0] * i0].re * den->data[i0].re -
        kern->data[kern->size[0] * i0].im * den->data[i0].im;
      p_re = kern->data[kern->size[0] * i0].re * den->data[i0].im + kern->
        data[kern->size[0] * i0].im * den->data[i0].re;
      B += kern_re;
      im += p_re;
    }

    b_y.re = B;
    b_y.im = im;
  } else {
    B = 0.0;
    im = 0.0;
    for (i0 = 0; i0 < kern->size[1]; i0++) {
      kern_re = kern->data[kern->size[0] * i0].re * den->data[i0].re -
        kern->data[kern->size[0] * i0].im * den->data[i0].im;
      p_re = kern->data[kern->size[0] * i0].re * den->data[i0].im + kern->
        data[kern->size[0] * i0].im * den->data[i0].re;
      B += kern_re;
      im += p_re;
    }

    b_y.re = B;
    b_y.im = im;
  }

  if ((kern->size[1] == 1) || (p->size[1] == 1)) {
    B = 0.0;
    im = 0.0;
    for (i0 = 0; i0 < kern->size[1]; i0++) {
      kern_re = kern->data[kern->size[0] * i0].re * p->data[i0].re - kern->
        data[kern->size[0] * i0].im * p->data[i0].im;
      p_re = kern->data[kern->size[0] * i0].re * p->data[i0].im + kern->
        data[kern->size[0] * i0].im * p->data[i0].re;
      B += kern_re;
      im += p_re;
    }

    d.re = B;
    d.im = im;
  } else {
    B = 0.0;
    im = 0.0;
    for (i0 = 0; i0 < kern->size[1]; i0++) {
      kern_re = kern->data[kern->size[0] * i0].re * p->data[i0].re - kern->
        data[kern->size[0] * i0].im * p->data[i0].im;
      p_re = kern->data[kern->size[0] * i0].re * p->data[i0].im + kern->
        data[kern->size[0] * i0].im * p->data[i0].re;
      B += kern_re;
      im += p_re;
    }

    d.re = B;
    d.im = im;
  }

  emxFree_creal_T(&kern);
  i0 = num->size[0] * num->size[1];
  num->size[0] = 1;
  num->size[1] = p->size[1];
  emxEnsureCapacity_real_T(num, i0);
  k = p->size[0] * p->size[1];
  for (i0 = 0; i0 < k; i0++) {
    p_re = p->data[i0].re * b_y.re - p->data[i0].im * b_y.im;
    kern_re = p->data[i0].re * b_y.im + p->data[i0].im * b_y.re;
    if (d.im == 0.0) {
      if (kern_re == 0.0) {
        p_re /= d.re;
      } else if (p_re == 0.0) {
        p_re = 0.0;
      } else {
        p_re /= d.re;
      }
    } else if (d.re == 0.0) {
      if (p_re == 0.0) {
        p_re = kern_re / d.im;
      } else if (kern_re == 0.0) {
        p_re = 0.0;
      } else {
        p_re = kern_re / d.im;
      }
    } else {
      B = fabs(d.re);
      im = fabs(d.im);
      if (B > im) {
        B = d.im / d.re;
        p_re = (p_re + B * kern_re) / (d.re + B * d.im);
      } else if (im == B) {
        if (d.re > 0.0) {
          b_d = 0.5;
        } else {
          b_d = -0.5;
        }

        if (d.im > 0.0) {
          c_d = 0.5;
        } else {
          c_d = -0.5;
        }

        p_re = (p_re * b_d + kern_re * c_d) / B;
      } else {
        B = d.re / d.im;
        p_re = (B * p_re + kern_re) / (d.im + B * d.re);
      }
    }

    num->data[i0] = p_re;
  }

  emxFree_creal_T(&p);

  /*  num = poly(a-b*c)+(d-1)*den; */
}

/*
 * File trailer for butterBandpassOnly.c
 *
 * [EOF]
 */
