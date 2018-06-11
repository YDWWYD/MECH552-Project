/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mybilinear.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 14-May-2018 23:15:09
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
#include "mybilinear.h"
#include "butterBandpassOnly_emxutil.h"
#include "mrdivide.h"
#include "xtrsm.h"
#include "xgetrf.h"
#include "qrsolve.h"
#include "xgeqp3.h"
#include "eye.h"

/* Function Definitions */

/*
 * BILINEAR Bilinear transformation with optional frequency prewarping.
 *    [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs) converts the s-domain transfer
 *    function specified by Z, P, and K to a z-transform discrete
 *    equivalent obtained from the bilinear transformation:
 *
 *       H(z) = H(s) |
 *                   | s = 2*Fs*(z-1)/(z+1)
 *
 *    where column vectors Z and P specify the zeros and poles, scalar
 *    K specifies the gain, and Fs is the sample frequency in Hz.
 *
 *    [NUMd,DENd] = BILINEAR(NUM,DEN,Fs), where NUM and DEN are
 *    row vectors containing numerator and denominator transfer
 *    function coefficients, NUM(s)/DEN(s), in descending powers of
 *    s, transforms to z-transform coefficients NUMd(z)/DENd(z).
 *
 *    [Ad,Bd,Cd,Dd] = BILINEAR(A,B,C,D,Fs) is a state-space version.
 *
 *    Each of the above three forms of BILINEAR accepts an optional
 *    additional input argument that specifies prewarping.
 *
 *    For example, [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs,Fp) applies prewarping
 *    before the bilinear transformation so that the frequency responses
 *    before and after mapping match exactly at frequency point Fp
 *    (match point Fp is specified in Hz).
 *
 *    % Example:
 *    %   Design a 6-th order Elliptic analog low pass filter and transform
 *    %   it to a Discrete-time representation.
 *
 *    Fs =0.5;                            % Sampling Frequency
 *    [z,p,k]=ellipap(6,5,90);            % Lowpass filter prototype
 *    [num,den]=zp2tf(z,p,k);             % Convert to transfer function form
 *    [numd,dend]=bilinear(num,den,Fs);   % Analog to Digital conversion
 *    fvtool(numd,dend)                   % Visualize the filter
 *
 *    See also IMPINVAR.
 * Arguments    : const emxArray_creal_T *z
 *                const emxArray_creal_T *p
 *                const emxArray_creal_T *k
 *                const creal_T fs
 *                emxArray_creal_T *zd
 *                emxArray_creal_T *pd
 *                emxArray_creal_T *kd
 *                creal_T *dd
 * Return Type  : void
 */
void mybilinear(const emxArray_creal_T *z, const emxArray_creal_T *p, const
                emxArray_creal_T *k, const creal_T fs, emxArray_creal_T *zd,
                emxArray_creal_T *pd, emxArray_creal_T *kd, creal_T *dd)
{
  int m;
  emxArray_creal_T *t1;
  double b_z[2];
  emxArray_real_T *r4;
  int jBcol;
  double z_im;
  double wj_im;
  emxArray_creal_T *t2;
  double wj_re;
  emxArray_creal_T *A;
  emxArray_creal_T *tau;
  emxArray_int32_T *jpvt;
  unsigned int unnamed_idx_0;
  int mn;
  unsigned int unnamed_idx_1;
  int rankA;
  int nb;
  int b_nb;
  int j;
  emxArray_creal_T *a;
  int b_k;
  double tauj_re;
  double tauj_im;
  boolean_T b_zd;
  int i;
  double A_re;
  double A_im;

  /*    Author(s): J.N. Little, 4-28-87 */
  /*    	   J.N. Little, 5-5-87, revised */
  /*    Copyright 1988-2006 The MathWorks, Inc. */
  /*    Gene Franklin, Stanford Univ., motivated the state-space */
  /*    approach to the bilinear transformation. */
  /* [mn,nn] = size(z); */
  /* [md,nd] = size(p); */
  /*  if (nd == 1 && nn < 2) && nargout ~= 4	% In zero-pole-gain form */
  /*      if mn > md */
  /*          error(message('signal:bilinear:InvalidRange')) */
  /*      end */
  /*      if nargin == 5 */
  /*        % Cast to enforce Precision Rules */
  /*        fp = signal.internal.sigcasttofloat(fp,'double','bilinear','',... */
  /*          'allownumeric'); */
  /*        fs = signal.internal.sigcasttofloat(fs,'double','bilinear','',... */
  /*          'allownumeric'); */
  /*        % Prewarp */
  /*        fp = 2*pi*fp; */
  /*        fs = fp/tan(fp/fs/2); */
  /*      else */
  /*        % Cast to enforce Precision Rules */
  /*        fs = 2*signal.internal.sigcasttofloat(fs,'double','bilinear','',... */
  /*          'allownumeric'); */
  /*      end */
  /*      % Cast to enforce Precision Rules */
  /*      if any([signal.internal.sigcheckfloattype(z,'single','bilinear','Z') ... */
  /*          signal.internal.sigcheckfloattype(p,'single','bilinear','P')... */
  /*          signal.internal.sigcheckfloattype(k,'single','bilinear','K')]) */
  /*        z = single(z); */
  /*        p = single(p); */
  /*        k = single(k); */
  /*      end */
  /*      z = z(isfinite(z));	 % Strip infinities from zeros */
  /*      pd = (1+p/fs)./(1-p/fs); % Do bilinear transformation */
  /*      zd = (1+z/fs)./(1-z/fs); */
  /*      % real(kd) or just kd? */
  /*      kd = (k*prod(fs-z)./prod(fs-p)); */
  /*      zd = [zd;-ones(length(pd)-length(zd),1)];  % Add extra zeros at -1 */
  /* elseif (md == 1 && mn == 1) || nargout == 4 % */
  /* if nargout == 4		% State-space case */
  /*  Cast to enforce Precision Rules */
  /*          error(abcdchk(a,b,c,d)); */
  /*          if nargin == 6			% Prewarp */
  /*            % Cast to enforce Precision Rules */
  /*            fp = signal.internal.sigcasttofloat(fp1,'double','bilinear','',... */
  /*              'allownumeric');		% Decode arguments */
  /*            fp = 2*pi*fp; */
  /*            fs = fp/tan(fp/fs/2)/2; */
  /*          end */
  /*      else			% Transfer function case */
  /*          if nn > nd */
  /*              error(message('signal:bilinear:InvalidRange')) */
  /*          end */
  /*          num = z; den = p;		% Decode arguments */
  /*          if any([signal.internal.sigcheckfloattype(num,'single','bilinear','NUM')... */
  /*              signal.internal.sigcheckfloattype(den,'single','bilinear','DEN')]) */
  /*            num = single(num); */
  /*            den = single(den); */
  /*          end */
  /*          if nargin == 4			% Prewarp */
  /*              % Cast to enforce Precision Rules */
  /*              fp = signal.internal.sigcasttofloat(fs,'double','bilinear','',... */
  /*                'allownumeric');  */
  /*              fs = signal.internal.sigcasttofloat(k,'double','bilinear','',... */
  /*                'allownumeric');	% Decode arguments */
  /*              fp = 2*pi*fp; */
  /*              fs = fp/tan(fp/fs/2)/2; */
  /*               */
  /*          else */
  /*              % Cast to enforce Precision Rules */
  /*              fs = signal.internal.sigcasttofloat(k,'double','bilinear','',... */
  /*                'allownumeric');	% Decode arguments */
  /*          end */
  /*           */
  /*          % Put num(s)/den(s) in state-space canonical form. */
  /*          [a,b,c,d] = tf2ss(num,den); */
  /*      end */
  /*  Now do state-space version of bilinear transformation: */
  for (m = 0; m < 2; m++) {
    b_z[m] = z->size[m];
  }

  emxInit_creal_T(&t1, 2);
  emxInit_real_T1(&r4, 2);
  b_eye(b_z, r4);
  m = t1->size[0] * t1->size[1];
  t1->size[0] = r4->size[0];
  t1->size[1] = r4->size[1];
  emxEnsureCapacity_creal_T(t1, m);
  jBcol = r4->size[0] * r4->size[1];
  for (m = 0; m < jBcol; m++) {
    z_im = z->data[m].re * 0.5;
    wj_im = z->data[m].im * 0.5;
    if (wj_im == 0.0) {
      wj_re = z_im / 2.0;
      z_im = 0.0;
    } else if (z_im == 0.0) {
      wj_re = 0.0;
      z_im = wj_im / 2.0;
    } else {
      wj_re = z_im / 2.0;
      z_im = wj_im / 2.0;
    }

    t1->data[m].re = r4->data[m] + wj_re;
    t1->data[m].im = z_im;
  }

  for (m = 0; m < 2; m++) {
    b_z[m] = z->size[m];
  }

  emxInit_creal_T(&t2, 2);
  b_eye(b_z, r4);
  m = t2->size[0] * t2->size[1];
  t2->size[0] = r4->size[0];
  t2->size[1] = r4->size[1];
  emxEnsureCapacity_creal_T(t2, m);
  jBcol = r4->size[0] * r4->size[1];
  for (m = 0; m < jBcol; m++) {
    z_im = z->data[m].re * 0.5;
    wj_im = z->data[m].im * 0.5;
    if (wj_im == 0.0) {
      wj_re = z_im / 2.0;
      z_im = 0.0;
    } else if (z_im == 0.0) {
      wj_re = 0.0;
      z_im = wj_im / 2.0;
    } else {
      wj_re = z_im / 2.0;
      z_im = wj_im / 2.0;
    }

    t2->data[m].re = r4->data[m] - wj_re;
    t2->data[m].im = 0.0 - z_im;
  }

  emxFree_real_T(&r4);
  emxInit_creal_T(&A, 2);
  emxInit_creal_T1(&tau, 1);
  emxInit_int32_T1(&jpvt, 2);
  if ((t2->size[0] == 0) || (t2->size[1] == 0) || ((t1->size[0] == 0) ||
       (t1->size[1] == 0))) {
    unnamed_idx_0 = (unsigned int)t2->size[1];
    unnamed_idx_1 = (unsigned int)t1->size[1];
    m = zd->size[0] * zd->size[1];
    zd->size[0] = (int)unnamed_idx_0;
    zd->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity_creal_T(zd, m);
    jBcol = (int)unnamed_idx_0 * (int)unnamed_idx_1;
    for (m = 0; m < jBcol; m++) {
      zd->data[m].re = 0.0;
      zd->data[m].im = 0.0;
    }
  } else if (t2->size[0] == t2->size[1]) {
    mn = t2->size[1];
    m = A->size[0] * A->size[1];
    A->size[0] = t2->size[0];
    A->size[1] = t2->size[1];
    emxEnsureCapacity_creal_T(A, m);
    jBcol = t2->size[0] * t2->size[1];
    for (m = 0; m < jBcol; m++) {
      A->data[m] = t2->data[m];
    }

    xgetrf(t2->size[1], t2->size[1], A, t2->size[1], jpvt, &jBcol);
    nb = t1->size[1];
    m = zd->size[0] * zd->size[1];
    zd->size[0] = t1->size[0];
    zd->size[1] = t1->size[1];
    emxEnsureCapacity_creal_T(zd, m);
    jBcol = t1->size[0] * t1->size[1];
    for (m = 0; m < jBcol; m++) {
      zd->data[m] = t1->data[m];
    }

    for (jBcol = 0; jBcol + 1 < mn; jBcol++) {
      if (jpvt->data[jBcol] != jBcol + 1) {
        b_nb = jpvt->data[jBcol] - 1;
        for (m = 0; m + 1 <= nb; m++) {
          tauj_re = zd->data[jBcol + zd->size[0] * m].re;
          tauj_im = zd->data[jBcol + zd->size[0] * m].im;
          zd->data[jBcol + zd->size[0] * m] = zd->data[b_nb + zd->size[0] * m];
          zd->data[b_nb + zd->size[0] * m].re = tauj_re;
          zd->data[b_nb + zd->size[0] * m].im = tauj_im;
        }
      }
    }

    for (j = 1; j <= nb; j++) {
      jBcol = mn * (j - 1);
      for (b_k = 0; b_k + 1 <= mn; b_k++) {
        b_nb = mn * b_k;
        b_zd = ((zd->data[b_k + jBcol].re != 0.0) || (zd->data[b_k + jBcol].im
                 != 0.0));
        if (b_zd) {
          for (i = b_k + 1; i + 1 <= mn; i++) {
            tauj_re = zd->data[b_k + jBcol].re * A->data[i + b_nb].re - zd->
              data[b_k + jBcol].im * A->data[i + b_nb].im;
            tauj_im = zd->data[b_k + jBcol].re * A->data[i + b_nb].im + zd->
              data[b_k + jBcol].im * A->data[i + b_nb].re;
            zd->data[i + jBcol].re -= tauj_re;
            zd->data[i + jBcol].im -= tauj_im;
          }
        }
      }
    }

    xtrsm(t2->size[1], t1->size[1], A, t2->size[1], zd, t2->size[1]);
  } else {
    m = A->size[0] * A->size[1];
    A->size[0] = t2->size[0];
    A->size[1] = t2->size[1];
    emxEnsureCapacity_creal_T(A, m);
    jBcol = t2->size[0] * t2->size[1];
    for (m = 0; m < jBcol; m++) {
      A->data[m] = t2->data[m];
    }

    xgeqp3(A, tau, jpvt);
    rankA = rankFromQR(A);
    nb = t1->size[1];
    jBcol = A->size[1];
    b_nb = t1->size[1];
    m = zd->size[0] * zd->size[1];
    zd->size[0] = jBcol;
    zd->size[1] = b_nb;
    emxEnsureCapacity_creal_T(zd, m);
    jBcol *= b_nb;
    for (m = 0; m < jBcol; m++) {
      zd->data[m].re = 0.0;
      zd->data[m].im = 0.0;
    }

    m = A->size[0];
    b_nb = t1->size[1];
    jBcol = A->size[0];
    mn = A->size[1];
    if (jBcol < mn) {
      mn = jBcol;
    }

    for (j = 0; j + 1 <= mn; j++) {
      tauj_re = tau->data[j].re;
      tauj_im = -tau->data[j].im;
      if ((tauj_re != 0.0) || (tauj_im != 0.0)) {
        for (b_k = 0; b_k + 1 <= b_nb; b_k++) {
          wj_re = t1->data[j + t1->size[0] * b_k].re;
          wj_im = t1->data[j + t1->size[0] * b_k].im;
          for (i = j + 1; i + 1 <= m; i++) {
            wj_re += A->data[i + A->size[0] * j].re * t1->data[i + t1->size[0] *
              b_k].re + A->data[i + A->size[0] * j].im * t1->data[i + t1->size[0]
              * b_k].im;
            wj_im += A->data[i + A->size[0] * j].re * t1->data[i + t1->size[0] *
              b_k].im - A->data[i + A->size[0] * j].im * t1->data[i + t1->size[0]
              * b_k].re;
          }

          z_im = wj_re;
          wj_re = tauj_re * wj_re - tauj_im * wj_im;
          wj_im = tauj_re * wj_im + tauj_im * z_im;
          if ((wj_re != 0.0) || (wj_im != 0.0)) {
            t1->data[j + t1->size[0] * b_k].re -= wj_re;
            t1->data[j + t1->size[0] * b_k].im -= wj_im;
            for (i = j + 1; i + 1 <= m; i++) {
              A_re = A->data[i + A->size[0] * j].re * wj_re - A->data[i +
                A->size[0] * j].im * wj_im;
              A_im = A->data[i + A->size[0] * j].re * wj_im + A->data[i +
                A->size[0] * j].im * wj_re;
              t1->data[i + t1->size[0] * b_k].re -= A_re;
              t1->data[i + t1->size[0] * b_k].im -= A_im;
            }
          }
        }
      }
    }

    for (b_k = 0; b_k + 1 <= nb; b_k++) {
      for (i = 0; i + 1 <= rankA; i++) {
        zd->data[(jpvt->data[i] + zd->size[0] * b_k) - 1] = t1->data[i +
          t1->size[0] * b_k];
      }

      for (j = rankA - 1; j + 1 > 0; j--) {
        tauj_re = zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re;
        tauj_im = zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im;
        A_re = A->data[j + A->size[0] * j].re;
        A_im = A->data[j + A->size[0] * j].im;
        if (A_im == 0.0) {
          if (tauj_im == 0.0) {
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = tauj_re /
              A_re;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = 0.0;
          } else if (tauj_re == 0.0) {
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = 0.0;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = tauj_im /
              A_re;
          } else {
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = tauj_re /
              A_re;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = tauj_im /
              A_re;
          }
        } else if (A_re == 0.0) {
          if (tauj_re == 0.0) {
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = tauj_im /
              A_im;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = 0.0;
          } else if (tauj_im == 0.0) {
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = 0.0;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = -(tauj_re /
              A_im);
          } else {
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = tauj_im /
              A_im;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = -(tauj_re /
              A_im);
          }
        } else {
          wj_im = fabs(A_re);
          z_im = fabs(A_im);
          if (wj_im > z_im) {
            wj_re = A_im / A_re;
            z_im = A_re + wj_re * A_im;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = (tauj_re +
              wj_re * tauj_im) / z_im;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = (tauj_im -
              wj_re * tauj_re) / z_im;
          } else if (z_im == wj_im) {
            if (A_re > 0.0) {
              wj_re = 0.5;
            } else {
              wj_re = -0.5;
            }

            if (A_im > 0.0) {
              z_im = 0.5;
            } else {
              z_im = -0.5;
            }

            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = (tauj_re *
              wj_re + tauj_im * z_im) / wj_im;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = (tauj_im *
              wj_re - tauj_re * z_im) / wj_im;
          } else {
            wj_re = A_re / A_im;
            z_im = A_im + wj_re * A_re;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re = (wj_re *
              tauj_re + tauj_im) / z_im;
            zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].im = (wj_re *
              tauj_im - tauj_re) / z_im;
          }
        }

        for (i = 0; i + 1 <= j; i++) {
          tauj_re = zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re *
            A->data[i + A->size[0] * j].re - zd->data[(jpvt->data[j] + zd->size
            [0] * b_k) - 1].im * A->data[i + A->size[0] * j].im;
          tauj_im = zd->data[(jpvt->data[j] + zd->size[0] * b_k) - 1].re *
            A->data[i + A->size[0] * j].im + zd->data[(jpvt->data[j] + zd->size
            [0] * b_k) - 1].im * A->data[i + A->size[0] * j].re;
          zd->data[(jpvt->data[i] + zd->size[0] * b_k) - 1].re -= tauj_re;
          zd->data[(jpvt->data[i] + zd->size[0] * b_k) - 1].im -= tauj_im;
        }
      }
    }
  }

  emxFree_creal_T(&tau);
  emxFree_creal_T(&t1);
  if ((t2->size[0] == 0) || (t2->size[1] == 0) || (p->size[0] == 0)) {
    unnamed_idx_0 = (unsigned int)t2->size[1];
    m = pd->size[0];
    pd->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity_creal_T1(pd, m);
    jBcol = (int)unnamed_idx_0;
    for (m = 0; m < jBcol; m++) {
      pd->data[m].re = 0.0;
      pd->data[m].im = 0.0;
    }
  } else if (t2->size[0] == t2->size[1]) {
    m = A->size[0] * A->size[1];
    A->size[0] = t2->size[0];
    A->size[1] = t2->size[1];
    emxEnsureCapacity_creal_T(A, m);
    jBcol = t2->size[0] * t2->size[1];
    for (m = 0; m < jBcol; m++) {
      A->data[m] = t2->data[m];
    }

    xgetrf(t2->size[1], t2->size[1], A, t2->size[1], jpvt, &jBcol);
    m = pd->size[0];
    pd->size[0] = p->size[0];
    emxEnsureCapacity_creal_T1(pd, m);
    jBcol = p->size[0];
    for (m = 0; m < jBcol; m++) {
      pd->data[m] = p->data[m];
    }

    for (jBcol = 0; jBcol + 1 < t2->size[1]; jBcol++) {
      if (jpvt->data[jBcol] != jBcol + 1) {
        tauj_re = pd->data[jBcol].re;
        tauj_im = pd->data[jBcol].im;
        pd->data[jBcol] = pd->data[jpvt->data[jBcol] - 1];
        pd->data[jpvt->data[jBcol] - 1].re = tauj_re;
        pd->data[jpvt->data[jBcol] - 1].im = tauj_im;
      }
    }

    b_xtrsm(t2->size[1], A, t2->size[1], pd);
    c_xtrsm(t2->size[1], A, t2->size[1], pd);
  } else {
    b_qrsolve(t2, p, pd);
  }

  emxFree_int32_T(&jpvt);
  emxFree_creal_T(&A);
  m = pd->size[0];
  emxEnsureCapacity_creal_T1(pd, m);
  jBcol = pd->size[0];
  for (m = 0; m < jBcol; m++) {
    pd->data[m].re *= 0.70710678118654746;
    pd->data[m].im *= 0.70710678118654746;
  }

  emxInit_creal_T(&a, 2);
  m = a->size[0] * a->size[1];
  a->size[0] = 1;
  a->size[1] = k->size[1];
  emxEnsureCapacity_creal_T(a, m);
  jBcol = k->size[0] * k->size[1];
  for (m = 0; m < jBcol; m++) {
    a->data[m].re = 0.70710678118654757 * k->data[m].re;
    a->data[m].im = 0.70710678118654757 * k->data[m].im;
  }

  b_mrdivide(a, t2, kd);
  b_mrdivide(k, t2, a);
  emxFree_creal_T(&t2);
  if ((a->size[1] == 1) || (p->size[0] == 1)) {
    tauj_re = 0.0;
    tauj_im = 0.0;
    for (m = 0; m < a->size[1]; m++) {
      z_im = a->data[a->size[0] * m].re * p->data[m].re - a->data[a->size[0] * m]
        .im * p->data[m].im;
      wj_re = a->data[a->size[0] * m].re * p->data[m].im + a->data[a->size[0] *
        m].im * p->data[m].re;
      tauj_re += z_im;
      tauj_im += wj_re;
    }
  } else {
    tauj_re = 0.0;
    tauj_im = 0.0;
    for (m = 0; m < a->size[1]; m++) {
      z_im = a->data[a->size[0] * m].re * p->data[m].re - a->data[a->size[0] * m]
        .im * p->data[m].im;
      wj_re = a->data[a->size[0] * m].re * p->data[m].im + a->data[a->size[0] *
        m].im * p->data[m].re;
      tauj_re += z_im;
      tauj_im += wj_re;
    }
  }

  emxFree_creal_T(&a);
  z_im = tauj_re * 0.5;
  wj_im = tauj_im * 0.5;
  if (wj_im == 0.0) {
    tauj_re = z_im / 2.0;
    tauj_im = 0.0;
  } else if (z_im == 0.0) {
    tauj_re = 0.0;
    tauj_im = wj_im / 2.0;
  } else {
    tauj_re = z_im / 2.0;
    tauj_im = wj_im / 2.0;
  }

  dd->re = tauj_re + fs.re;
  dd->im = tauj_im + fs.im;

  /* if nargout == 4 */
  /*      else */
  /*          % Convert back to transfer function form: */
  /*          p = poly(ad); */
  /*          zd = poly(ad-bd*cd)+(dd-1)*p; */
  /*          pd = p; */
  /*      end */
  /*  else */
  /*      error(message('signal:bilinear:SignalErr')) */
}

/*
 * File trailer for mybilinear.c
 *
 * [EOF]
 */
