/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: buttordBandpassOnly.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 11-May-2018 15:51:52
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "buttordBandpassOnly.h"
#include "sortButtord.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */

//int main(void)
//{
//	double wp[] = { 0.3781, 0.4031 };
//	double ws[] = { 0.3750, 0.4063 };
//	double rp = 2, rs = 5;
//	double n = 1;
//	double* order = (double *)malloc(1*sizeof(double));
//	double wn[2];
//	buttordBandpassOnly(wp, ws, rp, rs, order, wn);
//	free(order);
//	order = NULL;
//	return 0;
//}
/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = 1.0;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}
/*
 * Cast to enforce precision rules
 *  wp = signal.internal.sigcasttofloat(wp,'double','buttord','Wp',...
 *      'allownumeric');
 *  ws = signal.internal.sigcasttofloat(ws,'double','buttord','Ws',...
 *      'allownumeric');
 *  rp = signal.internal.sigcasttofloat(rp,'double','buttord','Rp',...
 *      'allownumeric');
 *  rs = signal.internal.sigcasttofloat(rs,'double','buttord','Rs',...
 *      'allownumeric');
 * Arguments    : const double wp[2]
 *                const double ws[2]
 *                double rp
 *                double rs
 *                double *order
 *                double wn[2]
 * Return Type  : void
 */
void buttordBandpassOnly(const double wp[2], const double ws[2], double rp,
  double rs, double *order, double wn[2])
{
  int ixstart;
  double b;
  double WP[2];
  double W0;
  double WA[2];
  double WS[2];
  int ix;
  boolean_T exitg1;
  double c;

  /*  first, prewarp frequencies from digital (unit circle) to analog (imag. axis): */
  /*  next, transform to low pass prototype with passband edge of 1 and stopband */
  /*  edges determined by the following: (see Rabiner and Gold, p.258) */
  for (ixstart = 0; ixstart < 2; ixstart++) {
    W0 = tan(3.1415926535897931 * ws[ixstart] / 2.0);
    WA[ixstart] = W0 * W0;
    WP[ixstart] = tan(3.1415926535897931 * wp[ixstart] / 2.0);
    WS[ixstart] = W0;
  }

  b = WP[0] - WP[1];
  W0 = WP[0] * WP[1];
  for (ixstart = 0; ixstart < 2; ixstart++) {
    WA[ixstart] = (WA[ixstart] - W0) / (WS[ixstart] * b);
  }

  /*  find the minimum order b'worth filter to meet the more demanding spec: */
  for (ixstart = 0; ixstart < 2; ixstart++) {
    WS[ixstart] = fabs(WA[ixstart]);
  }

  ixstart = 1;
  W0 = WS[0];
  if (rtIsNaN(WS[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 3)) {
      ixstart = 2;
      if (!rtIsNaN(WS[1])) {
        W0 = WS[1];
        exitg1 = true;
      } else {
        ix = 3;
      }
    }
  }

  if ((ixstart < 2) && (WS[1] < W0)) {
    W0 = WS[1];
  }

  *order = ceil(log10((rt_powd_snf(10.0, 0.1 * fabs(rs)) - 1.0) / (rt_powd_snf
    (10.0, 0.1 * fabs(rp)) - 1.0)) / (2.0 * log10(W0)));

  /*  next find the butterworth natural frequency W0 (or, the "3dB frequency") */
  /*  to give exactly rs dB at WA.  W0 will be between 1 and WA: */
  W0 /= rt_powd_snf(rt_powd_snf(10.0, 0.1 * fabs(rs)) - 1.0, 1.0 / (2.0 * fabs
    (*order)));

  /*  now convert this frequency back from lowpass prototype */
  /*  to the original analog filter: */
  wn[0] = -W0;
  wn[1] = W0;

  /*  need both left and right 3dB frequencies */
  b = WP[1] - WP[0];
  W0 = WP[1] - WP[0];
  c = W0 * W0;
  W0 = WP[0] * WP[1];
  for (ixstart = 0; ixstart < 2; ixstart++) {
    wn[ixstart] = -wn[ixstart] * b / 2.0 + sqrt(wn[ixstart] * wn[ixstart] / 4.0 *
      c + W0);
  }

  for (ixstart = 0; ixstart < 2; ixstart++) {
    wn[ixstart] = fabs(wn[ixstart]);
  }

  sortButtord(wn);

  /*  finally, transform frequencies from analog to digital if necessary: */
  for (ixstart = 0; ixstart < 2; ixstart++) {
    wn[ixstart] = 0.63661977236758138 * atan(wn[ixstart]);
  }

  /*  bilinear transform */
}

/*
 * File trailer for buttordBandpassOnly.c
 *
 * [EOF]
 */
