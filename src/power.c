/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: power.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 09-Aug-2018 12:35:13
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "sphereToimage.h"
#include "power.h"
#include "log.h"
#include "sphereToimage_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T a
 * Return Type  : creal_T
 */
creal_T b_power(const creal_T a)
{
  creal_T y;
  creal_T dc7;
  double r;
  if ((a.im == 0.0) && (a.re >= 0.0)) {
    y.re = rt_powd_snf(a.re, -2.0);
    y.im = 0.0;
  } else if (a.re == 0.0) {
    y.re = -rt_powd_snf(a.im, -2.0);
    y.im = 0.0;
  } else {
    dc7 = a;
    b_log(&dc7);
    y.re = -2.0 * dc7.re;
    y.im = -2.0 * dc7.im;
    if (y.im == 0.0) {
      y.re = exp(y.re);
      y.im = 0.0;
    } else if (rtIsInf(y.im) && rtIsInf(y.re) && (y.re < 0.0)) {
      y.re = 0.0;
      y.im = 0.0;
    } else {
      r = exp(y.re / 2.0);
      y.re = r * (r * cos(y.im));
      y.im = r * (r * sin(y.im));
    }
  }

  return y;
}

/*
 * Arguments    : const creal_T a
 * Return Type  : creal_T
 */
creal_T c_power(const creal_T a)
{
  creal_T y;
  creal_T dc8;
  double r;
  if ((a.im == 0.0) && (a.re >= 0.0)) {
    y.re = rt_powd_snf(a.re, 3.0);
    y.im = 0.0;
  } else if (a.re == 0.0) {
    y.re = 0.0;
    y.im = -rt_powd_snf(a.im, 3.0);
  } else {
    dc8 = a;
    b_log(&dc8);
    y.re = 3.0 * dc8.re;
    y.im = 3.0 * dc8.im;
    if (y.im == 0.0) {
      y.re = exp(y.re);
      y.im = 0.0;
    } else if (rtIsInf(y.im) && rtIsInf(y.re) && (y.re < 0.0)) {
      y.re = 0.0;
      y.im = 0.0;
    } else {
      r = exp(y.re / 2.0);
      y.re = r * (r * cos(y.im));
      y.im = r * (r * sin(y.im));
    }
  }

  return y;
}

/*
 * Arguments    : const creal_T a
 * Return Type  : creal_T
 */
creal_T d_power(const creal_T a)
{
  creal_T y;
  double absxi;
  double absxr;
  if (a.im == 0.0) {
    if (a.re < 0.0) {
      absxi = 0.0;
      absxr = sqrt(-a.re);
    } else {
      absxi = sqrt(a.re);
      absxr = 0.0;
    }
  } else if (a.re == 0.0) {
    if (a.im < 0.0) {
      absxi = sqrt(-a.im / 2.0);
      absxr = -absxi;
    } else {
      absxi = sqrt(a.im / 2.0);
      absxr = absxi;
    }
  } else if (rtIsNaN(a.re)) {
    absxi = a.re;
    absxr = a.re;
  } else if (rtIsNaN(a.im)) {
    absxi = a.im;
    absxr = a.im;
  } else if (rtIsInf(a.im)) {
    absxi = fabs(a.im);
    absxr = a.im;
  } else if (rtIsInf(a.re)) {
    if (a.re < 0.0) {
      absxi = 0.0;
      absxr = a.im * -a.re;
    } else {
      absxi = a.re;
      absxr = 0.0;
    }
  } else {
    absxr = fabs(a.re);
    absxi = fabs(a.im);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = sqrt(absxi) * sqrt(1.0 + absxr / absxi);
      } else {
        absxi = sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxi = sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (a.re > 0.0) {
      absxr = 0.5 * (a.im / absxi);
    } else {
      if (a.im < 0.0) {
        absxr = -absxi;
      } else {
        absxr = absxi;
      }

      absxi = 0.5 * (a.im / absxr);
    }
  }

  y.re = absxi;
  y.im = absxr;
  return y;
}

/*
 * Arguments    : const creal_T a
 * Return Type  : creal_T
 */
creal_T e_power(const creal_T a)
{
  creal_T y;
  creal_T dc9;
  double r;
  if ((a.im == 0.0) && (a.re >= 0.0)) {
    y.re = rt_powd_snf(a.re, -0.33333333333333331);
    y.im = 0.0;
  } else {
    dc9 = a;
    b_log(&dc9);
    y.re = -0.33333333333333331 * dc9.re;
    y.im = -0.33333333333333331 * dc9.im;
    if (y.im == 0.0) {
      y.re = exp(y.re);
      y.im = 0.0;
    } else if (rtIsInf(y.im) && rtIsInf(y.re) && (y.re < 0.0)) {
      y.re = 0.0;
      y.im = 0.0;
    } else {
      r = exp(y.re / 2.0);
      y.re = r * (r * cos(y.im));
      y.im = r * (r * sin(y.im));
    }
  }

  return y;
}

/*
 * Arguments    : const creal_T a
 * Return Type  : creal_T
 */
creal_T f_power(const creal_T a)
{
  creal_T y;
  creal_T dc10;
  double r;
  if ((a.im == 0.0) && (a.re >= 0.0)) {
    y.re = rt_powd_snf(a.re, 0.33333333333333331);
    y.im = 0.0;
  } else {
    dc10 = a;
    b_log(&dc10);
    y.re = 0.33333333333333331 * dc10.re;
    y.im = 0.33333333333333331 * dc10.im;
    if (y.im == 0.0) {
      y.re = exp(y.re);
      y.im = 0.0;
    } else if (rtIsInf(y.im) && rtIsInf(y.re) && (y.re < 0.0)) {
      y.re = 0.0;
      y.im = 0.0;
    } else {
      r = exp(y.re / 2.0);
      y.re = r * (r * cos(y.im));
      y.im = r * (r * sin(y.im));
    }
  }

  return y;
}

/*
 * Arguments    : const creal_T a
 * Return Type  : creal_T
 */
creal_T g_power(const creal_T a)
{
  creal_T y;
  creal_T dc11;
  double r;
  if ((a.im == 0.0) && (a.re >= 0.0)) {
    y.re = rt_powd_snf(a.re, -3.0);
    y.im = 0.0;
  } else if (a.re == 0.0) {
    y.re = 0.0;
    y.im = rt_powd_snf(a.im, -3.0);
  } else {
    dc11 = a;
    b_log(&dc11);
    y.re = -3.0 * dc11.re;
    y.im = -3.0 * dc11.im;
    if (y.im == 0.0) {
      y.re = exp(y.re);
      y.im = 0.0;
    } else if (rtIsInf(y.im) && rtIsInf(y.re) && (y.re < 0.0)) {
      y.re = 0.0;
      y.im = 0.0;
    } else {
      r = exp(y.re / 2.0);
      y.re = r * (r * cos(y.im));
      y.im = r * (r * sin(y.im));
    }
  }

  return y;
}

/*
 * Arguments    : const creal_T a
 * Return Type  : creal_T
 */
creal_T power(const creal_T a)
{
  creal_T y;
  double brm;
  double bim;
  double d;
  brm = fabs(a.re);
  bim = fabs(a.im);
  if (a.im == 0.0) {
    y.re = 1.0 / a.re;
    y.im = 0.0;
  } else if (a.re == 0.0) {
    y.re = 0.0;
    y.im = -1.0 / a.im;
  } else if (brm > bim) {
    bim = a.im / a.re;
    d = a.re + bim * a.im;
    y.re = 1.0 / d;
    y.im = -bim / d;
  } else if (brm == bim) {
    bim = 0.5;
    if (a.re < 0.0) {
      bim = -0.5;
    }

    d = 0.5;
    if (a.im < 0.0) {
      d = -0.5;
    }

    y.re = bim / brm;
    y.im = -d / brm;
  } else {
    bim = a.re / a.im;
    d = a.im + bim * a.re;
    y.re = bim / d;
    y.im = -1.0 / d;
  }

  return y;
}

/*
 * File trailer for power.c
 *
 * [EOF]
 */
