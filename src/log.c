#include <math.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "sphereToimage.h"
#include "log.h"
#include "power.h"
#include "sphereToimage_rtwutil.h"

static double rt_atan2d_snf(double u0, double u1);

static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(b_u0, b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

void b_log(creal_T *x)
{
  double x_im;
  double x_re;
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      x->re = log(fabs(x->re));
      x->im = 3.1415926535897931;
    } else {
      x->re = log(fabs(x->re));
      x->im = 0.0;
    }
  } else if ((fabs(x->re) > 8.9884656743115785E+307) || (fabs(x->im) >
              8.9884656743115785E+307)) {
    x_im = x->im;
    x_re = x->re;
    x->re = log(rt_hypotd_snf(x->re / 2.0, x->im / 2.0)) + 0.69314718055994529;
    x->im = rt_atan2d_snf(x_im, x_re);
  } else {
    x_im = x->im;
    x_re = x->re;
    x->re = log(rt_hypotd_snf(x->re, x->im));
    x->im = rt_atan2d_snf(x_im, x_re);
  }
}