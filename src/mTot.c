/*
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 09-Aug-2018 12:35:13
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "sphereToimage.h"
#include "mTot.h"
#include "power.h"
#include "log.h"
#include "sphereToimage_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : double m
 *                double d
 *                double rm
 *                double R
 *                double del
 * Return Type  : double
 */
double mTot(double m, double d, double rm, double R, double del)
{
  creal_T b_m;
  creal_T b_R;
  double re;
  double b_re;
  double c_re;
  double im;
  double m_re;
  double m_im;
  double b_m_re;
  double b_m_im;
  double b_im;
  double d_re;
  double c_m_re;
  double c_m_im;
  double e_re;
  double f_re;
  double c_im;
  double g_re;
  double d_im;
  double R_re;
  double R_im;
  double e_im;
  double h_re;
  double b_d_re;
  double b_d_im;
  double i_re;
  double j_re;
  double f_im;
  double g_im;
  double k_re;
  double a_re;
  double a_im;
  double b_a_re;
  double b_a_im;
  double c_a_re;
  double c_a_im;
  double d_a_re;
  double d_a_im;
  double e_a_re;
  double e_a_im;
  double f_a_re;
  double f_a_im;
  double g_a_re;
  double g_a_im;
  creal_T b_d;
  creal_T dc0;
  double c_d_re;
  double c_d_im;
  double b_R_re;
  double b_R_im;
  double h_a_re;
  double h_im;
  double d_d_re;
  double d_d_im;
  double l_re;
  double n_re;
  double i_im;
  double j_im;
  double o_re;
  double i_a_re;
  double h_a_im;
  double i_a_im;
  double j_a_re;
  double j_a_im;
  double k_a_re;
  double k_a_im;
  double l_a_re;
  double m_a_re;
  double l_a_im;
  double m_a_im;
  double n_a_re;
  double n_a_im;
  double o_a_re;
  double o_a_im;
  double p_a_re;
  double p_a_im;
  double q_a_re;
  double q_a_im;
  double r_a_re;
  double s_a_re;
  double r_a_im;
  double s_a_im;
  double t_a_re;
  double t_a_im;
  creal_T a;
  double u_a_re;
  double v_a_re;
  double u_a_im;
  double v_a_im;
  double w_a_re;
  double w_a_im;
  double x_a_re;
  double x_a_im;
  double y_a_re;
  double y_a_im;
  double ab_a_re;
  double ab_a_im;
  double bb_a_re;
  double cb_a_re;
  double bb_a_im;
  double cb_a_im;
  double db_a_re;
  double db_a_im;
  double eb_a_re;
  double eb_a_im;
  double fb_a_re;
  double r;
  creal_T dc1;
  creal_T dc2;
  creal_T dc3;
  double p_re;
  double k_im;
  double q_re;
  double r_re;
  double s_re;
  double l_im;
  double e_m_re;
  double d_m_im;
  double n_im;
  double t_re;
  double f_m_re;
  double f_m_im;
  double u_re;
  double v_re;
  double o_im;
  double w_re;
  double p_im;
  double q_im;
  double x_re;
  double y_re;
  double ab_re;
  double r_im;
  double s_im;
  double bb_re;
  double t_im;
  double cb_re;
  double u_im;
  double v_im;
  double db_re;
  double e_d_re;
  double e_d_im;
  double eb_re;
  double fb_re;
  double w_im;
  double x_im;
  double gb_re;
  creal_T dc4;
  double y_im;
  double hb_re;
  creal_T dc5;
  double ib_re;
  double ab_im;
  double bb_im;
  double jb_re;
  creal_T dc6;
  b_m.re = m;
  b_m.im = 0.0;
  b_R.re = R;
  b_R.im = 0.0;
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                   + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R -
    c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d *
              rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm -
            f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del)) +
                    (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im * b_m_re))
                  + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re *
    0.0 + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) +
             (h_re * 0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re *
            0.0 + f_im * rm)) + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  c_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  c_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  d_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  d_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  e_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  f_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  f_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  g_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  g_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  dc0 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  re = -m;
  d_re = re * R;
  c_im = re * 0.0 + -0.0 * R;
  re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (d_re * rm - c_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (d_re * 0.0 + c_im * rm);
  b_re = re * b_d_re - im * b_d_im;
  im = re * b_d_im + im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (c_re * 0.0 + b_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  h_a_re = f_a_re * f_a_re - f_a_im * f_a_im;
  f_a_im = f_a_re * f_a_im + f_a_im * f_a_re;
  re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  b_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  c_re = -del;
  d_re = c_re * m;
  c_im = c_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_re = 2.0 * d;
  e_re = 2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = -2.0 * (del * del);
  e_im = -2.0 * (del * 0.0 + 0.0 * del);
  g_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -4.0 * d;
  i_re = h_re * c_m_re - -0.0 * c_m_im;
  f_im = h_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 6.0 * del;
  j_re = h_re * c_m_re - 0.0 * c_m_im;
  g_im = h_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  h_re = -4.0 * (m * m);
  h_im = -4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (d_re * R - c_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  c_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (d_re * 0.0 + c_im * R)) + (m *
    R_im + 0.0 * R_re));
  d_re = ((((((((((((c_re * del + -2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * b_R_re - h_im * b_R_im)) + d * rm)
             + k_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((c_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + 0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * b_R_im + h_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (k_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (n_re *
            0.0 + i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  c_re = l_re * d_re - c_im * d_im;
  c_im = l_re * d_im + c_im * d_re;
  d_re = -2.0 * d;
  e_re = -2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = 2.0 * (del * del);
  e_im = 2.0 * (del * 0.0 + 0.0 * del);
  g_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 4.0 * d;
  i_re = h_re * c_m_re - 0.0 * c_m_im;
  f_im = h_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -6.0 * del;
  j_re = h_re * c_m_re - -0.0 * c_m_im;
  g_im = h_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = 4.0 * (m * m);
  h_im = 4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = ((((((((((((d_re * del + 2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * R_re - h_im * R_im)) + d * rm) +
             k_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((d_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + -0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * R_im + h_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (k_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (n_re * 0.0 +
            i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  e_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  d_re = -m;
  f_re = d_re * R;
  f_im = d_re * 0.0 + -0.0 * R;
  i_a_re = g_a_re * g_a_re - g_a_im * g_a_im;
  g_a_im = g_a_re * g_a_im + g_a_im * g_a_re;
  d_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (e_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (f_re * rm - f_im * 0.0));
  e_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (e_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re)) + (f_re * 0.0 + f_im * rm));
  f_a_re = (((-432.0 * dc0.re + (b_re * c_d_re - im * b_d_im)) + (re * h_a_re -
              b_im * f_a_im)) + (c_re * l_re - c_im * d_im)) + (d_re * i_a_re -
    e_im * g_a_im);
  f_a_im = (((-432.0 * dc0.im + (b_re * b_d_im + im * c_d_re)) + (re * f_a_im +
              b_im * h_a_re)) + (c_re * d_im + c_im * l_re)) + (d_re * g_a_im +
    e_im * i_a_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  g_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  g_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  h_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  h_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  i_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  i_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  j_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  j_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  k_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  k_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  dc0 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  re = -m;
  d_re = re * R;
  c_im = re * 0.0 + -0.0 * R;
  re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (d_re * rm - c_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (d_re * 0.0 + c_im * rm);
  b_re = re * b_d_re - im * b_d_im;
  im = re * b_d_im + im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (c_re * 0.0 + b_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  l_a_re = j_a_re * j_a_re - j_a_im * j_a_im;
  j_a_im = j_a_re * j_a_im + j_a_im * j_a_re;
  re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  b_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  c_re = -del;
  d_re = c_re * m;
  c_im = c_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_re = 2.0 * d;
  e_re = 2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = -2.0 * (del * del);
  e_im = -2.0 * (del * 0.0 + 0.0 * del);
  g_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -4.0 * d;
  i_re = h_re * c_m_re - -0.0 * c_m_im;
  f_im = h_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 6.0 * del;
  j_re = h_re * c_m_re - 0.0 * c_m_im;
  g_im = h_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  h_re = -4.0 * (m * m);
  h_im = -4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (d_re * R - c_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  c_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (d_re * 0.0 + c_im * R)) + (m *
    R_im + 0.0 * R_re));
  d_re = ((((((((((((c_re * del + -2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * b_R_re - h_im * b_R_im)) + d * rm)
             + k_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((c_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + 0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * b_R_im + h_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (k_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (n_re *
            0.0 + i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  c_re = l_re * d_re - c_im * d_im;
  c_im = l_re * d_im + c_im * d_re;
  d_re = -2.0 * d;
  e_re = -2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = 2.0 * (del * del);
  e_im = 2.0 * (del * 0.0 + 0.0 * del);
  g_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 4.0 * d;
  i_re = h_re * c_m_re - 0.0 * c_m_im;
  f_im = h_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -6.0 * del;
  j_re = h_re * c_m_re - -0.0 * c_m_im;
  g_im = h_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = 4.0 * (m * m);
  h_im = 4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = ((((((((((((d_re * del + 2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * R_re - h_im * R_im)) + d * rm) +
             k_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((d_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + -0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * R_im + h_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (k_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (n_re * 0.0 +
            i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  e_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  d_re = -m;
  f_re = d_re * R;
  f_im = d_re * 0.0 + -0.0 * R;
  m_a_re = k_a_re * k_a_re - k_a_im * k_a_im;
  k_a_im = k_a_re * k_a_im + k_a_im * k_a_re;
  d_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (e_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (f_re * rm - f_im * 0.0));
  e_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (e_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re)) + (f_re * 0.0 + f_im * rm));
  j_a_re = (((-432.0 * dc0.re + (b_re * c_d_re - im * b_d_im)) + (re * l_a_re -
              b_im * j_a_im)) + (c_re * l_re - c_im * d_im)) + (d_re * m_a_re -
    e_im * k_a_im);
  j_a_im = (((-432.0 * dc0.im + (b_re * b_d_im + im * c_d_re)) + (re * j_a_im +
              b_im * l_a_re)) + (c_re * d_im + c_im * l_re)) + (d_re * k_a_im +
    e_im * m_a_re);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  k_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  k_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  l_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  l_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  m_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  m_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  n_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  n_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  o_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  o_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  p_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  p_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  q_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  q_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  dc0 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  re = -m;
  d_re = re * R;
  c_im = re * 0.0 + -0.0 * R;
  re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (d_re * rm - c_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (d_re * 0.0 + c_im * rm);
  b_re = re * b_d_re - im * b_d_im;
  im = re * b_d_im + im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (c_re * 0.0 + b_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  r_a_re = p_a_re * p_a_re - p_a_im * p_a_im;
  p_a_im = p_a_re * p_a_im + p_a_im * p_a_re;
  re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  b_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  c_re = -del;
  d_re = c_re * m;
  c_im = c_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_re = 2.0 * d;
  e_re = 2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = -2.0 * (del * del);
  e_im = -2.0 * (del * 0.0 + 0.0 * del);
  g_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -4.0 * d;
  i_re = h_re * c_m_re - -0.0 * c_m_im;
  f_im = h_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 6.0 * del;
  j_re = h_re * c_m_re - 0.0 * c_m_im;
  g_im = h_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  h_re = -4.0 * (m * m);
  h_im = -4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (d_re * R - c_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  c_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (d_re * 0.0 + c_im * R)) + (m *
    R_im + 0.0 * R_re));
  d_re = ((((((((((((c_re * del + -2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * b_R_re - h_im * b_R_im)) + d * rm)
             + k_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((c_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + 0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * b_R_im + h_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (k_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (n_re *
            0.0 + i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  c_re = l_re * d_re - c_im * d_im;
  c_im = l_re * d_im + c_im * d_re;
  d_re = -2.0 * d;
  e_re = -2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = 2.0 * (del * del);
  e_im = 2.0 * (del * 0.0 + 0.0 * del);
  g_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 4.0 * d;
  i_re = h_re * c_m_re - 0.0 * c_m_im;
  f_im = h_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -6.0 * del;
  j_re = h_re * c_m_re - -0.0 * c_m_im;
  g_im = h_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = 4.0 * (m * m);
  h_im = 4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = ((((((((((((d_re * del + 2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * R_re - h_im * R_im)) + d * rm) +
             k_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((d_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + -0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * R_im + h_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (k_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (n_re * 0.0 +
            i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  e_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  d_re = -m;
  f_re = d_re * R;
  f_im = d_re * 0.0 + -0.0 * R;
  s_a_re = q_a_re * q_a_re - q_a_im * q_a_im;
  q_a_im = q_a_re * q_a_im + q_a_im * q_a_re;
  d_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (e_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (f_re * rm - f_im * 0.0));
  e_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (e_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re)) + (f_re * 0.0 + f_im * rm));
  p_a_re = (((-432.0 * dc0.re + (b_re * c_d_re - im * b_d_im)) + (re * r_a_re -
              b_im * p_a_im)) + (c_re * l_re - c_im * d_im)) + (d_re * s_a_re -
    e_im * q_a_im);
  p_a_im = (((-432.0 * dc0.im + (b_re * b_d_im + im * c_d_re)) + (re * p_a_im +
              b_im * r_a_re)) + (c_re * d_im + c_im * l_re)) + (d_re * q_a_im +
    e_im * s_a_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  q_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  q_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  r_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  r_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  s_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  s_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  t_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  t_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  a.re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                   + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R -
    c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d *
              rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm -
            f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  a.im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del)) +
                    (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im * b_m_re))
                  + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re *
    0.0 + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) +
             (h_re * 0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re *
            0.0 + f_im * rm)) + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  dc0 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  re = -m;
  d_re = re * R;
  c_im = re * 0.0 + -0.0 * R;
  re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (d_re * rm - c_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (d_re * 0.0 + c_im * rm);
  b_re = re * b_d_re - im * b_d_im;
  im = re * b_d_im + im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (c_re * 0.0 + b_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  u_a_re = t_a_re * t_a_re - t_a_im * t_a_im;
  t_a_im = t_a_re * t_a_im + t_a_im * t_a_re;
  re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  b_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  c_re = -del;
  d_re = c_re * m;
  c_im = c_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_re = 2.0 * d;
  e_re = 2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = -2.0 * (del * del);
  e_im = -2.0 * (del * 0.0 + 0.0 * del);
  g_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -4.0 * d;
  i_re = h_re * c_m_re - -0.0 * c_m_im;
  f_im = h_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 6.0 * del;
  j_re = h_re * c_m_re - 0.0 * c_m_im;
  g_im = h_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  h_re = -4.0 * (m * m);
  h_im = -4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (d_re * R - c_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  c_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (d_re * 0.0 + c_im * R)) + (m *
    R_im + 0.0 * R_re));
  d_re = ((((((((((((c_re * del + -2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * b_R_re - h_im * b_R_im)) + d * rm)
             + k_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((c_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + 0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * b_R_im + h_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (k_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (n_re *
            0.0 + i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  c_re = l_re * d_re - c_im * d_im;
  c_im = l_re * d_im + c_im * d_re;
  d_re = -2.0 * d;
  e_re = -2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = 2.0 * (del * del);
  e_im = 2.0 * (del * 0.0 + 0.0 * del);
  g_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 4.0 * d;
  i_re = h_re * c_m_re - 0.0 * c_m_im;
  f_im = h_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -6.0 * del;
  j_re = h_re * c_m_re - -0.0 * c_m_im;
  g_im = h_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = 4.0 * (m * m);
  h_im = 4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = ((((((((((((d_re * del + 2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * R_re - h_im * R_im)) + d * rm) +
             k_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((d_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + -0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * R_im + h_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (k_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (n_re * 0.0 +
            i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  e_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  d_re = -m;
  f_re = d_re * R;
  f_im = d_re * 0.0 + -0.0 * R;
  v_a_re = a.re * a.re - a.im * a.im;
  u_a_im = a.re * a.im + a.im * a.re;
  d_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (e_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (f_re * rm - f_im * 0.0));
  e_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (e_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re)) + (f_re * 0.0 + f_im * rm));
  t_a_re = (((-432.0 * dc0.re + (b_re * c_d_re - im * b_d_im)) + (re * u_a_re -
              b_im * t_a_im)) + (c_re * l_re - c_im * d_im)) + (d_re * v_a_re -
    e_im * u_a_im);
  t_a_im = (((-432.0 * dc0.im + (b_re * b_d_im + im * c_d_re)) + (re * t_a_im +
              b_im * u_a_re)) + (c_re * d_im + c_im * l_re)) + (d_re * u_a_im +
    e_im * v_a_re);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  a.re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                   + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R -
    c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d *
              rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm -
            f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  a.im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del)) +
                    (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im * b_m_re))
                  + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re *
    0.0 + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) +
             (h_re * 0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re *
            0.0 + f_im * rm)) + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  u_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  u_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  v_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  v_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  w_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  w_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  x_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  x_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  y_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im * m_im))
                     + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
              - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  y_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                      + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
            + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  ab_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                      + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                 * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
    - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  ab_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                       + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
             + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  dc0 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  re = -m;
  d_re = re * R;
  c_im = re * 0.0 + -0.0 * R;
  re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (d_re * rm - c_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (d_re * 0.0 + c_im * rm);
  b_re = re * b_d_re - im * b_d_im;
  im = re * b_d_im + im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (c_re * 0.0 + b_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  bb_a_re = y_a_re * y_a_re - y_a_im * y_a_im;
  y_a_im = y_a_re * y_a_im + y_a_im * y_a_re;
  re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  b_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  c_re = -del;
  d_re = c_re * m;
  c_im = c_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_re = 2.0 * d;
  e_re = 2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = -2.0 * (del * del);
  e_im = -2.0 * (del * 0.0 + 0.0 * del);
  g_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -4.0 * d;
  i_re = h_re * c_m_re - -0.0 * c_m_im;
  f_im = h_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 6.0 * del;
  j_re = h_re * c_m_re - 0.0 * c_m_im;
  g_im = h_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  h_re = -4.0 * (m * m);
  h_im = -4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (d_re * R - c_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  c_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (d_re * 0.0 + c_im * R)) + (m *
    R_im + 0.0 * R_re));
  d_re = ((((((((((((c_re * del + -2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * b_R_re - h_im * b_R_im)) + d * rm)
             + k_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((c_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + 0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * b_R_im + h_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (k_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (n_re *
            0.0 + i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  c_re = l_re * d_re - c_im * d_im;
  c_im = l_re * d_im + c_im * d_re;
  d_re = -2.0 * d;
  e_re = -2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = 2.0 * (del * del);
  e_im = 2.0 * (del * 0.0 + 0.0 * del);
  g_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 4.0 * d;
  i_re = h_re * c_m_re - 0.0 * c_m_im;
  f_im = h_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -6.0 * del;
  j_re = h_re * c_m_re - -0.0 * c_m_im;
  g_im = h_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = 4.0 * (m * m);
  h_im = 4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = ((((((((((((d_re * del + 2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * R_re - h_im * R_im)) + d * rm) +
             k_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((d_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + -0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * R_im + h_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (k_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (n_re * 0.0 +
            i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  e_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  d_re = -m;
  f_re = d_re * R;
  f_im = d_re * 0.0 + -0.0 * R;
  cb_a_re = ab_a_re * ab_a_re - ab_a_im * ab_a_im;
  ab_a_im = ab_a_re * ab_a_im + ab_a_im * ab_a_re;
  d_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (e_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (f_re * rm - f_im * 0.0));
  e_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (e_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re)) + (f_re * 0.0 + f_im * rm));
  y_a_re = (((-432.0 * dc0.re + (b_re * c_d_re - im * b_d_im)) + (re * bb_a_re -
              b_im * y_a_im)) + (c_re * l_re - c_im * d_im)) + (d_re * cb_a_re -
    e_im * ab_a_im);
  y_a_im = (((-432.0 * dc0.im + (b_re * b_d_im + im * c_d_re)) + (re * y_a_im +
              b_im * bb_a_re)) + (c_re * d_im + c_im * l_re)) + (d_re * ab_a_im
    + e_im * cb_a_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  ab_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im *
    m_im)) + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R - c_im *
    0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d * rm) +
                h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm - f_im *
    0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  ab_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                       + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
             + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  bb_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                      + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                 * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
    - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  bb_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                       + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
             + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  cb_a_re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  cb_a_im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  re = 2.0 * d;
  b_re = 2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = -2.0 * (del * del);
  b_im = -2.0 * (del * 0.0 + 0.0 * del);
  d_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -4.0 * d;
  f_re = e_re * c_m_re - -0.0 * c_m_im;
  c_im = e_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 6.0 * del;
  g_re = e_re * c_m_re - 0.0 * c_m_im;
  d_im = e_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -4.0 * (m * m);
  e_im = -4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  db_a_re = ((((((((((((re * del + -2.0 * (del * del)) + (c_re * m_re - im *
    m_im)) + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R - c_im *
    0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d * rm) +
                h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm - f_im *
    0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  db_a_im = (((((((((((((re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                       + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + 0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
             + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  re = -2.0 * d;
  b_re = -2.0 * d;
  c_re = b_re * del;
  im = b_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  b_re = 2.0 * (del * del);
  b_im = 2.0 * (del * 0.0 + 0.0 * del);
  d_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = 4.0 * d;
  f_re = e_re * c_m_re - 0.0 * c_m_im;
  c_im = e_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_re = -6.0 * del;
  g_re = e_re * c_m_re - -0.0 * c_m_im;
  d_im = e_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = 4.0 * (m * m);
  e_im = 4.0 * (m * 0.0 + 0.0 * m);
  h_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * del;
  j_re = i_re * c_m_re - -0.0 * c_m_im;
  f_im = i_re * c_m_im + -0.0 * c_m_re;
  i_re = 3.0 * (m * m);
  g_im = 3.0 * (m * 0.0 + 0.0 * m);
  k_re = i_re * R - g_im * 0.0;
  g_im = i_re * 0.0 + g_im * R;
  eb_a_re = ((((((((((((re * del + 2.0 * (del * del)) + (c_re * m_re - im * m_im))
                      + (b_re * b_m_re - b_im * b_m_im)) + d_re * R) + (f_re * R
    - c_im * 0.0)) + (g_re * R - d_im * 0.0)) + (e_re * R_re - e_im * R_im)) + d
                 * rm) + h_re * rm) + (b_d_re * rm - b_d_im * 0.0)) + (j_re * rm
    - f_im * 0.0)) + R * rm) + (k_re * rm - g_im * 0.0);
  eb_a_im = (((((((((((((re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                       + (c_re * m_im + im * m_re)) + (b_re * b_m_im + b_im *
    b_m_re)) + (d_re * 0.0 + -0.0 * R)) + (f_re * 0.0 + c_im * R)) + (g_re * 0.0
    + d_im * R)) + (e_re * R_im + e_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (h_re *
    0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (j_re * 0.0 + f_im * rm))
             + (R * 0.0 + 0.0 * rm)) + (k_re * 0.0 + g_im * rm);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m * R_im +
    0.0 * R_re);
  dc0 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  b_re = re * m;
  im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  re = -m;
  d_re = re * R;
  c_im = re * 0.0 + -0.0 * R;
  re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (b_re * R - im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (b_re * 0.0 + im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (d_re * rm - c_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (d_re * 0.0 + c_im * rm);
  b_re = re * b_d_re - im * b_d_im;
  im = re * b_d_im + im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (c_re * 0.0 + b_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  re = -del;
  c_re = re * m;
  b_im = re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  fb_a_re = db_a_re * db_a_re - db_a_im * db_a_im;
  db_a_im = db_a_re * db_a_im + db_a_im * db_a_re;
  re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (c_re * R - b_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  b_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (c_re * 0.0 + b_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  c_re = -del;
  d_re = c_re * m;
  c_im = c_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_re = 2.0 * d;
  e_re = 2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = -2.0 * (del * del);
  e_im = -2.0 * (del * 0.0 + 0.0 * del);
  g_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -4.0 * d;
  i_re = h_re * c_m_re - -0.0 * c_m_im;
  f_im = h_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 6.0 * del;
  j_re = h_re * c_m_re - 0.0 * c_m_im;
  g_im = h_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  h_re = -4.0 * (m * m);
  h_im = -4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (d_re * R - c_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  c_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (d_re * 0.0 + c_im * R)) + (m *
    R_im + 0.0 * R_re));
  d_re = ((((((((((((c_re * del + -2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * b_R_re - h_im * b_R_im)) + d * rm)
             + k_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((c_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + 0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * b_R_im + h_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (k_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (n_re *
            0.0 + i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  c_re = l_re * d_re - c_im * d_im;
  c_im = l_re * d_im + c_im * d_re;
  d_re = -2.0 * d;
  e_re = -2.0 * d;
  f_re = e_re * del;
  d_im = e_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  e_re = 2.0 * (del * del);
  e_im = 2.0 * (del * 0.0 + 0.0 * del);
  g_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = 4.0 * d;
  i_re = h_re * c_m_re - 0.0 * c_m_im;
  f_im = h_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  h_re = -6.0 * del;
  j_re = h_re * c_m_re - -0.0 * c_m_im;
  g_im = h_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = 4.0 * (m * m);
  h_im = 4.0 * (m * 0.0 + 0.0 * m);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * del;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  l_re = 3.0 * (m * m);
  j_im = 3.0 * (m * 0.0 + 0.0 * m);
  o_re = l_re * R - j_im * 0.0;
  j_im = l_re * 0.0 + j_im * R;
  l_re = ((((((((((((d_re * del + 2.0 * (del * del)) + (f_re * m_re - d_im *
    m_im)) + (e_re * b_m_re - e_im * b_m_im)) + g_re * R) + (i_re * R - f_im *
    0.0)) + (j_re * R - g_im * 0.0)) + (h_re * R_re - h_im * R_im)) + d * rm) +
             k_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (n_re * rm - i_im *
            0.0)) + R * rm) + (o_re * rm - j_im * 0.0);
  d_im = (((((((((((((d_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (f_re * m_im + d_im * m_re)) + (e_re * b_m_im + e_im *
    b_m_re)) + (g_re * 0.0 + -0.0 * R)) + (i_re * 0.0 + f_im * R)) + (j_re * 0.0
    + g_im * R)) + (h_re * R_im + h_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (k_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (n_re * 0.0 +
            i_im * rm)) + (R * 0.0 + 0.0 * rm)) + (o_re * 0.0 + j_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  e_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  d_re = -m;
  f_re = d_re * R;
  f_im = d_re * 0.0 + -0.0 * R;
  r = eb_a_re * eb_a_re - eb_a_im * eb_a_im;
  eb_a_im = eb_a_re * eb_a_im + eb_a_im * eb_a_re;
  d_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (e_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (f_re * rm - f_im * 0.0));
  e_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (e_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re)) + (f_re * 0.0 + f_im * rm));
  db_a_re = (((-432.0 * dc0.re + (b_re * c_d_re - im * b_d_im)) + (re * fb_a_re
    - b_im * db_a_im)) + (c_re * l_re - c_im * d_im)) + (d_re * r - e_im *
    eb_a_im);
  db_a_im = (((-432.0 * dc0.im + (b_re * b_d_im + im * c_d_re)) + (re * db_a_im
    + b_im * fb_a_re)) + (c_re * d_im + c_im * l_re)) + (d_re * eb_a_im + e_im *
    r);
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc0 = power(b_d);
  re = 6.0 * ((d + -del) + R);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  b_re = -del;
  c_re = b_re * m;
  im = b_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  b_re = -del;
  d_re = b_re * m;
  b_im = b_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  b_d.re = (((c_d_re * R - c_d_im * 0.0) + (d_re * R - b_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (m_re * rm - m_im * 0.0);
  b_d.im = (((c_d_re * 0.0 + c_d_im * R) + (d_re * 0.0 + b_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (m_re * 0.0 + m_im * rm);
  dc1 = power(b_d);
  b_re = -2.0 * (((b_d_re * R - b_d_im * 0.0) + (c_re * R - im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  im = -2.0 * (((b_d_re * 0.0 + b_d_im * R) + (c_re * 0.0 + im * R)) + (m * R_im
    + 0.0 * R_re));
  dc2 = b_power(b_m);
  dc3 = b_power(b_R);
  c_re = 0.25 * dc2.re;
  b_im = 0.25 * dc2.im;
  d_re = c_re * dc3.re - b_im * dc3.im;
  b_im = c_re * dc3.im + b_im * dc3.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc2 = b_power(b_d);
  c_re = d_re * dc2.re - b_im * dc2.im;
  b_im = d_re * dc2.im + b_im * dc2.re;
  eb_a_re = a.re * a.re - a.im * a.im;
  eb_a_im = a.re * a.im + a.im * a.re;
  dc2 = power(b_m);
  dc3 = power(b_R);
  d_re = 0.41997368329829105 * dc2.re;
  c_im = 0.41997368329829105 * dc2.im;
  e_re = d_re * dc3.re - c_im * dc3.im;
  c_im = d_re * dc3.im + c_im * dc3.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc2 = power(b_d);
  d_re = e_re * dc2.re - c_im * dc2.im;
  c_im = e_re * dc2.im + c_im * dc2.re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  e_re = -del;
  f_re = e_re * m;
  d_im = e_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  e_re = -m;
  g_re = e_re * R;
  e_im = e_re * 0.0 + -0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  e_re = -del;
  h_re = e_re * m;
  f_im = e_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  e_re = 12.0 * ((((b_d_re * R - b_d_im * 0.0) + (f_re * R - d_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (g_re * rm - e_im * 0.0));
  d_im = 12.0 * ((((b_d_re * 0.0 + b_d_im * R) + (f_re * 0.0 + d_im * R)) + (m *
    R_im + 0.0 * R_re)) + (g_re * 0.0 + e_im * rm));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (h_re * R - f_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (h_re * 0.0 + f_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (m_re * 0.0 + m_im * rm);
  f_re = 2.0 * d;
  g_re = 2.0 * d;
  h_re = g_re * del;
  e_im = g_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  g_re = -2.0 * (del * del);
  f_im = -2.0 * (del * 0.0 + 0.0 * del);
  i_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  j_re = -4.0 * d;
  k_re = j_re * c_m_re - -0.0 * c_m_im;
  g_im = j_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  j_re = 6.0 * del;
  l_re = j_re * c_m_re - 0.0 * c_m_im;
  h_im = j_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  j_re = -4.0 * (m * m);
  i_im = -4.0 * (m * 0.0 + 0.0 * m);
  n_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  o_re = -2.0 * del;
  p_re = o_re * c_m_re - -0.0 * c_m_im;
  j_im = o_re * c_m_im + -0.0 * c_m_re;
  o_re = 3.0 * (m * m);
  k_im = 3.0 * (m * 0.0 + 0.0 * m);
  q_re = o_re * R - k_im * 0.0;
  k_im = o_re * 0.0 + k_im * R;
  o_re = -2.0 * d;
  r_re = -2.0 * d;
  s_re = r_re * del;
  l_im = r_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  r_re = 2.0 * (del * del);
  n_im = 2.0 * (del * 0.0 + 0.0 * del);
  t_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  u_re = 4.0 * d;
  v_re = u_re * f_m_re - 0.0 * f_m_im;
  o_im = u_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  u_re = -6.0 * del;
  w_re = u_re * f_m_re - -0.0 * f_m_im;
  p_im = u_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  u_re = 4.0 * (m * m);
  q_im = 4.0 * (m * 0.0 + 0.0 * m);
  x_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * f_m_re - 0.0 * f_m_im;
  d_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = -2.0 * del;
  ab_re = y_re * f_m_re - -0.0 * f_m_im;
  r_im = y_re * f_m_im + -0.0 * f_m_re;
  y_re = 3.0 * (m * m);
  s_im = 3.0 * (m * 0.0 + 0.0 * m);
  bb_re = y_re * R - s_im * 0.0;
  s_im = y_re * 0.0 + s_im * R;
  y_re = -3.0 * (((((((((((((f_re * del + -2.0 * (del * del)) + (h_re * m_re -
    e_im * m_im)) + (g_re * b_m_re - f_im * b_m_im)) + i_re * R) + (k_re * R -
    g_im * 0.0)) + (l_re * R - h_im * 0.0)) + (j_re * R_re - i_im * R_im)) + d *
                      rm) + n_re * rm) + (c_d_re * rm - c_d_im * 0.0)) + (p_re *
    rm - j_im * 0.0)) + R * rm) + (q_re * rm - k_im * 0.0));
  e_im = -3.0 * ((((((((((((((f_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (h_re * m_im + e_im * m_re)) + (g_re * b_m_im + f_im * b_m_re)) +
    (i_re * 0.0 + 0.0 * R)) + (k_re * 0.0 + g_im * R)) + (l_re * 0.0 + h_im * R))
                       + (j_re * R_im + i_im * R_re)) + (d * 0.0 + 0.0 * rm)) +
                     (n_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + c_d_im * rm)) +
                   (p_re * 0.0 + j_im * rm)) + (R * 0.0 + 0.0 * rm)) + (q_re *
    0.0 + k_im * rm));
  f_re = ((((((((((((o_re * del + 2.0 * (del * del)) + (s_re * c_m_re - l_im *
    c_m_im)) + (r_re * e_m_re - n_im * d_m_im)) + t_re * R) + (v_re * R - o_im *
    0.0)) + (w_re * R - p_im * 0.0)) + (u_re * b_R_re - q_im * b_R_im)) + d * rm)
             + x_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (ab_re * rm - r_im *
            0.0)) + R * rm) + (bb_re * rm - s_im * 0.0);
  f_im = (((((((((((((o_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (s_re * c_m_im + l_im * c_m_re)) + (r_re * d_m_im + n_im *
    e_m_re)) + (t_re * 0.0 + -0.0 * R)) + (v_re * 0.0 + o_im * R)) + (w_re * 0.0
    + p_im * R)) + (u_re * b_R_im + q_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (x_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (ab_re *
            0.0 + r_im * rm)) + (R * 0.0 + 0.0 * rm)) + (bb_re * 0.0 + s_im * rm);
  g_re = (36.0 * (u_a_re * u_a_re - u_a_im * u_a_im) + (e_re * b_d_re - d_im *
           b_d_im)) + (y_re * f_re - e_im * f_im);
  d_im = (36.0 * (u_a_re * u_a_im + u_a_im * u_a_re) + (e_re * b_d_im + d_im *
           b_d_re)) + (y_re * f_im + e_im * f_re);
  e_re = d_re * g_re - c_im * d_im;
  c_im = d_re * d_im + c_im * g_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  f_re = d_re * m;
  d_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (f_re * R - d_im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (f_re * 0.0 + d_im * R)) + (m * R_im +
    0.0 * R_re);
  dc2 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  f_re = d_re * m;
  d_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  g_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  d_re = -m;
  h_re = d_re * R;
  f_im = d_re * 0.0 + -0.0 * R;
  d_re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (f_re * R - d_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  d_im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (f_re * 0.0 + d_im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (g_re * R - e_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (h_re * rm - f_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (g_re * 0.0 + e_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (h_re * 0.0 + f_im * rm);
  f_re = d_re * b_d_re - d_im * b_d_im;
  d_im = d_re * b_d_im + d_im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  g_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (g_re * R - e_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (g_re * 0.0 + e_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  d_re = -del;
  g_re = d_re * m;
  e_im = d_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  u_a_re = v_a_re * v_a_re - v_a_im * v_a_im;
  u_a_im = v_a_re * v_a_im + v_a_im * v_a_re;
  d_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (g_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  e_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (g_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  g_re = -del;
  h_re = g_re * m;
  f_im = g_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  g_re = 2.0 * d;
  i_re = 2.0 * d;
  j_re = i_re * del;
  g_im = i_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  i_re = -2.0 * (del * del);
  h_im = -2.0 * (del * 0.0 + 0.0 * del);
  k_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -4.0 * d;
  n_re = l_re * c_m_re - -0.0 * c_m_im;
  i_im = l_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = 6.0 * del;
  o_re = l_re * c_m_re - 0.0 * c_m_im;
  j_im = l_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  l_re = -4.0 * (m * m);
  k_im = -4.0 * (m * 0.0 + 0.0 * m);
  p_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  q_re = -2.0 * del;
  r_re = q_re * c_m_re - -0.0 * c_m_im;
  l_im = q_re * c_m_im + -0.0 * c_m_re;
  q_re = 3.0 * (m * m);
  n_im = 3.0 * (m * 0.0 + 0.0 * m);
  s_re = q_re * R - n_im * 0.0;
  n_im = q_re * 0.0 + n_im * R;
  q_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (h_re * R - f_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  f_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (h_re * 0.0 + f_im * R)) + (m *
    R_im + 0.0 * R_re));
  h_re = ((((((((((((g_re * del + -2.0 * (del * del)) + (j_re * m_re - g_im *
    m_im)) + (i_re * b_m_re - h_im * b_m_im)) + k_re * R) + (n_re * R - i_im *
    0.0)) + (o_re * R - j_im * 0.0)) + (l_re * b_R_re - k_im * b_R_im)) + d * rm)
             + p_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (r_re * rm - l_im *
            0.0)) + R * rm) + (s_re * rm - n_im * 0.0);
  g_im = (((((((((((((g_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (j_re * m_im + g_im * m_re)) + (i_re * b_m_im + h_im *
    b_m_re)) + (k_re * 0.0 + 0.0 * R)) + (n_re * 0.0 + i_im * R)) + (o_re * 0.0
    + j_im * R)) + (l_re * b_R_im + k_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (p_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (r_re *
            0.0 + l_im * rm)) + (R * 0.0 + 0.0 * rm)) + (s_re * 0.0 + n_im * rm);
  g_re = q_re * h_re - f_im * g_im;
  f_im = q_re * g_im + f_im * h_re;
  h_re = -2.0 * d;
  i_re = -2.0 * d;
  j_re = i_re * del;
  g_im = i_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  i_re = 2.0 * (del * del);
  h_im = 2.0 * (del * 0.0 + 0.0 * del);
  k_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = 4.0 * d;
  n_re = l_re * c_m_re - 0.0 * c_m_im;
  i_im = l_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  l_re = -6.0 * del;
  o_re = l_re * c_m_re - -0.0 * c_m_im;
  j_im = l_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  l_re = 4.0 * (m * m);
  k_im = 4.0 * (m * 0.0 + 0.0 * m);
  p_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  q_re = -2.0 * del;
  r_re = q_re * c_m_re - -0.0 * c_m_im;
  l_im = q_re * c_m_im + -0.0 * c_m_re;
  q_re = 3.0 * (m * m);
  n_im = 3.0 * (m * 0.0 + 0.0 * m);
  s_re = q_re * R - n_im * 0.0;
  n_im = q_re * 0.0 + n_im * R;
  q_re = ((((((((((((h_re * del + 2.0 * (del * del)) + (j_re * m_re - g_im *
    m_im)) + (i_re * b_m_re - h_im * b_m_im)) + k_re * R) + (n_re * R - i_im *
    0.0)) + (o_re * R - j_im * 0.0)) + (l_re * R_re - k_im * R_im)) + d * rm) +
             p_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (r_re * rm - l_im *
            0.0)) + R * rm) + (s_re * rm - n_im * 0.0);
  g_im = (((((((((((((h_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (j_re * m_im + g_im * m_re)) + (i_re * b_m_im + h_im *
    b_m_re)) + (k_re * 0.0 + -0.0 * R)) + (n_re * 0.0 + i_im * R)) + (o_re * 0.0
    + j_im * R)) + (l_re * R_im + k_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (p_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (r_re * 0.0 +
            l_im * rm)) + (R * 0.0 + 0.0 * rm)) + (s_re * 0.0 + n_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  h_re = -del;
  i_re = h_re * m;
  h_im = h_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = -m;
  j_re = h_re * R;
  i_im = h_re * 0.0 + -0.0 * R;
  v_a_re = w_a_re * w_a_re - w_a_im * w_a_im;
  v_a_im = w_a_re * w_a_im + w_a_im * w_a_re;
  h_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (i_re * R - h_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (j_re * rm - i_im * 0.0));
  h_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (i_re * 0.0 + h_im * R)) + (m *
    R_im + 0.0 * R_re)) + (j_re * 0.0 + i_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  i_re = -del;
  j_re = i_re * m;
  i_im = i_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  i_re = -m;
  k_re = i_re * R;
  j_im = i_re * 0.0 + -0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  i_re = -del;
  l_re = i_re * m;
  k_im = i_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  i_re = 12.0 * ((((b_d_re * R - c_d_im * 0.0) + (j_re * R - i_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (k_re * rm - j_im * 0.0));
  i_im = 12.0 * ((((b_d_re * 0.0 + c_d_im * R) + (j_re * 0.0 + i_im * R)) + (m *
    R_im + 0.0 * R_re)) + (k_re * 0.0 + j_im * rm));
  b_d_re = (((d_d_re * R - d_d_im * 0.0) + (l_re * R - k_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (m_re * rm - m_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (l_re * 0.0 + k_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (m_re * 0.0 + m_im * rm);
  j_re = 2.0 * d;
  k_re = 2.0 * d;
  l_re = k_re * del;
  j_im = k_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  k_re = -2.0 * (del * del);
  k_im = -2.0 * (del * 0.0 + 0.0 * del);
  n_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  o_re = -4.0 * d;
  p_re = o_re * c_m_re - -0.0 * c_m_im;
  l_im = o_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  o_re = 6.0 * del;
  r_re = o_re * c_m_re - 0.0 * c_m_im;
  n_im = o_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  o_re = -4.0 * (m * m);
  o_im = -4.0 * (m * 0.0 + 0.0 * m);
  s_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  t_re = -2.0 * del;
  u_re = t_re * c_m_re - -0.0 * c_m_im;
  p_im = t_re * c_m_im + -0.0 * c_m_re;
  t_re = 3.0 * (m * m);
  q_im = 3.0 * (m * 0.0 + 0.0 * m);
  v_re = t_re * R - q_im * 0.0;
  q_im = t_re * 0.0 + q_im * R;
  t_re = -2.0 * d;
  w_re = -2.0 * d;
  x_re = w_re * del;
  r_im = w_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  w_re = 2.0 * (del * del);
  s_im = 2.0 * (del * 0.0 + 0.0 * del);
  y_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = 4.0 * d;
  bb_re = ab_re * f_m_re - 0.0 * f_m_im;
  t_im = ab_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = -6.0 * del;
  cb_re = ab_re * f_m_re - -0.0 * f_m_im;
  u_im = ab_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  ab_re = 4.0 * (m * m);
  v_im = 4.0 * (m * 0.0 + 0.0 * m);
  db_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  eb_re = -2.0 * del;
  fb_re = eb_re * f_m_re - -0.0 * f_m_im;
  w_im = eb_re * f_m_im + -0.0 * f_m_re;
  eb_re = 3.0 * (m * m);
  x_im = 3.0 * (m * 0.0 + 0.0 * m);
  gb_re = eb_re * R - x_im * 0.0;
  x_im = eb_re * 0.0 + x_im * R;
  eb_re = -3.0 * (((((((((((((j_re * del + -2.0 * (del * del)) + (l_re * m_re -
    j_im * m_im)) + (k_re * b_m_re - k_im * b_m_im)) + n_re * R) + (p_re * R -
    l_im * 0.0)) + (r_re * R - n_im * 0.0)) + (o_re * R_re - o_im * R_im)) + d *
                       rm) + s_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (u_re *
    rm - p_im * 0.0)) + R * rm) + (v_re * rm - q_im * 0.0));
  j_im = -3.0 * ((((((((((((((j_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (l_re * m_im + j_im * m_re)) + (k_re * b_m_im + k_im * b_m_re)) +
    (n_re * 0.0 + 0.0 * R)) + (p_re * 0.0 + l_im * R)) + (r_re * 0.0 + n_im * R))
                       + (o_re * R_im + o_im * R_re)) + (d * 0.0 + 0.0 * rm)) +
                     (s_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) +
                   (u_re * 0.0 + p_im * rm)) + (R * 0.0 + 0.0 * rm)) + (v_re *
    0.0 + q_im * rm));
  j_re = ((((((((((((t_re * del + 2.0 * (del * del)) + (x_re * c_m_re - r_im *
    c_m_im)) + (w_re * e_m_re - s_im * d_m_im)) + y_re * R) + (bb_re * R - t_im *
    0.0)) + (cb_re * R - u_im * 0.0)) + (ab_re * b_R_re - v_im * b_R_im)) + d *
              rm) + db_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (fb_re * rm -
            w_im * 0.0)) + R * rm) + (gb_re * rm - x_im * 0.0);
  k_im = (((((((((((((t_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (x_re * c_m_im + r_im * c_m_re)) + (w_re * d_m_im + s_im *
    e_m_re)) + (y_re * 0.0 + -0.0 * R)) + (bb_re * 0.0 + t_im * R)) + (cb_re *
    0.0 + u_im * R)) + (ab_re * b_R_im + v_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
             + (db_re * 0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) +
           (fb_re * 0.0 + w_im * rm)) + (R * 0.0 + 0.0 * rm)) + (gb_re * 0.0 +
    x_im * rm);
  dc3.re = (36.0 * (x_a_re * x_a_re - x_a_im * x_a_im) + (i_re * b_d_re - i_im *
             c_d_im)) + (eb_re * j_re - j_im * k_im);
  dc3.im = (36.0 * (x_a_re * x_a_im + x_a_im * x_a_re) + (i_re * c_d_im + i_im *
             b_d_re)) + (eb_re * k_im + j_im * j_re);
  dc3 = c_power(dc3);
  dc4.re = -4.0 * dc3.re + (y_a_re * y_a_re - y_a_im * y_a_im);
  dc4.im = -4.0 * dc3.im + (y_a_re * y_a_im + y_a_im * y_a_re);
  dc3 = d_power(dc4);
  dc4.re = ((((-432.0 * dc2.re + (f_re * c_d_re - d_im * b_d_im)) + (d_re *
    u_a_re - e_im * u_a_im)) + (g_re * q_re - f_im * g_im)) + (h_re * v_a_re -
             h_im * v_a_im)) + dc3.re;
  dc4.im = ((((-432.0 * dc2.im + (f_re * b_d_im + d_im * c_d_re)) + (d_re *
    u_a_im + e_im * u_a_re)) + (g_re * g_im + f_im * q_re)) + (h_re * v_a_im +
             h_im * v_a_re)) + dc3.im;
  dc2 = e_power(dc4);
  dc3 = power(b_m);
  dc4 = power(b_R);
  d_re = 0.26456684199469993 * dc3.re;
  d_im = 0.26456684199469993 * dc3.im;
  f_re = d_re * dc4.re - d_im * dc4.im;
  d_im = d_re * dc4.im + d_im * dc4.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc3 = power(b_d);
  d_re = f_re * dc3.re - d_im * dc3.im;
  d_im = f_re * dc3.im + d_im * dc3.re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  f_re = -del;
  g_re = f_re * m;
  e_im = f_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((b_d_re * R - b_d_im * 0.0) + (g_re * R - e_im * 0.0)) + (m * R_re -
    0.0 * R_im);
  b_d.im = ((b_d_re * 0.0 + b_d_im * R) + (g_re * 0.0 + e_im * R)) + (m * R_im +
    0.0 * R_re);
  dc3 = c_power(b_d);
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  f_re = -del;
  g_re = f_re * m;
  e_im = f_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  f_re = -del;
  h_re = f_re * m;
  f_im = f_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  f_re = -m;
  i_re = f_re * R;
  g_im = f_re * 0.0 + -0.0 * R;
  f_re = 432.0 * (((b_d_re * R - b_d_im * 0.0) + (g_re * R - e_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  e_im = 432.0 * (((b_d_re * 0.0 + b_d_im * R) + (g_re * 0.0 + e_im * R)) + (m *
    R_im + 0.0 * R_re));
  b_d_re = (((c_d_re * R - c_d_im * 0.0) + (h_re * R - f_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (i_re * rm - g_im * 0.0);
  b_d_im = (((c_d_re * 0.0 + c_d_im * R) + (h_re * 0.0 + f_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (i_re * 0.0 + g_im * rm);
  g_re = f_re * b_d_re - e_im * b_d_im;
  e_im = f_re * b_d_im + e_im * b_d_re;
  b_d_re = d * m;
  b_d_im = d * 0.0 + 0.0 * m;
  f_re = -del;
  h_re = f_re * m;
  f_im = f_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  c_d_re = (((b_d_re * R - b_d_im * 0.0) + (h_re * R - f_im * 0.0)) + (m * R_re
             - 0.0 * R_im)) + (m_re * rm - m_im * 0.0);
  b_d_im = (((b_d_re * 0.0 + b_d_im * R) + (h_re * 0.0 + f_im * R)) + (m * R_im
             + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  f_re = -del;
  h_re = f_re * m;
  f_im = f_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  u_a_re = ab_a_re * ab_a_re - ab_a_im * ab_a_im;
  u_a_im = ab_a_re * ab_a_im + ab_a_im * ab_a_re;
  f_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (h_re * R - f_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (m_re * rm - m_im * 0.0));
  f_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (h_re * 0.0 + f_im * R)) + (m *
    R_im + 0.0 * R_re)) + (m_re * 0.0 + m_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  h_re = -del;
  i_re = h_re * m;
  g_im = h_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  h_re = 2.0 * d;
  j_re = 2.0 * d;
  k_re = j_re * del;
  h_im = j_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  j_re = -2.0 * (del * del);
  i_im = -2.0 * (del * 0.0 + 0.0 * del);
  l_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  n_re = -4.0 * d;
  o_re = n_re * c_m_re - -0.0 * c_m_im;
  j_im = n_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  n_re = 6.0 * del;
  p_re = n_re * c_m_re - 0.0 * c_m_im;
  k_im = n_re * c_m_im + 0.0 * c_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  n_re = -4.0 * (m * m);
  l_im = -4.0 * (m * 0.0 + 0.0 * m);
  q_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  r_re = -2.0 * del;
  s_re = r_re * c_m_re - -0.0 * c_m_im;
  n_im = r_re * c_m_im + -0.0 * c_m_re;
  r_re = 3.0 * (m * m);
  o_im = 3.0 * (m * 0.0 + 0.0 * m);
  t_re = r_re * R - o_im * 0.0;
  o_im = r_re * 0.0 + o_im * R;
  r_re = 54.0 * (((b_d_re * R - c_d_im * 0.0) + (i_re * R - g_im * 0.0)) + (m *
    R_re - 0.0 * R_im));
  g_im = 54.0 * (((b_d_re * 0.0 + c_d_im * R) + (i_re * 0.0 + g_im * R)) + (m *
    R_im + 0.0 * R_re));
  i_re = ((((((((((((h_re * del + -2.0 * (del * del)) + (k_re * m_re - h_im *
    m_im)) + (j_re * b_m_re - i_im * b_m_im)) + l_re * R) + (o_re * R - j_im *
    0.0)) + (p_re * R - k_im * 0.0)) + (n_re * b_R_re - l_im * b_R_im)) + d * rm)
             + q_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (s_re * rm - n_im *
            0.0)) + R * rm) + (t_re * rm - o_im * 0.0);
  h_im = (((((((((((((h_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (k_re * m_im + h_im * m_re)) + (j_re * b_m_im + i_im *
    b_m_re)) + (l_re * 0.0 + 0.0 * R)) + (o_re * 0.0 + j_im * R)) + (p_re * 0.0
    + k_im * R)) + (n_re * b_R_im + l_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (q_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (s_re *
            0.0 + n_im * rm)) + (R * 0.0 + 0.0 * rm)) + (t_re * 0.0 + o_im * rm);
  h_re = r_re * i_re - g_im * h_im;
  g_im = r_re * h_im + g_im * i_re;
  i_re = -2.0 * d;
  j_re = -2.0 * d;
  k_re = j_re * del;
  h_im = j_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  j_re = 2.0 * (del * del);
  i_im = 2.0 * (del * 0.0 + 0.0 * del);
  l_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  n_re = 4.0 * d;
  o_re = n_re * c_m_re - 0.0 * c_m_im;
  j_im = n_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  n_re = -6.0 * del;
  p_re = n_re * c_m_re - -0.0 * c_m_im;
  k_im = n_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  n_re = 4.0 * (m * m);
  l_im = 4.0 * (m * 0.0 + 0.0 * m);
  q_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  c_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  r_re = -2.0 * del;
  s_re = r_re * c_m_re - -0.0 * c_m_im;
  n_im = r_re * c_m_im + -0.0 * c_m_re;
  r_re = 3.0 * (m * m);
  o_im = 3.0 * (m * 0.0 + 0.0 * m);
  t_re = r_re * R - o_im * 0.0;
  o_im = r_re * 0.0 + o_im * R;
  r_re = ((((((((((((i_re * del + 2.0 * (del * del)) + (k_re * m_re - h_im *
    m_im)) + (j_re * b_m_re - i_im * b_m_im)) + l_re * R) + (o_re * R - j_im *
    0.0)) + (p_re * R - k_im * 0.0)) + (n_re * R_re - l_im * R_im)) + d * rm) +
             q_re * rm) + (b_d_re * rm - c_d_im * 0.0)) + (s_re * rm - n_im *
            0.0)) + R * rm) + (t_re * rm - o_im * 0.0);
  h_im = (((((((((((((i_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (k_re * m_im + h_im * m_re)) + (j_re * b_m_im + i_im *
    b_m_re)) + (l_re * 0.0 + -0.0 * R)) + (o_re * 0.0 + j_im * R)) + (p_re * 0.0
    + k_im * R)) + (n_re * R_im + l_im * R_re)) + (d * 0.0 + 0.0 * rm)) + (q_re *
              0.0 + -0.0 * rm)) + (b_d_re * 0.0 + c_d_im * rm)) + (s_re * 0.0 +
            n_im * rm)) + (R * 0.0 + 0.0 * rm)) + (t_re * 0.0 + o_im * rm);
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  i_re = -del;
  j_re = i_re * m;
  i_im = i_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  i_re = -m;
  k_re = i_re * R;
  j_im = i_re * 0.0 + -0.0 * R;
  v_a_re = bb_a_re * bb_a_re - bb_a_im * bb_a_im;
  v_a_im = bb_a_re * bb_a_im + bb_a_im * bb_a_re;
  i_re = 27.0 * ((((b_d_re * R - c_d_im * 0.0) + (j_re * R - i_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (k_re * rm - j_im * 0.0));
  i_im = 27.0 * ((((b_d_re * 0.0 + c_d_im * R) + (j_re * 0.0 + i_im * R)) + (m *
    R_im + 0.0 * R_re)) + (k_re * 0.0 + j_im * rm));
  b_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  j_re = -del;
  k_re = j_re * m;
  j_im = j_re * 0.0 + -0.0 * m;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  j_re = -m;
  l_re = j_re * R;
  k_im = j_re * 0.0 + -0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  j_re = -del;
  n_re = j_re * m;
  l_im = j_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  m_re = m * R;
  m_im = m * 0.0 + 0.0 * R;
  j_re = 12.0 * ((((b_d_re * R - c_d_im * 0.0) + (k_re * R - j_im * 0.0)) + (m *
    R_re - 0.0 * R_im)) + (l_re * rm - k_im * 0.0));
  j_im = 12.0 * ((((b_d_re * 0.0 + c_d_im * R) + (k_re * 0.0 + j_im * R)) + (m *
    R_im + 0.0 * R_re)) + (l_re * 0.0 + k_im * rm));
  b_d_re = (((d_d_re * R - d_d_im * 0.0) + (n_re * R - l_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (m_re * rm - m_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (n_re * 0.0 + l_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (m_re * 0.0 + m_im * rm);
  k_re = 2.0 * d;
  l_re = 2.0 * d;
  n_re = l_re * del;
  k_im = l_re * 0.0 + 0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  l_re = -2.0 * (del * del);
  l_im = -2.0 * (del * 0.0 + 0.0 * del);
  o_re = 2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  p_re = -4.0 * d;
  q_re = p_re * c_m_re - -0.0 * c_m_im;
  n_im = p_re * c_m_im + -0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  p_re = 6.0 * del;
  s_re = p_re * c_m_re - 0.0 * c_m_im;
  o_im = p_re * c_m_im + 0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  p_re = -4.0 * (m * m);
  p_im = -4.0 * (m * 0.0 + 0.0 * m);
  t_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * c_m_re - 0.0 * c_m_im;
  d_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  u_re = -2.0 * del;
  v_re = u_re * c_m_re - -0.0 * c_m_im;
  q_im = u_re * c_m_im + -0.0 * c_m_re;
  u_re = 3.0 * (m * m);
  r_im = 3.0 * (m * 0.0 + 0.0 * m);
  w_re = u_re * R - r_im * 0.0;
  r_im = u_re * 0.0 + r_im * R;
  u_re = -2.0 * d;
  x_re = -2.0 * d;
  y_re = x_re * del;
  s_im = x_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  x_re = 2.0 * (del * del);
  t_im = 2.0 * (del * 0.0 + 0.0 * del);
  ab_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  bb_re = 4.0 * d;
  cb_re = bb_re * f_m_re - 0.0 * f_m_im;
  u_im = bb_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  bb_re = -6.0 * del;
  db_re = bb_re * f_m_re - -0.0 * f_m_im;
  v_im = bb_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  bb_re = 4.0 * (m * m);
  w_im = 4.0 * (m * 0.0 + 0.0 * m);
  eb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  fb_re = -2.0 * del;
  gb_re = fb_re * f_m_re - -0.0 * f_m_im;
  x_im = fb_re * f_m_im + -0.0 * f_m_re;
  fb_re = 3.0 * (m * m);
  y_im = 3.0 * (m * 0.0 + 0.0 * m);
  hb_re = fb_re * R - y_im * 0.0;
  y_im = fb_re * 0.0 + y_im * R;
  fb_re = -3.0 * (((((((((((((k_re * del + -2.0 * (del * del)) + (n_re * m_re -
    k_im * m_im)) + (l_re * b_m_re - l_im * b_m_im)) + o_re * R) + (q_re * R -
    n_im * 0.0)) + (s_re * R - o_im * 0.0)) + (p_re * R_re - p_im * R_im)) + d *
                       rm) + t_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (v_re *
    rm - q_im * 0.0)) + R * rm) + (w_re * rm - r_im * 0.0));
  k_im = -3.0 * ((((((((((((((k_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (n_re * m_im + k_im * m_re)) + (l_re * b_m_im + l_im * b_m_re)) +
    (o_re * 0.0 + 0.0 * R)) + (q_re * 0.0 + n_im * R)) + (s_re * 0.0 + o_im * R))
                       + (p_re * R_im + p_im * R_re)) + (d * 0.0 + 0.0 * rm)) +
                     (t_re * 0.0 + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) +
                   (v_re * 0.0 + q_im * rm)) + (R * 0.0 + 0.0 * rm)) + (w_re *
    0.0 + r_im * rm));
  k_re = ((((((((((((u_re * del + 2.0 * (del * del)) + (y_re * c_m_re - s_im *
    c_m_im)) + (x_re * e_m_re - t_im * d_m_im)) + ab_re * R) + (cb_re * R - u_im
    * 0.0)) + (db_re * R - v_im * 0.0)) + (bb_re * b_R_re - w_im * b_R_im)) + d *
              rm) + eb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (gb_re * rm -
            x_im * 0.0)) + R * rm) + (hb_re * rm - y_im * 0.0);
  l_im = (((((((((((((u_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (y_re * c_m_im + s_im * c_m_re)) + (x_re * d_m_im + t_im *
    e_m_re)) + (ab_re * 0.0 + -0.0 * R)) + (cb_re * 0.0 + u_im * R)) + (db_re *
    0.0 + v_im * R)) + (bb_re * b_R_im + w_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
             + (eb_re * 0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) +
           (gb_re * 0.0 + x_im * rm)) + (R * 0.0 + 0.0 * rm)) + (hb_re * 0.0 +
    y_im * rm);
  dc4.re = (36.0 * (cb_a_re * cb_a_re - cb_a_im * cb_a_im) + (j_re * b_d_re -
             j_im * c_d_im)) + (fb_re * k_re - k_im * l_im);
  dc4.im = (36.0 * (cb_a_re * cb_a_im + cb_a_im * cb_a_re) + (j_re * c_d_im +
             j_im * b_d_re)) + (fb_re * l_im + k_im * k_re);
  dc4 = c_power(dc4);
  dc5.re = -4.0 * dc4.re + (db_a_re * db_a_re - db_a_im * db_a_im);
  dc5.im = -4.0 * dc4.im + (db_a_re * db_a_im + db_a_im * db_a_re);
  dc4 = d_power(dc5);
  dc5.re = ((((-432.0 * dc3.re + (g_re * c_d_re - e_im * b_d_im)) + (f_re *
    u_a_re - f_im * u_a_im)) + (h_re * r_re - g_im * h_im)) + (i_re * v_a_re -
             i_im * v_a_im)) + dc4.re;
  dc5.im = ((((-432.0 * dc3.im + (g_re * b_d_im + e_im * c_d_re)) + (f_re *
    u_a_im + f_im * u_a_re)) + (h_re * h_im + g_im * r_re)) + (i_re * v_a_im +
             i_im * v_a_re)) + dc4.im;
  dc3 = f_power(dc5);
  a.re = ((((re * dc0.re - 0.0 * dc0.im) + (b_re * dc1.re - im * dc1.im)) +
           (c_re * eb_a_re - b_im * eb_a_im)) + (e_re * dc2.re - c_im * dc2.im))
    + (d_re * dc3.re - d_im * dc3.im);
  a.im = ((((re * dc0.im + 0.0 * dc0.re) + (b_re * dc1.im + im * dc1.re)) +
           (c_re * eb_a_im + b_im * eb_a_re)) + (e_re * dc2.im + c_im * dc2.re))
    + (d_re * dc3.im + d_im * dc3.re);
  if ((a.im == 0.0) && (a.re >= 0.0)) {
    a.re = rt_powd_snf(a.re, -0.5);
    a.im = 0.0;
  } else {
    b_log(&a);
    a.re *= -0.5;
    a.im *= -0.5;
    if (a.im == 0.0) {
      a.re = exp(a.re);
      a.im = 0.0;
    } else if (rtIsInf(a.im) && rtIsInf(a.re) && (a.re < 0.0)) {
      a.re = 0.0;
      a.im = 0.0;
    } else {
      r = exp(a.re / 2.0);
      a.re = r * (r * cos(a.im));
      a.im = r * (r * sin(a.im));
    }
  }

  dc0 = power(b_m);
  dc1 = power(b_R);
  re = -0.25 * dc0.re;
  im = -0.25 * dc0.im;
  b_re = re * dc1.re - im * dc1.im;
  im = re * dc1.im + im * dc1.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc0 = power(b_d);
  re = -2.0 * d;
  c_re = -2.0 * d;
  d_re = c_re * del;
  b_im = c_re * 0.0 + -0.0 * del;
  m_re = m * m;
  m_im = m * 0.0 + 0.0 * m;
  b_m_re = m * m;
  b_m_im = m * 0.0 + 0.0 * m;
  c_re = 2.0 * (del * del);
  c_im = 2.0 * (del * 0.0 + 0.0 * del);
  e_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  f_re = 4.0 * d;
  g_re = f_re * c_m_re - 0.0 * c_m_im;
  d_im = f_re * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  f_re = -6.0 * del;
  h_re = f_re * c_m_re - -0.0 * c_m_im;
  e_im = f_re * c_m_im + -0.0 * c_m_re;
  R_re = R * R;
  R_im = R * 0.0 + 0.0 * R;
  f_re = 4.0 * (m * m);
  f_im = 4.0 * (m * 0.0 + 0.0 * m);
  i_re = -2.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  b_d_re = d * c_m_re - 0.0 * c_m_im;
  b_d_im = d * c_m_im + 0.0 * c_m_re;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  j_re = -2.0 * del;
  k_re = j_re * c_m_re - -0.0 * c_m_im;
  g_im = j_re * c_m_im + -0.0 * c_m_re;
  j_re = 3.0 * (m * m);
  h_im = 3.0 * (m * 0.0 + 0.0 * m);
  l_re = j_re * R - h_im * 0.0;
  h_im = j_re * 0.0 + h_im * R;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc1 = power(b_d);
  j_re = 6.0 * ((d + -del) + R);
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  n_re = -del;
  o_re = n_re * m;
  i_im = n_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  n_re = -del;
  p_re = n_re * m;
  j_im = n_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  b_d.re = (((d_d_re * R - d_d_im * 0.0) + (p_re * R - j_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  b_d.im = (((d_d_re * 0.0 + d_d_im * R) + (p_re * 0.0 + j_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  dc2 = power(b_d);
  n_re = -2.0 * (((c_d_re * R - c_d_im * 0.0) + (o_re * R - i_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  i_im = -2.0 * (((c_d_re * 0.0 + c_d_im * R) + (o_re * 0.0 + i_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  dc3 = b_power(b_m);
  dc4 = b_power(b_R);
  o_re = 0.25 * dc3.re;
  j_im = 0.25 * dc3.im;
  p_re = o_re * dc4.re - j_im * dc4.im;
  j_im = o_re * dc4.im + j_im * dc4.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc3 = b_power(b_d);
  o_re = p_re * dc3.re - j_im * dc3.im;
  j_im = p_re * dc3.im + j_im * dc3.re;
  u_a_re = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  dc3 = power(b_m);
  dc4 = power(b_R);
  p_re = 0.41997368329829105 * dc3.re;
  k_im = 0.41997368329829105 * dc3.im;
  q_re = p_re * dc4.re - k_im * dc4.im;
  k_im = p_re * dc4.im + k_im * dc4.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc3 = power(b_d);
  p_re = q_re * dc3.re - k_im * dc3.im;
  k_im = q_re * dc3.im + k_im * dc3.re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  q_re = -del;
  r_re = q_re * m;
  l_im = q_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  q_re = -m;
  s_re = q_re * R;
  n_im = q_re * 0.0 + -0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  q_re = -del;
  t_re = q_re * m;
  o_im = q_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  q_re = 12.0 * ((((c_d_re * R - c_d_im * 0.0) + (r_re * R - l_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (s_re * rm - n_im * 0.0));
  l_im = 12.0 * ((((c_d_re * 0.0 + c_d_im * R) + (r_re * 0.0 + l_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (s_re * 0.0 + n_im * rm));
  c_d_re = (((d_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  r_re = 2.0 * d;
  s_re = 2.0 * d;
  t_re = s_re * del;
  n_im = s_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  s_re = -2.0 * (del * del);
  o_im = -2.0 * (del * 0.0 + 0.0 * del);
  u_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  v_re = -4.0 * d;
  w_re = v_re * f_m_re - -0.0 * f_m_im;
  p_im = v_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  v_re = 6.0 * del;
  x_re = v_re * f_m_re - 0.0 * f_m_im;
  q_im = v_re * f_m_im + 0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  v_re = -4.0 * (m * m);
  r_im = -4.0 * (m * 0.0 + 0.0 * m);
  y_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * f_m_re - 0.0 * f_m_im;
  d_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = -2.0 * del;
  bb_re = ab_re * f_m_re - -0.0 * f_m_im;
  s_im = ab_re * f_m_im + -0.0 * f_m_re;
  ab_re = 3.0 * (m * m);
  t_im = 3.0 * (m * 0.0 + 0.0 * m);
  cb_re = ab_re * R - t_im * 0.0;
  t_im = ab_re * 0.0 + t_im * R;
  ab_re = -2.0 * d;
  db_re = -2.0 * d;
  eb_re = db_re * del;
  u_im = db_re * 0.0 + -0.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_a_im = m * m;
  db_a_re = m * 0.0 + 0.0 * m;
  db_re = 2.0 * (del * del);
  v_im = 2.0 * (del * 0.0 + 0.0 * del);
  fb_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  gb_re = 4.0 * d;
  hb_re = gb_re * r - 0.0 * fb_a_re;
  w_im = gb_re * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  gb_re = -6.0 * del;
  db_a_im = gb_re * r - -0.0 * fb_a_re;
  x_im = gb_re * fb_a_re + -0.0 * r;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  gb_re = 4.0 * (m * m);
  y_im = 4.0 * (m * 0.0 + 0.0 * m);
  eb_a_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  e_d_re = d * r - 0.0 * fb_a_re;
  e_d_im = d * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  eb_a_im = -2.0 * del;
  ib_re = eb_a_im * r - -0.0 * fb_a_re;
  ab_im = eb_a_im * fb_a_re + -0.0 * r;
  eb_a_im = 3.0 * (m * m);
  bb_im = 3.0 * (m * 0.0 + 0.0 * m);
  jb_re = eb_a_im * R - bb_im * 0.0;
  bb_im = eb_a_im * 0.0 + bb_im * R;
  eb_a_im = -3.0 * (((((((((((((r_re * del + -2.0 * (del * del)) + (t_re *
    c_m_re - n_im * c_m_im)) + (s_re * e_m_re - o_im * d_m_im)) + u_re * R) +
    (w_re * R - p_im * 0.0)) + (x_re * R - q_im * 0.0)) + (v_re * b_R_re - r_im *
    b_R_im)) + d * rm) + y_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (bb_re *
    rm - s_im * 0.0)) + R * rm) + (cb_re * rm - t_im * 0.0));
  n_im = -3.0 * ((((((((((((((r_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (t_re * c_m_im + n_im * c_m_re)) + (s_re * d_m_im + o_im * e_m_re))
    + (u_re * 0.0 + 0.0 * R)) + (w_re * 0.0 + p_im * R)) + (x_re * 0.0 + q_im *
    R)) + (v_re * b_R_im + r_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) + (y_re * 0.0
    + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (bb_re * 0.0 + s_im * rm)) +
                  (R * 0.0 + 0.0 * rm)) + (cb_re * 0.0 + t_im * rm));
  r_re = ((((((((((((ab_re * del + 2.0 * (del * del)) + (eb_re * f_m_re - u_im *
    f_m_im)) + (db_re * cb_a_im - v_im * db_a_re)) + fb_re * R) + (hb_re * R -
    w_im * 0.0)) + (db_a_im * R - x_im * 0.0)) + (gb_re * v_a_im - y_im *
    cb_a_re)) + d * rm) + eb_a_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (ib_re
            * rm - ab_im * 0.0)) + R * rm) + (jb_re * rm - bb_im * 0.0);
  o_im = (((((((((((((ab_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (eb_re * f_m_im + u_im * f_m_re)) + (db_re * db_a_re +
    v_im * cb_a_im)) + (fb_re * 0.0 + -0.0 * R)) + (hb_re * 0.0 + w_im * R)) +
                (db_a_im * 0.0 + x_im * R)) + (gb_re * cb_a_re + y_im * v_a_im))
              + (d * 0.0 + 0.0 * rm)) + (eb_a_re * 0.0 + -0.0 * rm)) + (e_d_re *
             0.0 + e_d_im * rm)) + (ib_re * 0.0 + ab_im * rm)) + (R * 0.0 + 0.0 *
           rm)) + (jb_re * 0.0 + bb_im * rm);
  s_re = (36.0 * (b_a_re * b_a_re - b_a_im * b_a_im) + (q_re * c_d_re - l_im *
           c_d_im)) + (eb_a_im * r_re - n_im * o_im);
  l_im = (36.0 * (b_a_re * b_a_im + b_a_im * b_a_re) + (q_re * c_d_im + l_im *
           c_d_re)) + (eb_a_im * o_im + n_im * r_re);
  q_re = p_re * s_re - k_im * l_im;
  k_im = p_re * l_im + k_im * s_re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  r_re = p_re * m;
  l_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((c_d_re * R - c_d_im * 0.0) + (r_re * R - l_im * 0.0)) + (m * b_R_re
    - 0.0 * b_R_im);
  b_d.im = ((c_d_re * 0.0 + c_d_im * R) + (r_re * 0.0 + l_im * R)) + (m * b_R_im
    + 0.0 * b_R_re);
  dc3 = c_power(b_d);
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  r_re = p_re * m;
  l_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  s_re = p_re * m;
  n_im = p_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  p_re = -m;
  t_re = p_re * R;
  o_im = p_re * 0.0 + -0.0 * R;
  p_re = 432.0 * (((c_d_re * R - c_d_im * 0.0) + (r_re * R - l_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  l_im = 432.0 * (((c_d_re * 0.0 + c_d_im * R) + (r_re * 0.0 + l_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  c_d_re = (((d_d_re * R - d_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (t_re * rm - o_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (t_re * 0.0 + o_im * rm);
  r_re = p_re * c_d_re - l_im * c_d_im;
  l_im = p_re * c_d_im + l_im * c_d_re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  s_re = p_re * m;
  n_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  d_d_re = (((c_d_re * R - c_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0);
  c_d_im = (((c_d_re * 0.0 + c_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  s_re = p_re * m;
  n_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  a_re = c_a_re * c_a_re - c_a_im * c_a_im;
  b_a_im = c_a_re * c_a_im + c_a_im * c_a_re;
  p_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0));
  n_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  s_re = -del;
  t_re = s_re * m;
  o_im = s_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  s_re = 2.0 * d;
  u_re = 2.0 * d;
  v_re = u_re * del;
  p_im = u_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  u_re = -2.0 * (del * del);
  q_im = -2.0 * (del * 0.0 + 0.0 * del);
  w_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = -4.0 * d;
  y_re = x_re * f_m_re - -0.0 * f_m_im;
  r_im = x_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = 6.0 * del;
  ab_re = x_re * f_m_re - 0.0 * f_m_im;
  s_im = x_re * f_m_im + 0.0 * f_m_re;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  x_re = -4.0 * (m * m);
  t_im = -4.0 * (m * 0.0 + 0.0 * m);
  bb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_re = -2.0 * del;
  db_re = cb_re * f_m_re - -0.0 * f_m_im;
  u_im = cb_re * f_m_im + -0.0 * f_m_re;
  cb_re = 3.0 * (m * m);
  v_im = 3.0 * (m * 0.0 + 0.0 * m);
  eb_re = cb_re * R - v_im * 0.0;
  v_im = cb_re * 0.0 + v_im * R;
  cb_re = 54.0 * (((c_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  o_im = 54.0 * (((c_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  t_re = ((((((((((((s_re * del + -2.0 * (del * del)) + (v_re * c_m_re - p_im *
    c_m_im)) + (u_re * e_m_re - q_im * d_m_im)) + w_re * R) + (y_re * R - r_im *
    0.0)) + (ab_re * R - s_im * 0.0)) + (x_re * v_a_im - t_im * cb_a_re)) + d *
              rm) + bb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (db_re * rm -
            u_im * 0.0)) + R * rm) + (eb_re * rm - v_im * 0.0);
  p_im = (((((((((((((s_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (v_re * c_m_im + p_im * c_m_re)) + (u_re * d_m_im + q_im *
    e_m_re)) + (w_re * 0.0 + 0.0 * R)) + (y_re * 0.0 + r_im * R)) + (ab_re * 0.0
    + s_im * R)) + (x_re * cb_a_re + t_im * v_a_im)) + (d * 0.0 + 0.0 * rm)) +
             (bb_re * 0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) + (db_re
            * 0.0 + u_im * rm)) + (R * 0.0 + 0.0 * rm)) + (eb_re * 0.0 + v_im *
    rm);
  s_re = cb_re * t_re - o_im * p_im;
  o_im = cb_re * p_im + o_im * t_re;
  t_re = -2.0 * d;
  u_re = -2.0 * d;
  v_re = u_re * del;
  p_im = u_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  u_re = 2.0 * (del * del);
  q_im = 2.0 * (del * 0.0 + 0.0 * del);
  w_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = 4.0 * d;
  y_re = x_re * f_m_re - 0.0 * f_m_im;
  r_im = x_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = -6.0 * del;
  ab_re = x_re * f_m_re - -0.0 * f_m_im;
  s_im = x_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  x_re = 4.0 * (m * m);
  t_im = 4.0 * (m * 0.0 + 0.0 * m);
  bb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * f_m_re - 0.0 * f_m_im;
  d_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_re = -2.0 * del;
  db_re = cb_re * f_m_re - -0.0 * f_m_im;
  u_im = cb_re * f_m_im + -0.0 * f_m_re;
  cb_re = 3.0 * (m * m);
  v_im = 3.0 * (m * 0.0 + 0.0 * m);
  eb_re = cb_re * R - v_im * 0.0;
  v_im = cb_re * 0.0 + v_im * R;
  cb_re = ((((((((((((t_re * del + 2.0 * (del * del)) + (v_re * c_m_re - p_im *
    c_m_im)) + (u_re * e_m_re - q_im * d_m_im)) + w_re * R) + (y_re * R - r_im *
    0.0)) + (ab_re * R - s_im * 0.0)) + (x_re * b_R_re - t_im * b_R_im)) + d *
               rm) + bb_re * rm) + (c_d_re * rm - d_d_im * 0.0)) + (db_re * rm -
             u_im * 0.0)) + R * rm) + (eb_re * rm - v_im * 0.0);
  p_im = (((((((((((((t_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (v_re * c_m_im + p_im * c_m_re)) + (u_re * d_m_im + q_im *
    e_m_re)) + (w_re * 0.0 + -0.0 * R)) + (y_re * 0.0 + r_im * R)) + (ab_re *
    0.0 + s_im * R)) + (x_re * b_R_im + t_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
             + (bb_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + d_d_im * rm)) +
           (db_re * 0.0 + u_im * rm)) + (R * 0.0 + 0.0 * rm)) + (eb_re * 0.0 +
    v_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  t_re = -del;
  u_re = t_re * m;
  q_im = t_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  t_re = -m;
  v_re = t_re * R;
  r_im = t_re * 0.0 + -0.0 * R;
  b_a_re = d_a_re * d_a_re - d_a_im * d_a_im;
  c_a_im = d_a_re * d_a_im + d_a_im * d_a_re;
  t_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (u_re * R - q_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (v_re * rm - r_im * 0.0));
  q_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (u_re * 0.0 + q_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (v_re * 0.0 + r_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  u_re = -del;
  v_re = u_re * m;
  r_im = u_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  u_re = -m;
  w_re = u_re * R;
  s_im = u_re * 0.0 + -0.0 * R;
  e_d_re = d * m;
  e_d_im = d * 0.0 + 0.0 * m;
  u_re = -del;
  x_re = u_re * m;
  t_im = u_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  u_re = 12.0 * ((((c_d_re * R - d_d_im * 0.0) + (v_re * R - r_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (w_re * rm - s_im * 0.0));
  r_im = 12.0 * ((((c_d_re * 0.0 + d_d_im * R) + (v_re * 0.0 + r_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (w_re * 0.0 + s_im * rm));
  c_d_re = (((e_d_re * R - e_d_im * 0.0) + (x_re * R - t_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  d_d_im = (((e_d_re * 0.0 + e_d_im * R) + (x_re * 0.0 + t_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  v_re = 2.0 * d;
  w_re = 2.0 * d;
  x_re = w_re * del;
  s_im = w_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  w_re = -2.0 * (del * del);
  t_im = -2.0 * (del * 0.0 + 0.0 * del);
  y_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = -4.0 * d;
  bb_re = ab_re * f_m_re - -0.0 * f_m_im;
  u_im = ab_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = 6.0 * del;
  db_re = ab_re * f_m_re - 0.0 * f_m_im;
  v_im = ab_re * f_m_im + 0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  ab_re = -4.0 * (m * m);
  w_im = -4.0 * (m * 0.0 + 0.0 * m);
  eb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  fb_re = -2.0 * del;
  gb_re = fb_re * f_m_re - -0.0 * f_m_im;
  x_im = fb_re * f_m_im + -0.0 * f_m_re;
  fb_re = 3.0 * (m * m);
  y_im = 3.0 * (m * 0.0 + 0.0 * m);
  hb_re = fb_re * R - y_im * 0.0;
  y_im = fb_re * 0.0 + y_im * R;
  fb_re = -2.0 * d;
  db_a_im = -2.0 * d;
  eb_a_re = db_a_im * del;
  ab_im = db_a_im * 0.0 + -0.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_a_im = m * m;
  db_a_re = m * 0.0 + 0.0 * m;
  db_a_im = 2.0 * (del * del);
  bb_im = 2.0 * (del * 0.0 + 0.0 * del);
  eb_a_im = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ib_re = 4.0 * d;
  jb_re = ib_re * r - 0.0 * fb_a_re;
  w_a_im = ib_re * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ib_re = -6.0 * del;
  x_a_re = ib_re * r - -0.0 * fb_a_re;
  x_a_im = ib_re * fb_a_re + -0.0 * r;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  ib_re = 4.0 * (m * m);
  y_a_re = 4.0 * (m * 0.0 + 0.0 * m);
  y_a_im = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ab_a_re = d * r - 0.0 * fb_a_re;
  ab_a_im = d * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  bb_a_re = -2.0 * del;
  bb_a_im = bb_a_re * r - -0.0 * fb_a_re;
  v_a_re = bb_a_re * fb_a_re + -0.0 * r;
  bb_a_re = 3.0 * (m * m);
  w_a_re = 3.0 * (m * 0.0 + 0.0 * m);
  u_a_im = bb_a_re * R - w_a_re * 0.0;
  w_a_re = bb_a_re * 0.0 + w_a_re * R;
  bb_a_re = -3.0 * (((((((((((((v_re * del + -2.0 * (del * del)) + (x_re *
    c_m_re - s_im * c_m_im)) + (w_re * e_m_re - t_im * d_m_im)) + y_re * R) +
    (bb_re * R - u_im * 0.0)) + (db_re * R - v_im * 0.0)) + (ab_re * b_R_re -
    w_im * b_R_im)) + d * rm) + eb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) +
                      (gb_re * rm - x_im * 0.0)) + R * rm) + (hb_re * rm - y_im *
    0.0));
  s_im = -3.0 * ((((((((((((((v_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (x_re * c_m_im + s_im * c_m_re)) + (w_re * d_m_im + t_im * e_m_re))
    + (y_re * 0.0 + 0.0 * R)) + (bb_re * 0.0 + u_im * R)) + (db_re * 0.0 + v_im *
    R)) + (ab_re * b_R_im + w_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) + (eb_re *
    0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) + (gb_re * 0.0 + x_im * rm))
                  + (R * 0.0 + 0.0 * rm)) + (hb_re * 0.0 + y_im * rm));
  v_re = ((((((((((((fb_re * del + 2.0 * (del * del)) + (eb_a_re * f_m_re -
    ab_im * f_m_im)) + (db_a_im * cb_a_im - bb_im * db_a_re)) + eb_a_im * R) +
                 (jb_re * R - w_a_im * 0.0)) + (x_a_re * R - x_a_im * 0.0)) +
               (ib_re * v_a_im - y_a_re * cb_a_re)) + d * rm) + y_a_im * rm) +
            (ab_a_re * rm - ab_a_im * 0.0)) + (bb_a_im * rm - v_a_re * 0.0)) + R
          * rm) + (u_a_im * rm - w_a_re * 0.0);
  t_im = (((((((((((((fb_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (eb_a_re * f_m_im + ab_im * f_m_re)) + (db_a_im * db_a_re
    + bb_im * cb_a_im)) + (eb_a_im * 0.0 + -0.0 * R)) + (jb_re * 0.0 + w_a_im *
    R)) + (x_a_re * 0.0 + x_a_im * R)) + (ib_re * cb_a_re + y_a_re * v_a_im)) +
              (d * 0.0 + 0.0 * rm)) + (y_a_im * 0.0 + -0.0 * rm)) + (ab_a_re *
             0.0 + ab_a_im * rm)) + (bb_a_im * 0.0 + v_a_re * rm)) + (R * 0.0 +
           0.0 * rm)) + (u_a_im * 0.0 + w_a_re * rm);
  dc4.re = (36.0 * (e_a_re * e_a_re - e_a_im * e_a_im) + (u_re * c_d_re - r_im *
             d_d_im)) + (bb_a_re * v_re - s_im * t_im);
  dc4.im = (36.0 * (e_a_re * e_a_im + e_a_im * e_a_re) + (u_re * d_d_im + r_im *
             c_d_re)) + (bb_a_re * t_im + s_im * v_re);
  dc4 = c_power(dc4);
  dc5.re = -4.0 * dc4.re + (f_a_re * f_a_re - f_a_im * f_a_im);
  dc5.im = -4.0 * dc4.im + (f_a_re * f_a_im + f_a_im * f_a_re);
  dc4 = d_power(dc5);
  dc5.re = ((((-432.0 * dc3.re + (r_re * d_d_re - l_im * c_d_im)) + (p_re * a_re
    - n_im * b_a_im)) + (s_re * cb_re - o_im * p_im)) + (t_re * b_a_re - q_im *
             c_a_im)) + dc4.re;
  dc5.im = ((((-432.0 * dc3.im + (r_re * c_d_im + l_im * d_d_re)) + (p_re *
    b_a_im + n_im * a_re)) + (s_re * p_im + o_im * cb_re)) + (t_re * c_a_im +
             q_im * b_a_re)) + dc4.im;
  dc3 = e_power(dc5);
  dc4 = power(b_m);
  dc5 = power(b_R);
  p_re = 0.26456684199469993 * dc4.re;
  l_im = 0.26456684199469993 * dc4.im;
  r_re = p_re * dc5.re - l_im * dc5.im;
  l_im = p_re * dc5.im + l_im * dc5.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc4 = power(b_d);
  p_re = r_re * dc4.re - l_im * dc4.im;
  l_im = r_re * dc4.im + l_im * dc4.re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  s_re = r_re * m;
  n_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((c_d_re * R - c_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m * b_R_re
    - 0.0 * b_R_im);
  b_d.im = ((c_d_re * 0.0 + c_d_im * R) + (s_re * 0.0 + n_im * R)) + (m * b_R_im
    + 0.0 * b_R_re);
  dc4 = c_power(b_d);
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  s_re = r_re * m;
  n_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  t_re = r_re * m;
  o_im = r_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  r_re = -m;
  u_re = r_re * R;
  p_im = r_re * 0.0 + -0.0 * R;
  r_re = 432.0 * (((c_d_re * R - c_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  n_im = 432.0 * (((c_d_re * 0.0 + c_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  c_d_re = (((d_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (u_re * rm - p_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (u_re * 0.0 + p_im * rm);
  s_re = r_re * c_d_re - n_im * c_d_im;
  n_im = r_re * c_d_im + n_im * c_d_re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  t_re = r_re * m;
  o_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  d_d_re = (((c_d_re * R - c_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0);
  c_d_im = (((c_d_re * 0.0 + c_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  t_re = r_re * m;
  o_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  a_re = g_a_re * g_a_re - g_a_im * g_a_im;
  b_a_im = g_a_re * g_a_im + g_a_im * g_a_re;
  r_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0));
  o_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  t_re = -del;
  u_re = t_re * m;
  p_im = t_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  t_re = 2.0 * d;
  v_re = 2.0 * d;
  w_re = v_re * del;
  q_im = v_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  v_re = -2.0 * (del * del);
  r_im = -2.0 * (del * 0.0 + 0.0 * del);
  x_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = -4.0 * d;
  ab_re = y_re * f_m_re - -0.0 * f_m_im;
  s_im = y_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = 6.0 * del;
  bb_re = y_re * f_m_re - 0.0 * f_m_im;
  t_im = y_re * f_m_im + 0.0 * f_m_re;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  y_re = -4.0 * (m * m);
  u_im = -4.0 * (m * 0.0 + 0.0 * m);
  cb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  db_re = -2.0 * del;
  eb_re = db_re * f_m_re - -0.0 * f_m_im;
  v_im = db_re * f_m_im + -0.0 * f_m_re;
  db_re = 3.0 * (m * m);
  w_im = 3.0 * (m * 0.0 + 0.0 * m);
  fb_re = db_re * R - w_im * 0.0;
  w_im = db_re * 0.0 + w_im * R;
  db_re = 54.0 * (((c_d_re * R - d_d_im * 0.0) + (u_re * R - p_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  p_im = 54.0 * (((c_d_re * 0.0 + d_d_im * R) + (u_re * 0.0 + p_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  u_re = ((((((((((((t_re * del + -2.0 * (del * del)) + (w_re * c_m_re - q_im *
    c_m_im)) + (v_re * e_m_re - r_im * d_m_im)) + x_re * R) + (ab_re * R - s_im *
    0.0)) + (bb_re * R - t_im * 0.0)) + (y_re * v_a_im - u_im * cb_a_re)) + d *
              rm) + cb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (eb_re * rm -
            v_im * 0.0)) + R * rm) + (fb_re * rm - w_im * 0.0);
  q_im = (((((((((((((t_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (w_re * c_m_im + q_im * c_m_re)) + (v_re * d_m_im + r_im *
    e_m_re)) + (x_re * 0.0 + 0.0 * R)) + (ab_re * 0.0 + s_im * R)) + (bb_re *
    0.0 + t_im * R)) + (y_re * cb_a_re + u_im * v_a_im)) + (d * 0.0 + 0.0 * rm))
             + (cb_re * 0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) +
           (eb_re * 0.0 + v_im * rm)) + (R * 0.0 + 0.0 * rm)) + (fb_re * 0.0 +
    w_im * rm);
  t_re = db_re * u_re - p_im * q_im;
  p_im = db_re * q_im + p_im * u_re;
  u_re = -2.0 * d;
  v_re = -2.0 * d;
  w_re = v_re * del;
  q_im = v_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  v_re = 2.0 * (del * del);
  r_im = 2.0 * (del * 0.0 + 0.0 * del);
  x_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = 4.0 * d;
  ab_re = y_re * f_m_re - 0.0 * f_m_im;
  s_im = y_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = -6.0 * del;
  bb_re = y_re * f_m_re - -0.0 * f_m_im;
  t_im = y_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  y_re = 4.0 * (m * m);
  u_im = 4.0 * (m * 0.0 + 0.0 * m);
  cb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * f_m_re - 0.0 * f_m_im;
  d_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  db_re = -2.0 * del;
  eb_re = db_re * f_m_re - -0.0 * f_m_im;
  v_im = db_re * f_m_im + -0.0 * f_m_re;
  db_re = 3.0 * (m * m);
  w_im = 3.0 * (m * 0.0 + 0.0 * m);
  fb_re = db_re * R - w_im * 0.0;
  w_im = db_re * 0.0 + w_im * R;
  db_re = ((((((((((((u_re * del + 2.0 * (del * del)) + (w_re * c_m_re - q_im *
    c_m_im)) + (v_re * e_m_re - r_im * d_m_im)) + x_re * R) + (ab_re * R - s_im *
    0.0)) + (bb_re * R - t_im * 0.0)) + (y_re * b_R_re - u_im * b_R_im)) + d *
               rm) + cb_re * rm) + (c_d_re * rm - d_d_im * 0.0)) + (eb_re * rm -
             v_im * 0.0)) + R * rm) + (fb_re * rm - w_im * 0.0);
  q_im = (((((((((((((u_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (w_re * c_m_im + q_im * c_m_re)) + (v_re * d_m_im + r_im *
    e_m_re)) + (x_re * 0.0 + -0.0 * R)) + (ab_re * 0.0 + s_im * R)) + (bb_re *
    0.0 + t_im * R)) + (y_re * b_R_im + u_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
             + (cb_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + d_d_im * rm)) +
           (eb_re * 0.0 + v_im * rm)) + (R * 0.0 + 0.0 * rm)) + (fb_re * 0.0 +
    w_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  u_re = -del;
  v_re = u_re * m;
  r_im = u_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  u_re = -m;
  w_re = u_re * R;
  s_im = u_re * 0.0 + -0.0 * R;
  b_a_re = h_a_re * h_a_re - h_a_im * h_a_im;
  c_a_im = h_a_re * h_a_im + h_a_im * h_a_re;
  u_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (v_re * R - r_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (w_re * rm - s_im * 0.0));
  r_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (v_re * 0.0 + r_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (w_re * 0.0 + s_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  v_re = -del;
  w_re = v_re * m;
  s_im = v_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  v_re = -m;
  x_re = v_re * R;
  t_im = v_re * 0.0 + -0.0 * R;
  e_d_re = d * m;
  e_d_im = d * 0.0 + 0.0 * m;
  v_re = -del;
  y_re = v_re * m;
  u_im = v_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  v_re = 12.0 * ((((c_d_re * R - d_d_im * 0.0) + (w_re * R - s_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (x_re * rm - t_im * 0.0));
  s_im = 12.0 * ((((c_d_re * 0.0 + d_d_im * R) + (w_re * 0.0 + s_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (x_re * 0.0 + t_im * rm));
  c_d_re = (((e_d_re * R - e_d_im * 0.0) + (y_re * R - u_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  d_d_im = (((e_d_re * 0.0 + e_d_im * R) + (y_re * 0.0 + u_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  w_re = 2.0 * d;
  x_re = 2.0 * d;
  y_re = x_re * del;
  t_im = x_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  x_re = -2.0 * (del * del);
  u_im = -2.0 * (del * 0.0 + 0.0 * del);
  ab_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  bb_re = -4.0 * d;
  cb_re = bb_re * f_m_re - -0.0 * f_m_im;
  v_im = bb_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  bb_re = 6.0 * del;
  eb_re = bb_re * f_m_re - 0.0 * f_m_im;
  w_im = bb_re * f_m_im + 0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  bb_re = -4.0 * (m * m);
  x_im = -4.0 * (m * 0.0 + 0.0 * m);
  fb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  gb_re = -2.0 * del;
  hb_re = gb_re * f_m_re - -0.0 * f_m_im;
  y_im = gb_re * f_m_im + -0.0 * f_m_re;
  gb_re = 3.0 * (m * m);
  ab_im = 3.0 * (m * 0.0 + 0.0 * m);
  db_a_im = gb_re * R - ab_im * 0.0;
  ab_im = gb_re * 0.0 + ab_im * R;
  gb_re = -2.0 * d;
  eb_a_re = -2.0 * d;
  eb_a_im = eb_a_re * del;
  bb_im = eb_a_re * 0.0 + -0.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_a_im = m * m;
  db_a_re = m * 0.0 + 0.0 * m;
  eb_a_re = 2.0 * (del * del);
  w_a_im = 2.0 * (del * 0.0 + 0.0 * del);
  ib_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  jb_re = 4.0 * d;
  x_a_re = jb_re * r - 0.0 * fb_a_re;
  x_a_im = jb_re * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  jb_re = -6.0 * del;
  y_a_im = jb_re * r - -0.0 * fb_a_re;
  y_a_re = jb_re * fb_a_re + -0.0 * r;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  jb_re = 4.0 * (m * m);
  v_a_re = 4.0 * (m * 0.0 + 0.0 * m);
  bb_a_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ab_a_re = d * r - 0.0 * fb_a_re;
  ab_a_im = d * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  bb_a_im = -2.0 * del;
  u_a_im = bb_a_im * r - -0.0 * fb_a_re;
  w_a_re = bb_a_im * fb_a_re + -0.0 * r;
  bb_a_im = 3.0 * (m * m);
  r = 3.0 * (m * 0.0 + 0.0 * m);
  fb_a_re = bb_a_im * R - r * 0.0;
  r = bb_a_im * 0.0 + r * R;
  bb_a_im = -3.0 * (((((((((((((w_re * del + -2.0 * (del * del)) + (y_re *
    c_m_re - t_im * c_m_im)) + (x_re * e_m_re - u_im * d_m_im)) + ab_re * R) +
    (cb_re * R - v_im * 0.0)) + (eb_re * R - w_im * 0.0)) + (bb_re * b_R_re -
    x_im * b_R_im)) + d * rm) + fb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) +
                      (hb_re * rm - y_im * 0.0)) + R * rm) + (db_a_im * rm -
    ab_im * 0.0));
  t_im = -3.0 * ((((((((((((((w_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (y_re * c_m_im + t_im * c_m_re)) + (x_re * d_m_im + u_im * e_m_re))
    + (ab_re * 0.0 + 0.0 * R)) + (cb_re * 0.0 + v_im * R)) + (eb_re * 0.0 + w_im
    * R)) + (bb_re * b_R_im + x_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) + (fb_re *
    0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) + (hb_re * 0.0 + y_im * rm))
                  + (R * 0.0 + 0.0 * rm)) + (db_a_im * 0.0 + ab_im * rm));
  w_re = ((((((((((((gb_re * del + 2.0 * (del * del)) + (eb_a_im * f_m_re -
    bb_im * f_m_im)) + (eb_a_re * cb_a_im - w_a_im * db_a_re)) + ib_re * R) +
                 (x_a_re * R - x_a_im * 0.0)) + (y_a_im * R - y_a_re * 0.0)) +
               (jb_re * v_a_im - v_a_re * cb_a_re)) + d * rm) + bb_a_re * rm) +
            (ab_a_re * rm - ab_a_im * 0.0)) + (u_a_im * rm - w_a_re * 0.0)) + R *
          rm) + (fb_a_re * rm - r * 0.0);
  u_im = (((((((((((((gb_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (eb_a_im * f_m_im + bb_im * f_m_re)) + (eb_a_re * db_a_re
    + w_a_im * cb_a_im)) + (ib_re * 0.0 + -0.0 * R)) + (x_a_re * 0.0 + x_a_im *
    R)) + (y_a_im * 0.0 + y_a_re * R)) + (jb_re * cb_a_re + v_a_re * v_a_im)) +
              (d * 0.0 + 0.0 * rm)) + (bb_a_re * 0.0 + -0.0 * rm)) + (ab_a_re *
             0.0 + ab_a_im * rm)) + (u_a_im * 0.0 + w_a_re * rm)) + (R * 0.0 +
           0.0 * rm)) + (fb_a_re * 0.0 + r * rm);
  dc5.re = (36.0 * (i_a_re * i_a_re - i_a_im * i_a_im) + (v_re * c_d_re - s_im *
             d_d_im)) + (bb_a_im * w_re - t_im * u_im);
  dc5.im = (36.0 * (i_a_re * i_a_im + i_a_im * i_a_re) + (v_re * d_d_im + s_im *
             c_d_re)) + (bb_a_im * u_im + t_im * w_re);
  dc5 = c_power(dc5);
  b_d.re = -4.0 * dc5.re + (j_a_re * j_a_re - j_a_im * j_a_im);
  b_d.im = -4.0 * dc5.im + (j_a_re * j_a_im + j_a_im * j_a_re);
  dc5 = d_power(b_d);
  b_d.re = ((((-432.0 * dc4.re + (s_re * d_d_re - n_im * c_d_im)) + (r_re * a_re
    - o_im * b_a_im)) + (t_re * db_re - p_im * q_im)) + (u_re * b_a_re - r_im *
             c_a_im)) + dc5.re;
  b_d.im = ((((-432.0 * dc4.im + (s_re * c_d_im + n_im * d_d_re)) + (r_re *
    b_a_im + o_im * a_re)) + (t_re * q_im + p_im * db_re)) + (u_re * c_a_im +
             r_im * b_a_re)) + dc5.im;
  dc4 = f_power(b_d);
  dc5.re = ((((j_re * dc1.re - 0.0 * dc1.im) + (n_re * dc2.re - i_im * dc2.im))
             + (o_re * u_a_re - j_im * a_im)) + (q_re * dc3.re - k_im * dc3.im))
    + (p_re * dc4.re - l_im * dc4.im);
  dc5.im = ((((j_re * dc1.im + 0.0 * dc1.re) + (n_re * dc2.im + i_im * dc2.re))
             + (o_re * a_im + j_im * u_a_re)) + (q_re * dc3.im + k_im * dc3.re))
    + (p_re * dc4.im + l_im * dc4.re);
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc1 = power(b_d);
  j_re = 6.0 * ((d + -del) + R);
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  n_re = -del;
  o_re = n_re * m;
  i_im = n_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  n_re = -del;
  p_re = n_re * m;
  j_im = n_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  b_d.re = (((d_d_re * R - d_d_im * 0.0) + (p_re * R - j_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  b_d.im = (((d_d_re * 0.0 + d_d_im * R) + (p_re * 0.0 + j_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  dc2 = power(b_d);
  n_re = 2.0 * (((c_d_re * R - c_d_im * 0.0) + (o_re * R - i_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  i_im = 2.0 * (((c_d_re * 0.0 + c_d_im * R) + (o_re * 0.0 + i_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  dc3 = b_power(b_m);
  dc4 = b_power(b_R);
  o_re = 0.5 * dc3.re;
  j_im = 0.5 * dc3.im;
  p_re = o_re * dc4.re - j_im * dc4.im;
  j_im = o_re * dc4.im + j_im * dc4.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc3 = b_power(b_d);
  o_re = p_re * dc3.re - j_im * dc3.im;
  j_im = p_re * dc3.im + j_im * dc3.re;
  a_re = k_a_re * k_a_re - k_a_im * k_a_im;
  a_im = k_a_re * k_a_im + k_a_im * k_a_re;
  dc3 = power(b_m);
  dc4 = power(b_R);
  p_re = -0.41997368329829105 * dc3.re;
  k_im = -0.41997368329829105 * dc3.im;
  q_re = p_re * dc4.re - k_im * dc4.im;
  k_im = p_re * dc4.im + k_im * dc4.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc3 = power(b_d);
  p_re = q_re * dc3.re - k_im * dc3.im;
  k_im = q_re * dc3.im + k_im * dc3.re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  q_re = -del;
  r_re = q_re * m;
  l_im = q_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  q_re = -m;
  s_re = q_re * R;
  n_im = q_re * 0.0 + -0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  q_re = -del;
  t_re = q_re * m;
  o_im = q_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  q_re = 12.0 * ((((c_d_re * R - c_d_im * 0.0) + (r_re * R - l_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (s_re * rm - n_im * 0.0));
  l_im = 12.0 * ((((c_d_re * 0.0 + c_d_im * R) + (r_re * 0.0 + l_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (s_re * 0.0 + n_im * rm));
  c_d_re = (((d_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  r_re = 2.0 * d;
  s_re = 2.0 * d;
  t_re = s_re * del;
  n_im = s_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  s_re = -2.0 * (del * del);
  o_im = -2.0 * (del * 0.0 + 0.0 * del);
  u_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  v_re = -4.0 * d;
  w_re = v_re * f_m_re - -0.0 * f_m_im;
  p_im = v_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  v_re = 6.0 * del;
  x_re = v_re * f_m_re - 0.0 * f_m_im;
  q_im = v_re * f_m_im + 0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  v_re = -4.0 * (m * m);
  r_im = -4.0 * (m * 0.0 + 0.0 * m);
  y_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  d_d_re = d * f_m_re - 0.0 * f_m_im;
  d_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = -2.0 * del;
  bb_re = ab_re * f_m_re - -0.0 * f_m_im;
  s_im = ab_re * f_m_im + -0.0 * f_m_re;
  ab_re = 3.0 * (m * m);
  t_im = 3.0 * (m * 0.0 + 0.0 * m);
  cb_re = ab_re * R - t_im * 0.0;
  t_im = ab_re * 0.0 + t_im * R;
  ab_re = -2.0 * d;
  db_re = -2.0 * d;
  eb_re = db_re * del;
  u_im = db_re * 0.0 + -0.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_a_im = m * m;
  db_a_re = m * 0.0 + 0.0 * m;
  db_re = 2.0 * (del * del);
  v_im = 2.0 * (del * 0.0 + 0.0 * del);
  fb_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  gb_re = 4.0 * d;
  hb_re = gb_re * r - 0.0 * fb_a_re;
  w_im = gb_re * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  gb_re = -6.0 * del;
  db_a_im = gb_re * r - -0.0 * fb_a_re;
  x_im = gb_re * fb_a_re + -0.0 * r;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  gb_re = 4.0 * (m * m);
  y_im = 4.0 * (m * 0.0 + 0.0 * m);
  eb_a_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  e_d_re = d * r - 0.0 * fb_a_re;
  e_d_im = d * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  eb_a_im = -2.0 * del;
  ib_re = eb_a_im * r - -0.0 * fb_a_re;
  ab_im = eb_a_im * fb_a_re + -0.0 * r;
  eb_a_im = 3.0 * (m * m);
  bb_im = 3.0 * (m * 0.0 + 0.0 * m);
  jb_re = eb_a_im * R - bb_im * 0.0;
  bb_im = eb_a_im * 0.0 + bb_im * R;
  eb_a_im = -3.0 * (((((((((((((r_re * del + -2.0 * (del * del)) + (t_re *
    c_m_re - n_im * c_m_im)) + (s_re * e_m_re - o_im * d_m_im)) + u_re * R) +
    (w_re * R - p_im * 0.0)) + (x_re * R - q_im * 0.0)) + (v_re * b_R_re - r_im *
    b_R_im)) + d * rm) + y_re * rm) + (d_d_re * rm - d_d_im * 0.0)) + (bb_re *
    rm - s_im * 0.0)) + R * rm) + (cb_re * rm - t_im * 0.0));
  n_im = -3.0 * ((((((((((((((r_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (t_re * c_m_im + n_im * c_m_re)) + (s_re * d_m_im + o_im * e_m_re))
    + (u_re * 0.0 + 0.0 * R)) + (w_re * 0.0 + p_im * R)) + (x_re * 0.0 + q_im *
    R)) + (v_re * b_R_im + r_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) + (y_re * 0.0
    + -0.0 * rm)) + (d_d_re * 0.0 + d_d_im * rm)) + (bb_re * 0.0 + s_im * rm)) +
                  (R * 0.0 + 0.0 * rm)) + (cb_re * 0.0 + t_im * rm));
  r_re = ((((((((((((ab_re * del + 2.0 * (del * del)) + (eb_re * f_m_re - u_im *
    f_m_im)) + (db_re * cb_a_im - v_im * db_a_re)) + fb_re * R) + (hb_re * R -
    w_im * 0.0)) + (db_a_im * R - x_im * 0.0)) + (gb_re * v_a_im - y_im *
    cb_a_re)) + d * rm) + eb_a_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (ib_re
            * rm - ab_im * 0.0)) + R * rm) + (jb_re * rm - bb_im * 0.0);
  o_im = (((((((((((((ab_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (eb_re * f_m_im + u_im * f_m_re)) + (db_re * db_a_re +
    v_im * cb_a_im)) + (fb_re * 0.0 + -0.0 * R)) + (hb_re * 0.0 + w_im * R)) +
                (db_a_im * 0.0 + x_im * R)) + (gb_re * cb_a_re + y_im * v_a_im))
              + (d * 0.0 + 0.0 * rm)) + (eb_a_re * 0.0 + -0.0 * rm)) + (e_d_re *
             0.0 + e_d_im * rm)) + (ib_re * 0.0 + ab_im * rm)) + (R * 0.0 + 0.0 *
           rm)) + (jb_re * 0.0 + bb_im * rm);
  s_re = (36.0 * (l_a_re * l_a_re - l_a_im * l_a_im) + (q_re * c_d_re - l_im *
           c_d_im)) + (eb_a_im * r_re - n_im * o_im);
  l_im = (36.0 * (l_a_re * l_a_im + l_a_im * l_a_re) + (q_re * c_d_im + l_im *
           c_d_re)) + (eb_a_im * o_im + n_im * r_re);
  q_re = p_re * s_re - k_im * l_im;
  k_im = p_re * l_im + k_im * s_re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  r_re = p_re * m;
  l_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((c_d_re * R - c_d_im * 0.0) + (r_re * R - l_im * 0.0)) + (m * b_R_re
    - 0.0 * b_R_im);
  b_d.im = ((c_d_re * 0.0 + c_d_im * R) + (r_re * 0.0 + l_im * R)) + (m * b_R_im
    + 0.0 * b_R_re);
  dc3 = c_power(b_d);
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  r_re = p_re * m;
  l_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  s_re = p_re * m;
  n_im = p_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  p_re = -m;
  t_re = p_re * R;
  o_im = p_re * 0.0 + -0.0 * R;
  p_re = 432.0 * (((c_d_re * R - c_d_im * 0.0) + (r_re * R - l_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  l_im = 432.0 * (((c_d_re * 0.0 + c_d_im * R) + (r_re * 0.0 + l_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  c_d_re = (((d_d_re * R - d_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (t_re * rm - o_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (t_re * 0.0 + o_im * rm);
  r_re = p_re * c_d_re - l_im * c_d_im;
  l_im = p_re * c_d_im + l_im * c_d_re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  s_re = p_re * m;
  n_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  d_d_re = (((c_d_re * R - c_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0);
  c_d_im = (((c_d_re * 0.0 + c_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  p_re = -del;
  s_re = p_re * m;
  n_im = p_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  b_a_re = m_a_re * m_a_re - m_a_im * m_a_im;
  b_a_im = m_a_re * m_a_im + m_a_im * m_a_re;
  p_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0));
  n_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  s_re = -del;
  t_re = s_re * m;
  o_im = s_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  s_re = 2.0 * d;
  u_re = 2.0 * d;
  v_re = u_re * del;
  p_im = u_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  u_re = -2.0 * (del * del);
  q_im = -2.0 * (del * 0.0 + 0.0 * del);
  w_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = -4.0 * d;
  y_re = x_re * f_m_re - -0.0 * f_m_im;
  r_im = x_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = 6.0 * del;
  ab_re = x_re * f_m_re - 0.0 * f_m_im;
  s_im = x_re * f_m_im + 0.0 * f_m_re;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  x_re = -4.0 * (m * m);
  t_im = -4.0 * (m * 0.0 + 0.0 * m);
  bb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_re = -2.0 * del;
  db_re = cb_re * f_m_re - -0.0 * f_m_im;
  u_im = cb_re * f_m_im + -0.0 * f_m_re;
  cb_re = 3.0 * (m * m);
  v_im = 3.0 * (m * 0.0 + 0.0 * m);
  eb_re = cb_re * R - v_im * 0.0;
  v_im = cb_re * 0.0 + v_im * R;
  cb_re = 54.0 * (((c_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  o_im = 54.0 * (((c_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  t_re = ((((((((((((s_re * del + -2.0 * (del * del)) + (v_re * c_m_re - p_im *
    c_m_im)) + (u_re * e_m_re - q_im * d_m_im)) + w_re * R) + (y_re * R - r_im *
    0.0)) + (ab_re * R - s_im * 0.0)) + (x_re * v_a_im - t_im * cb_a_re)) + d *
              rm) + bb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (db_re * rm -
            u_im * 0.0)) + R * rm) + (eb_re * rm - v_im * 0.0);
  p_im = (((((((((((((s_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (v_re * c_m_im + p_im * c_m_re)) + (u_re * d_m_im + q_im *
    e_m_re)) + (w_re * 0.0 + 0.0 * R)) + (y_re * 0.0 + r_im * R)) + (ab_re * 0.0
    + s_im * R)) + (x_re * cb_a_re + t_im * v_a_im)) + (d * 0.0 + 0.0 * rm)) +
             (bb_re * 0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) + (db_re
            * 0.0 + u_im * rm)) + (R * 0.0 + 0.0 * rm)) + (eb_re * 0.0 + v_im *
    rm);
  s_re = cb_re * t_re - o_im * p_im;
  o_im = cb_re * p_im + o_im * t_re;
  t_re = -2.0 * d;
  u_re = -2.0 * d;
  v_re = u_re * del;
  p_im = u_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  u_re = 2.0 * (del * del);
  q_im = 2.0 * (del * 0.0 + 0.0 * del);
  w_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = 4.0 * d;
  y_re = x_re * f_m_re - 0.0 * f_m_im;
  r_im = x_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = -6.0 * del;
  ab_re = x_re * f_m_re - -0.0 * f_m_im;
  s_im = x_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  x_re = 4.0 * (m * m);
  t_im = 4.0 * (m * 0.0 + 0.0 * m);
  bb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * f_m_re - 0.0 * f_m_im;
  d_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_re = -2.0 * del;
  db_re = cb_re * f_m_re - -0.0 * f_m_im;
  u_im = cb_re * f_m_im + -0.0 * f_m_re;
  cb_re = 3.0 * (m * m);
  v_im = 3.0 * (m * 0.0 + 0.0 * m);
  eb_re = cb_re * R - v_im * 0.0;
  v_im = cb_re * 0.0 + v_im * R;
  cb_re = ((((((((((((t_re * del + 2.0 * (del * del)) + (v_re * c_m_re - p_im *
    c_m_im)) + (u_re * e_m_re - q_im * d_m_im)) + w_re * R) + (y_re * R - r_im *
    0.0)) + (ab_re * R - s_im * 0.0)) + (x_re * b_R_re - t_im * b_R_im)) + d *
               rm) + bb_re * rm) + (c_d_re * rm - d_d_im * 0.0)) + (db_re * rm -
             u_im * 0.0)) + R * rm) + (eb_re * rm - v_im * 0.0);
  p_im = (((((((((((((t_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (v_re * c_m_im + p_im * c_m_re)) + (u_re * d_m_im + q_im *
    e_m_re)) + (w_re * 0.0 + -0.0 * R)) + (y_re * 0.0 + r_im * R)) + (ab_re *
    0.0 + s_im * R)) + (x_re * b_R_im + t_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
             + (bb_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + d_d_im * rm)) +
           (db_re * 0.0 + u_im * rm)) + (R * 0.0 + 0.0 * rm)) + (eb_re * 0.0 +
    v_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  t_re = -del;
  u_re = t_re * m;
  q_im = t_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  t_re = -m;
  v_re = t_re * R;
  r_im = t_re * 0.0 + -0.0 * R;
  c_a_re = n_a_re * n_a_re - n_a_im * n_a_im;
  c_a_im = n_a_re * n_a_im + n_a_im * n_a_re;
  t_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (u_re * R - q_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (v_re * rm - r_im * 0.0));
  q_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (u_re * 0.0 + q_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (v_re * 0.0 + r_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  u_re = -del;
  v_re = u_re * m;
  r_im = u_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  u_re = -m;
  w_re = u_re * R;
  s_im = u_re * 0.0 + -0.0 * R;
  e_d_re = d * m;
  e_d_im = d * 0.0 + 0.0 * m;
  u_re = -del;
  x_re = u_re * m;
  t_im = u_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  u_re = 12.0 * ((((c_d_re * R - d_d_im * 0.0) + (v_re * R - r_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (w_re * rm - s_im * 0.0));
  r_im = 12.0 * ((((c_d_re * 0.0 + d_d_im * R) + (v_re * 0.0 + r_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (w_re * 0.0 + s_im * rm));
  c_d_re = (((e_d_re * R - e_d_im * 0.0) + (x_re * R - t_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  d_d_im = (((e_d_re * 0.0 + e_d_im * R) + (x_re * 0.0 + t_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  v_re = 2.0 * d;
  w_re = 2.0 * d;
  x_re = w_re * del;
  s_im = w_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  w_re = -2.0 * (del * del);
  t_im = -2.0 * (del * 0.0 + 0.0 * del);
  y_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = -4.0 * d;
  bb_re = ab_re * f_m_re - -0.0 * f_m_im;
  u_im = ab_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  ab_re = 6.0 * del;
  db_re = ab_re * f_m_re - 0.0 * f_m_im;
  v_im = ab_re * f_m_im + 0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  ab_re = -4.0 * (m * m);
  w_im = -4.0 * (m * 0.0 + 0.0 * m);
  eb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  fb_re = -2.0 * del;
  gb_re = fb_re * f_m_re - -0.0 * f_m_im;
  x_im = fb_re * f_m_im + -0.0 * f_m_re;
  fb_re = 3.0 * (m * m);
  y_im = 3.0 * (m * 0.0 + 0.0 * m);
  hb_re = fb_re * R - y_im * 0.0;
  y_im = fb_re * 0.0 + y_im * R;
  fb_re = -2.0 * d;
  db_a_im = -2.0 * d;
  eb_a_re = db_a_im * del;
  ab_im = db_a_im * 0.0 + -0.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_a_im = m * m;
  db_a_re = m * 0.0 + 0.0 * m;
  db_a_im = 2.0 * (del * del);
  bb_im = 2.0 * (del * 0.0 + 0.0 * del);
  eb_a_im = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ib_re = 4.0 * d;
  jb_re = ib_re * r - 0.0 * fb_a_re;
  w_a_im = ib_re * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ib_re = -6.0 * del;
  x_a_re = ib_re * r - -0.0 * fb_a_re;
  x_a_im = ib_re * fb_a_re + -0.0 * r;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  ib_re = 4.0 * (m * m);
  y_a_re = 4.0 * (m * 0.0 + 0.0 * m);
  y_a_im = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ab_a_re = d * r - 0.0 * fb_a_re;
  ab_a_im = d * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  bb_a_re = -2.0 * del;
  bb_a_im = bb_a_re * r - -0.0 * fb_a_re;
  v_a_re = bb_a_re * fb_a_re + -0.0 * r;
  bb_a_re = 3.0 * (m * m);
  w_a_re = 3.0 * (m * 0.0 + 0.0 * m);
  u_a_im = bb_a_re * R - w_a_re * 0.0;
  w_a_re = bb_a_re * 0.0 + w_a_re * R;
  bb_a_re = -3.0 * (((((((((((((v_re * del + -2.0 * (del * del)) + (x_re *
    c_m_re - s_im * c_m_im)) + (w_re * e_m_re - t_im * d_m_im)) + y_re * R) +
    (bb_re * R - u_im * 0.0)) + (db_re * R - v_im * 0.0)) + (ab_re * b_R_re -
    w_im * b_R_im)) + d * rm) + eb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) +
                      (gb_re * rm - x_im * 0.0)) + R * rm) + (hb_re * rm - y_im *
    0.0));
  s_im = -3.0 * ((((((((((((((v_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (x_re * c_m_im + s_im * c_m_re)) + (w_re * d_m_im + t_im * e_m_re))
    + (y_re * 0.0 + 0.0 * R)) + (bb_re * 0.0 + u_im * R)) + (db_re * 0.0 + v_im *
    R)) + (ab_re * b_R_im + w_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) + (eb_re *
    0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) + (gb_re * 0.0 + x_im * rm))
                  + (R * 0.0 + 0.0 * rm)) + (hb_re * 0.0 + y_im * rm));
  v_re = ((((((((((((fb_re * del + 2.0 * (del * del)) + (eb_a_re * f_m_re -
    ab_im * f_m_im)) + (db_a_im * cb_a_im - bb_im * db_a_re)) + eb_a_im * R) +
                 (jb_re * R - w_a_im * 0.0)) + (x_a_re * R - x_a_im * 0.0)) +
               (ib_re * v_a_im - y_a_re * cb_a_re)) + d * rm) + y_a_im * rm) +
            (ab_a_re * rm - ab_a_im * 0.0)) + (bb_a_im * rm - v_a_re * 0.0)) + R
          * rm) + (u_a_im * rm - w_a_re * 0.0);
  t_im = (((((((((((((fb_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (eb_a_re * f_m_im + ab_im * f_m_re)) + (db_a_im * db_a_re
    + bb_im * cb_a_im)) + (eb_a_im * 0.0 + -0.0 * R)) + (jb_re * 0.0 + w_a_im *
    R)) + (x_a_re * 0.0 + x_a_im * R)) + (ib_re * cb_a_re + y_a_re * v_a_im)) +
              (d * 0.0 + 0.0 * rm)) + (y_a_im * 0.0 + -0.0 * rm)) + (ab_a_re *
             0.0 + ab_a_im * rm)) + (bb_a_im * 0.0 + v_a_re * rm)) + (R * 0.0 +
           0.0 * rm)) + (u_a_im * 0.0 + w_a_re * rm);
  dc4.re = (36.0 * (o_a_re * o_a_re - o_a_im * o_a_im) + (u_re * c_d_re - r_im *
             d_d_im)) + (bb_a_re * v_re - s_im * t_im);
  dc4.im = (36.0 * (o_a_re * o_a_im + o_a_im * o_a_re) + (u_re * d_d_im + r_im *
             c_d_re)) + (bb_a_re * t_im + s_im * v_re);
  dc4 = c_power(dc4);
  b_d.re = -4.0 * dc4.re + (p_a_re * p_a_re - p_a_im * p_a_im);
  b_d.im = -4.0 * dc4.im + (p_a_re * p_a_im + p_a_im * p_a_re);
  dc4 = d_power(b_d);
  b_d.re = ((((-432.0 * dc3.re + (r_re * d_d_re - l_im * c_d_im)) + (p_re *
    b_a_re - n_im * b_a_im)) + (s_re * cb_re - o_im * p_im)) + (t_re * c_a_re -
             q_im * c_a_im)) + dc4.re;
  b_d.im = ((((-432.0 * dc3.im + (r_re * c_d_im + l_im * d_d_re)) + (p_re *
    b_a_im + n_im * b_a_re)) + (s_re * p_im + o_im * cb_re)) + (t_re * c_a_im +
             q_im * c_a_re)) + dc4.im;
  dc3 = e_power(b_d);
  dc4 = power(b_m);
  b_d = power(b_R);
  p_re = -0.26456684199469993 * dc4.re;
  l_im = -0.26456684199469993 * dc4.im;
  r_re = p_re * b_d.re - l_im * b_d.im;
  l_im = p_re * b_d.im + l_im * b_d.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  dc4 = power(b_d);
  p_re = r_re * dc4.re - l_im * dc4.im;
  l_im = r_re * dc4.im + l_im * dc4.re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  s_re = r_re * m;
  n_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  b_d.re = ((c_d_re * R - c_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m * b_R_re
    - 0.0 * b_R_im);
  b_d.im = ((c_d_re * 0.0 + c_d_im * R) + (s_re * 0.0 + n_im * R)) + (m * b_R_im
    + 0.0 * b_R_re);
  dc4 = c_power(b_d);
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  s_re = r_re * m;
  n_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  d_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  t_re = r_re * m;
  o_im = r_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  r_re = -m;
  u_re = r_re * R;
  p_im = r_re * 0.0 + -0.0 * R;
  r_re = 432.0 * (((c_d_re * R - c_d_im * 0.0) + (s_re * R - n_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  n_im = 432.0 * (((c_d_re * 0.0 + c_d_im * R) + (s_re * 0.0 + n_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  c_d_re = (((d_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (u_re * rm - p_im * 0.0);
  c_d_im = (((d_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (u_re * 0.0 + p_im * rm);
  s_re = r_re * c_d_re - n_im * c_d_im;
  n_im = r_re * c_d_im + n_im * c_d_re;
  c_d_re = d * m;
  c_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  t_re = r_re * m;
  o_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  d_d_re = (((c_d_re * R - c_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
             b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0);
  c_d_im = (((c_d_re * 0.0 + c_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
             b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  r_re = -del;
  t_re = r_re * m;
  o_im = r_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  b_a_re = q_a_re * q_a_re - q_a_im * q_a_im;
  b_a_im = q_a_re * q_a_im + q_a_im * q_a_re;
  r_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (t_re * R - o_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (c_m_re * rm - c_m_im * 0.0));
  o_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (t_re * 0.0 + o_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (c_m_re * 0.0 + c_m_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  t_re = -del;
  u_re = t_re * m;
  p_im = t_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  t_re = 2.0 * d;
  v_re = 2.0 * d;
  w_re = v_re * del;
  q_im = v_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  v_re = -2.0 * (del * del);
  r_im = -2.0 * (del * 0.0 + 0.0 * del);
  x_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = -4.0 * d;
  ab_re = y_re * f_m_re - -0.0 * f_m_im;
  s_im = y_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = 6.0 * del;
  bb_re = y_re * f_m_re - 0.0 * f_m_im;
  t_im = y_re * f_m_im + 0.0 * f_m_re;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  y_re = -4.0 * (m * m);
  u_im = -4.0 * (m * 0.0 + 0.0 * m);
  cb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  db_re = -2.0 * del;
  eb_re = db_re * f_m_re - -0.0 * f_m_im;
  v_im = db_re * f_m_im + -0.0 * f_m_re;
  db_re = 3.0 * (m * m);
  w_im = 3.0 * (m * 0.0 + 0.0 * m);
  fb_re = db_re * R - w_im * 0.0;
  w_im = db_re * 0.0 + w_im * R;
  db_re = 54.0 * (((c_d_re * R - d_d_im * 0.0) + (u_re * R - p_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im));
  p_im = 54.0 * (((c_d_re * 0.0 + d_d_im * R) + (u_re * 0.0 + p_im * R)) + (m *
    b_R_im + 0.0 * b_R_re));
  u_re = ((((((((((((t_re * del + -2.0 * (del * del)) + (w_re * c_m_re - q_im *
    c_m_im)) + (v_re * e_m_re - r_im * d_m_im)) + x_re * R) + (ab_re * R - s_im *
    0.0)) + (bb_re * R - t_im * 0.0)) + (y_re * v_a_im - u_im * cb_a_re)) + d *
              rm) + cb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) + (eb_re * rm -
            v_im * 0.0)) + R * rm) + (fb_re * rm - w_im * 0.0);
  q_im = (((((((((((((t_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (w_re * c_m_im + q_im * c_m_re)) + (v_re * d_m_im + r_im *
    e_m_re)) + (x_re * 0.0 + 0.0 * R)) + (ab_re * 0.0 + s_im * R)) + (bb_re *
    0.0 + t_im * R)) + (y_re * cb_a_re + u_im * v_a_im)) + (d * 0.0 + 0.0 * rm))
             + (cb_re * 0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) +
           (eb_re * 0.0 + v_im * rm)) + (R * 0.0 + 0.0 * rm)) + (fb_re * 0.0 +
    w_im * rm);
  t_re = db_re * u_re - p_im * q_im;
  p_im = db_re * q_im + p_im * u_re;
  u_re = -2.0 * d;
  v_re = -2.0 * d;
  w_re = v_re * del;
  q_im = v_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  v_re = 2.0 * (del * del);
  r_im = 2.0 * (del * 0.0 + 0.0 * del);
  x_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = 4.0 * d;
  ab_re = y_re * f_m_re - 0.0 * f_m_im;
  s_im = y_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = -6.0 * del;
  bb_re = y_re * f_m_re - -0.0 * f_m_im;
  t_im = y_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  y_re = 4.0 * (m * m);
  u_im = 4.0 * (m * 0.0 + 0.0 * m);
  cb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * f_m_re - 0.0 * f_m_im;
  d_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  db_re = -2.0 * del;
  eb_re = db_re * f_m_re - -0.0 * f_m_im;
  v_im = db_re * f_m_im + -0.0 * f_m_re;
  db_re = 3.0 * (m * m);
  w_im = 3.0 * (m * 0.0 + 0.0 * m);
  fb_re = db_re * R - w_im * 0.0;
  w_im = db_re * 0.0 + w_im * R;
  db_re = ((((((((((((u_re * del + 2.0 * (del * del)) + (w_re * c_m_re - q_im *
    c_m_im)) + (v_re * e_m_re - r_im * d_m_im)) + x_re * R) + (ab_re * R - s_im *
    0.0)) + (bb_re * R - t_im * 0.0)) + (y_re * b_R_re - u_im * b_R_im)) + d *
               rm) + cb_re * rm) + (c_d_re * rm - d_d_im * 0.0)) + (eb_re * rm -
             v_im * 0.0)) + R * rm) + (fb_re * rm - w_im * 0.0);
  q_im = (((((((((((((u_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (w_re * c_m_im + q_im * c_m_re)) + (v_re * d_m_im + r_im *
    e_m_re)) + (x_re * 0.0 + -0.0 * R)) + (ab_re * 0.0 + s_im * R)) + (bb_re *
    0.0 + t_im * R)) + (y_re * b_R_im + u_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
             + (cb_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + d_d_im * rm)) +
           (eb_re * 0.0 + v_im * rm)) + (R * 0.0 + 0.0 * rm)) + (fb_re * 0.0 +
    w_im * rm);
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  u_re = -del;
  v_re = u_re * m;
  r_im = u_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  u_re = -m;
  w_re = u_re * R;
  s_im = u_re * 0.0 + -0.0 * R;
  c_a_re = r_a_re * r_a_re - r_a_im * r_a_im;
  c_a_im = r_a_re * r_a_im + r_a_im * r_a_re;
  u_re = 27.0 * ((((c_d_re * R - d_d_im * 0.0) + (v_re * R - r_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (w_re * rm - s_im * 0.0));
  r_im = 27.0 * ((((c_d_re * 0.0 + d_d_im * R) + (v_re * 0.0 + r_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (w_re * 0.0 + s_im * rm));
  c_d_re = d * m;
  d_d_im = d * 0.0 + 0.0 * m;
  v_re = -del;
  w_re = v_re * m;
  s_im = v_re * 0.0 + -0.0 * m;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  v_re = -m;
  x_re = v_re * R;
  t_im = v_re * 0.0 + -0.0 * R;
  e_d_re = d * m;
  e_d_im = d * 0.0 + 0.0 * m;
  v_re = -del;
  y_re = v_re * m;
  u_im = v_re * 0.0 + -0.0 * m;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  c_m_re = m * R;
  c_m_im = m * 0.0 + 0.0 * R;
  v_re = 12.0 * ((((c_d_re * R - d_d_im * 0.0) + (w_re * R - s_im * 0.0)) + (m *
    b_R_re - 0.0 * b_R_im)) + (x_re * rm - t_im * 0.0));
  s_im = 12.0 * ((((c_d_re * 0.0 + d_d_im * R) + (w_re * 0.0 + s_im * R)) + (m *
    b_R_im + 0.0 * b_R_re)) + (x_re * 0.0 + t_im * rm));
  c_d_re = (((e_d_re * R - e_d_im * 0.0) + (y_re * R - u_im * 0.0)) + (m *
             v_a_im - 0.0 * cb_a_re)) + (c_m_re * rm - c_m_im * 0.0);
  d_d_im = (((e_d_re * 0.0 + e_d_im * R) + (y_re * 0.0 + u_im * R)) + (m *
             cb_a_re + 0.0 * v_a_im)) + (c_m_re * 0.0 + c_m_im * rm);
  w_re = 2.0 * d;
  x_re = 2.0 * d;
  y_re = x_re * del;
  t_im = x_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  x_re = -2.0 * (del * del);
  u_im = -2.0 * (del * 0.0 + 0.0 * del);
  ab_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  bb_re = -4.0 * d;
  cb_re = bb_re * f_m_re - -0.0 * f_m_im;
  v_im = bb_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  bb_re = 6.0 * del;
  eb_re = bb_re * f_m_re - 0.0 * f_m_im;
  w_im = bb_re * f_m_im + 0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  bb_re = -4.0 * (m * m);
  x_im = -4.0 * (m * 0.0 + 0.0 * m);
  fb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  e_d_re = d * f_m_re - 0.0 * f_m_im;
  e_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  gb_re = -2.0 * del;
  hb_re = gb_re * f_m_re - -0.0 * f_m_im;
  y_im = gb_re * f_m_im + -0.0 * f_m_re;
  gb_re = 3.0 * (m * m);
  ab_im = 3.0 * (m * 0.0 + 0.0 * m);
  db_a_im = gb_re * R - ab_im * 0.0;
  ab_im = gb_re * 0.0 + ab_im * R;
  gb_re = -2.0 * d;
  eb_a_re = -2.0 * d;
  eb_a_im = eb_a_re * del;
  bb_im = eb_a_re * 0.0 + -0.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  cb_a_im = m * m;
  db_a_re = m * 0.0 + 0.0 * m;
  eb_a_re = 2.0 * (del * del);
  w_a_im = 2.0 * (del * 0.0 + 0.0 * del);
  ib_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  jb_re = 4.0 * d;
  x_a_re = jb_re * r - 0.0 * fb_a_re;
  x_a_im = jb_re * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  jb_re = -6.0 * del;
  y_a_im = jb_re * r - -0.0 * fb_a_re;
  y_a_re = jb_re * fb_a_re + -0.0 * r;
  v_a_im = R * R;
  cb_a_re = R * 0.0 + 0.0 * R;
  jb_re = 4.0 * (m * m);
  v_a_re = 4.0 * (m * 0.0 + 0.0 * m);
  bb_a_re = -2.0 * del;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  ab_a_re = d * r - 0.0 * fb_a_re;
  ab_a_im = d * fb_a_re + 0.0 * r;
  r = m * m;
  fb_a_re = m * 0.0 + 0.0 * m;
  bb_a_im = -2.0 * del;
  u_a_im = bb_a_im * r - -0.0 * fb_a_re;
  w_a_re = bb_a_im * fb_a_re + -0.0 * r;
  bb_a_im = 3.0 * (m * m);
  r = 3.0 * (m * 0.0 + 0.0 * m);
  fb_a_re = bb_a_im * R - r * 0.0;
  r = bb_a_im * 0.0 + r * R;
  bb_a_im = -3.0 * (((((((((((((w_re * del + -2.0 * (del * del)) + (y_re *
    c_m_re - t_im * c_m_im)) + (x_re * e_m_re - u_im * d_m_im)) + ab_re * R) +
    (cb_re * R - v_im * 0.0)) + (eb_re * R - w_im * 0.0)) + (bb_re * b_R_re -
    x_im * b_R_im)) + d * rm) + fb_re * rm) + (e_d_re * rm - e_d_im * 0.0)) +
                      (hb_re * rm - y_im * 0.0)) + R * rm) + (db_a_im * rm -
    ab_im * 0.0));
  t_im = -3.0 * ((((((((((((((w_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 *
    del)) + (y_re * c_m_im + t_im * c_m_re)) + (x_re * d_m_im + u_im * e_m_re))
    + (ab_re * 0.0 + 0.0 * R)) + (cb_re * 0.0 + v_im * R)) + (eb_re * 0.0 + w_im
    * R)) + (bb_re * b_R_im + x_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) + (fb_re *
    0.0 + -0.0 * rm)) + (e_d_re * 0.0 + e_d_im * rm)) + (hb_re * 0.0 + y_im * rm))
                  + (R * 0.0 + 0.0 * rm)) + (db_a_im * 0.0 + ab_im * rm));
  w_re = ((((((((((((gb_re * del + 2.0 * (del * del)) + (eb_a_im * f_m_re -
    bb_im * f_m_im)) + (eb_a_re * cb_a_im - w_a_im * db_a_re)) + ib_re * R) +
                 (x_a_re * R - x_a_im * 0.0)) + (y_a_im * R - y_a_re * 0.0)) +
               (jb_re * v_a_im - v_a_re * cb_a_re)) + d * rm) + bb_a_re * rm) +
            (ab_a_re * rm - ab_a_im * 0.0)) + (u_a_im * rm - w_a_re * 0.0)) + R *
          rm) + (fb_a_re * rm - r * 0.0);
  u_im = (((((((((((((gb_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (eb_a_im * f_m_im + bb_im * f_m_re)) + (eb_a_re * db_a_re
    + w_a_im * cb_a_im)) + (ib_re * 0.0 + -0.0 * R)) + (x_a_re * 0.0 + x_a_im *
    R)) + (y_a_im * 0.0 + y_a_re * R)) + (jb_re * cb_a_re + v_a_re * v_a_im)) +
              (d * 0.0 + 0.0 * rm)) + (bb_a_re * 0.0 + -0.0 * rm)) + (ab_a_re *
             0.0 + ab_a_im * rm)) + (u_a_im * 0.0 + w_a_re * rm)) + (R * 0.0 +
           0.0 * rm)) + (fb_a_re * 0.0 + r * rm);
  b_d.re = (36.0 * (s_a_re * s_a_re - s_a_im * s_a_im) + (v_re * c_d_re - s_im *
             d_d_im)) + (bb_a_im * w_re - t_im * u_im);
  b_d.im = (36.0 * (s_a_re * s_a_im + s_a_im * s_a_re) + (v_re * d_d_im + s_im *
             c_d_re)) + (bb_a_im * u_im + t_im * w_re);
  b_d = c_power(b_d);
  dc6.re = -4.0 * b_d.re + (t_a_re * t_a_re - t_a_im * t_a_im);
  dc6.im = -4.0 * b_d.im + (t_a_re * t_a_im + t_a_im * t_a_re);
  b_d = d_power(dc6);
  dc6.re = ((((-432.0 * dc4.re + (s_re * d_d_re - n_im * c_d_im)) + (r_re *
    b_a_re - o_im * b_a_im)) + (t_re * db_re - p_im * q_im)) + (u_re * c_a_re -
             r_im * c_a_im)) + b_d.re;
  dc6.im = ((((-432.0 * dc4.im + (s_re * c_d_im + n_im * d_d_re)) + (r_re *
    b_a_im + o_im * b_a_re)) + (t_re * q_im + p_im * db_re)) + (u_re * c_a_im +
             r_im * c_a_re)) + b_d.im;
  dc4 = f_power(dc6);
  b_d = power(b_m);
  dc6 = power(b_R);
  r_re = -8.0 * b_d.re;
  n_im = -8.0 * b_d.im;
  s_re = r_re * dc6.re - n_im * dc6.im;
  n_im = r_re * dc6.im + n_im * dc6.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  b_d = power(b_d);
  r_re = s_re * b_d.re - n_im * b_d.im;
  n_im = s_re * b_d.im + n_im * b_d.re;
  s_re = 2.0 * d;
  t_re = 2.0 * d;
  u_re = t_re * del;
  o_im = t_re * 0.0 + 0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  t_re = -2.0 * (del * del);
  p_im = -2.0 * (del * 0.0 + 0.0 * del);
  v_re = 2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  w_re = -4.0 * d;
  x_re = w_re * f_m_re - -0.0 * f_m_im;
  q_im = w_re * f_m_im + -0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  w_re = 6.0 * del;
  y_re = w_re * f_m_re - 0.0 * f_m_im;
  r_im = w_re * f_m_im + 0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  w_re = -4.0 * (m * m);
  s_im = -4.0 * (m * 0.0 + 0.0 * m);
  ab_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * f_m_re - 0.0 * f_m_im;
  c_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  bb_re = -2.0 * del;
  cb_re = bb_re * f_m_re - -0.0 * f_m_im;
  t_im = bb_re * f_m_im + -0.0 * f_m_re;
  bb_re = 3.0 * (m * m);
  u_im = 3.0 * (m * 0.0 + 0.0 * m);
  db_re = bb_re * R - u_im * 0.0;
  u_im = bb_re * 0.0 + u_im * R;
  bb_re = ((((((((((((s_re * del + -2.0 * (del * del)) + (u_re * c_m_re - o_im *
    c_m_im)) + (t_re * e_m_re - p_im * d_m_im)) + v_re * R) + (x_re * R - q_im *
    0.0)) + (y_re * R - r_im * 0.0)) + (w_re * b_R_re - s_im * b_R_im)) + d * rm)
              + ab_re * rm) + (c_d_re * rm - c_d_im * 0.0)) + (cb_re * rm - t_im
             * 0.0)) + R * rm) + (db_re * rm - u_im * 0.0);
  o_im = (((((((((((((s_re * 0.0 + 0.0 * del) + -2.0 * (del * 0.0 + 0.0 * del))
                    + (u_re * c_m_im + o_im * c_m_re)) + (t_re * d_m_im + p_im *
    e_m_re)) + (v_re * 0.0 + 0.0 * R)) + (x_re * 0.0 + q_im * R)) + (y_re * 0.0
    + r_im * R)) + (w_re * b_R_im + s_im * b_R_re)) + (d * 0.0 + 0.0 * rm)) +
             (ab_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + c_d_im * rm)) + (cb_re
            * 0.0 + t_im * rm)) + (R * 0.0 + 0.0 * rm)) + (db_re * 0.0 + u_im *
    rm);
  b_d = power(b_m);
  dc6 = power(b_R);
  s_re = -24.0 * b_d.re;
  p_im = -24.0 * b_d.im;
  t_re = s_re * dc6.re - p_im * dc6.im;
  p_im = s_re * dc6.im + p_im * dc6.re;
  c_d_re = (d + -del) + R;
  s_re = t_re * c_d_re - p_im * 0.0;
  p_im = t_re * 0.0 + p_im * c_d_re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  b_d = b_power(b_d);
  t_re = s_re * b_d.re - p_im * b_d.im;
  p_im = s_re * b_d.im + p_im * b_d.re;
  s_re = -2.0 * d;
  u_re = -2.0 * d;
  v_re = u_re * del;
  q_im = u_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  u_re = 2.0 * (del * del);
  r_im = 2.0 * (del * 0.0 + 0.0 * del);
  w_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = 4.0 * d;
  y_re = x_re * f_m_re - 0.0 * f_m_im;
  s_im = x_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  x_re = -6.0 * del;
  ab_re = x_re * f_m_re - -0.0 * f_m_im;
  t_im = x_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  x_re = 4.0 * (m * m);
  u_im = 4.0 * (m * 0.0 + 0.0 * m);
  cb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * f_m_re - 0.0 * f_m_im;
  c_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  db_re = -2.0 * del;
  eb_re = db_re * f_m_re - -0.0 * f_m_im;
  v_im = db_re * f_m_im + -0.0 * f_m_re;
  db_re = 3.0 * (m * m);
  w_im = 3.0 * (m * 0.0 + 0.0 * m);
  fb_re = db_re * R - w_im * 0.0;
  w_im = db_re * 0.0 + w_im * R;
  db_re = ((((((((((((s_re * del + 2.0 * (del * del)) + (v_re * c_m_re - q_im *
    c_m_im)) + (u_re * e_m_re - r_im * d_m_im)) + w_re * R) + (y_re * R - s_im *
    0.0)) + (ab_re * R - t_im * 0.0)) + (x_re * b_R_re - u_im * b_R_im)) + d *
               rm) + cb_re * rm) + (c_d_re * rm - c_d_im * 0.0)) + (eb_re * rm -
             v_im * 0.0)) + R * rm) + (fb_re * rm - w_im * 0.0);
  q_im = (((((((((((((s_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                    + (v_re * c_m_im + q_im * c_m_re)) + (u_re * d_m_im + r_im *
    e_m_re)) + (w_re * 0.0 + -0.0 * R)) + (y_re * 0.0 + s_im * R)) + (ab_re *
    0.0 + t_im * R)) + (x_re * b_R_im + u_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
             + (cb_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + c_d_im * rm)) +
           (eb_re * 0.0 + v_im * rm)) + (R * 0.0 + 0.0 * rm)) + (fb_re * 0.0 +
    w_im * rm);
  b_d = g_power(b_m);
  dc6 = g_power(b_R);
  s_re = -b_d.re;
  r_im = -b_d.im;
  u_re = s_re * dc6.re - r_im * dc6.im;
  r_im = s_re * dc6.im + r_im * dc6.re;
  b_d.re = ((d + -del) + R) + rm;
  b_d.im = 0.0;
  b_d = g_power(b_d);
  s_re = u_re * b_d.re - r_im * b_d.im;
  r_im = u_re * b_d.im + r_im * b_d.re;
  u_re = -2.0 * d;
  v_re = -2.0 * d;
  w_re = v_re * del;
  s_im = v_re * 0.0 + -0.0 * del;
  c_m_re = m * m;
  c_m_im = m * 0.0 + 0.0 * m;
  e_m_re = m * m;
  d_m_im = m * 0.0 + 0.0 * m;
  v_re = 2.0 * (del * del);
  t_im = 2.0 * (del * 0.0 + 0.0 * del);
  x_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = 4.0 * d;
  ab_re = y_re * f_m_re - 0.0 * f_m_im;
  u_im = y_re * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  y_re = -6.0 * del;
  cb_re = y_re * f_m_re - -0.0 * f_m_im;
  v_im = y_re * f_m_im + -0.0 * f_m_re;
  b_R_re = R * R;
  b_R_im = R * 0.0 + 0.0 * R;
  y_re = 4.0 * (m * m);
  w_im = 4.0 * (m * 0.0 + 0.0 * m);
  eb_re = -2.0 * del;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  c_d_re = d * f_m_re - 0.0 * f_m_im;
  c_d_im = d * f_m_im + 0.0 * f_m_re;
  f_m_re = m * m;
  f_m_im = m * 0.0 + 0.0 * m;
  fb_re = -2.0 * del;
  gb_re = fb_re * f_m_re - -0.0 * f_m_im;
  x_im = fb_re * f_m_im + -0.0 * f_m_re;
  fb_re = 3.0 * (m * m);
  y_im = 3.0 * (m * 0.0 + 0.0 * m);
  hb_re = fb_re * R - y_im * 0.0;
  y_im = fb_re * 0.0 + y_im * R;
  b_d.re = ((((((((((((u_re * del + 2.0 * (del * del)) + (w_re * c_m_re - s_im *
    c_m_im)) + (v_re * e_m_re - t_im * d_m_im)) + x_re * R) + (ab_re * R - u_im *
    0.0)) + (cb_re * R - v_im * 0.0)) + (y_re * b_R_re - w_im * b_R_im)) + d *
                rm) + eb_re * rm) + (c_d_re * rm - c_d_im * 0.0)) + (gb_re * rm
              - x_im * 0.0)) + R * rm) + (hb_re * rm - y_im * 0.0);
  b_d.im = (((((((((((((u_re * 0.0 + -0.0 * del) + 2.0 * (del * 0.0 + 0.0 * del))
                      + (w_re * c_m_im + s_im * c_m_re)) + (v_re * d_m_im + t_im
    * e_m_re)) + (x_re * 0.0 + -0.0 * R)) + (ab_re * 0.0 + u_im * R)) + (cb_re *
    0.0 + v_im * R)) + (y_re * b_R_im + w_im * b_R_re)) + (d * 0.0 + 0.0 * rm))
               + (eb_re * 0.0 + -0.0 * rm)) + (c_d_re * 0.0 + c_d_im * rm)) +
             (gb_re * 0.0 + x_im * rm)) + (R * 0.0 + 0.0 * rm)) + (hb_re * 0.0 +
    y_im * rm);
  b_d = c_power(b_d);
  u_re = 0.25 * (((r_re * bb_re - n_im * o_im) + (t_re * db_re - p_im * q_im)) +
                 (s_re * b_d.re - r_im * b_d.im));
  n_im = 0.25 * (((r_re * o_im + n_im * bb_re) + (t_re * q_im + p_im * db_re)) +
                 (s_re * b_d.im + r_im * b_d.re));
  b_d.re = (((((j_re * dc1.re - 0.0 * dc1.im) + (n_re * dc2.re - i_im * dc2.im))
              + (o_re * a_re - j_im * a_im)) + (q_re * dc3.re - k_im * dc3.im))
            + (p_re * dc4.re - l_im * dc4.im)) + (u_re * a.re - n_im * a.im);
  b_d.im = (((((j_re * dc1.im + 0.0 * dc1.re) + (n_re * dc2.im + i_im * dc2.re))
              + (o_re * a_im + j_im * a_re)) + (q_re * dc3.im + k_im * dc3.re))
            + (p_re * dc4.im + l_im * dc4.re)) + (u_re * a.im + n_im * a.re);
  return (((b_re * dc0.re - im * dc0.im) * (((((((((((((re * del + 2.0 * (del *
    del)) + (d_re * m_re - b_im * m_im)) + (c_re * b_m_re - c_im * b_m_im)) +
    e_re * R) + (g_re * R - d_im * 0.0)) + (h_re * R - e_im * 0.0)) + (f_re *
    R_re - f_im * R_im)) + d * rm) + i_re * rm) + (b_d_re * rm - b_d_im * 0.0))
              + (k_re * rm - g_im * 0.0)) + R * rm) + (l_re * rm - h_im * 0.0))
           - (b_re * dc0.im + im * dc0.re) * ((((((((((((((re * 0.0 + -0.0 * del)
    + 2.0 * (del * 0.0 + 0.0 * del)) + (d_re * m_im + b_im * m_re)) + (c_re *
    b_m_im + c_im * b_m_re)) + (e_re * 0.0 + -0.0 * R)) + (g_re * 0.0 + d_im * R))
    + (h_re * 0.0 + e_im * R)) + (f_re * R_im + f_im * R_re)) + (d * 0.0 + 0.0 *
    rm)) + (i_re * 0.0 + -0.0 * rm)) + (b_d_re * 0.0 + b_d_im * rm)) + (k_re *
    0.0 + g_im * rm)) + (R * 0.0 + 0.0 * rm)) + (l_re * 0.0 + h_im * rm))) + 0.5
          * (d_power(dc5)).re) + -0.5 * (d_power(b_d)).re;
}

/*
 * File trailer for mTot.c
 *
 * [EOF]
 */
