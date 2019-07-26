/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: sphereToplane.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 09-Aug-2018 12:35:13
 */

#include <math.h>
#include "rt_nonfinite.h"
#include "sphereToimage.h"
#include "sphereToplane.h"
#include "mTot.h"

void sphereToplane(double theta, double psi, double R, double rm, double ro,
                   double h, double epsilon, double L, double *xp, double *yp)
{
    double d;
    double del;
    double x;
    double phi;
    double l1;
    double l2;
    d = h + sqrt(R * R - ro * ro);
    del = R - sqrt(R * R - rm * rm);
    x = mTot(tan(theta / 2.0), d, rm, R, del);
    phi = 2.0 * atan(x);
    phi = atan(sin(phi) / (((R - del) + d) / rm - cos(phi)));
    l1 = L * tan(epsilon - phi);
    l2 = L * tan(epsilon + phi);
    del = (-l1 - (-l2)) / 2.0;
    phi = del * (sin(epsilon) / cos(phi));
    phi = sqrt(del * del - phi * phi);
    x = cos(psi - 1.5707963267948966);
    d = sin(psi - 1.5707963267948966);
    phi = sqrt(1.0 / (x * x / (del * del) + d * d / (phi * phi)));
    *xp = (-l1 + -l2) / 2.0 + phi * cos(psi - 1.5707963267948966);
    *yp = phi * sin(psi - 1.5707963267948966);
}