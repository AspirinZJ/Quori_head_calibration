#include "rt_nonfinite.h"
#include "sphereToimage.h"
#include "sphereToplane.h"

void sphereToimage(double theta, double psi, double H, double W, double *px,
                   double *py, double dx, double dy)
{
  double xp1;
  double yp1;
  double ax;
  double ay;
  double xp;
  double yp;
  double s;

  sphereToplane(0.52359877559829882, 3.1415926535897931, 4.0, 1.5, 2.0, 1.0,
                0.059068559067049511, 8.4476252817457738, &xp1, &yp1);
  sphereToplane(1.5707963267948966, 2.0943951023931953, 4.0, 1.5, 2.0, 1.0,
                0.059068559067049511, 8.4476252817457738, &ax, &ay);
  sphereToplane(theta, psi, 4.0, 1.5, 2.0, 1.0, 0.059068559067049511,
                8.4476252817457738, &xp, &yp);

  s = 10.8/7.5;
  ax = -159.78000000000009 / (xp1 - ax);
  ay = -248.13 / (yp1 - ay);

  *px = s * (ax * xp + (898.62 - ax * xp1)) + dx * W / 13.333333333333333;
  *py = s * (ay * yp + (182.37 - ay * yp1)) + dy * H / 7.5;
}
