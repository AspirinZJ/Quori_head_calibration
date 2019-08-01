#ifndef SPHERETOPLANE_H
#define SPHERETOPLANE_H

// #include <stddef.h>
// #include <stdlib.h>
// #include "tmwtypes.h"
#include <math.h>
// #include "rt_nonfinite.h"
#include "sphereToimage.h"
#include "mTot.h"

extern void sphereToplane(double theta, double psi, double R, double rm, double ro, double h, double epsilon, double L, double *xp, double *yp);

#endif
