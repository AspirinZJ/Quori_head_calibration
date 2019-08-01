#ifndef SPHERETOIMAGE_H
#define SPHERETOIMAGE_H

#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"

extern void sphereToimage(float R, double theta, double psi, double H, double W, double *px, double *py, double dx, double dy);

#endif
