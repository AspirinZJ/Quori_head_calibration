#include <math.h>
#include "rt_nonfinite.h"
#include "rotation_theta.h"

void rotation_theta(double theta0, double psi0, double psiO, double dtheta,
                    double *theta, double *psi)
{
    double psi_prime;
    double dv0[9];
    int i0;
    static const signed char iv0[3] = {0, 1, 0};

    double dv1[3];
    double r[3];
    int i1;
    psi_prime = psi0 - psiO;
    dv0[0] = cos(dtheta);
    dv0[3] = 0.0;
    dv0[6] = sin(dtheta);
    for (i0 = 0; i0 < 3; i0++)
    {
        dv0[1 + 3 * i0] = iv0[i0];
    }

    dv0[2] = -sin(dtheta);
    dv0[5] = 0.0;
    dv0[8] = cos(dtheta);
    dv1[0] = sin(theta0) * cos(psi_prime);
    dv1[1] = sin(theta0) * sin(psi_prime);
    dv1[2] = cos(theta0);
    for (i0 = 0; i0 < 3; i0++)
    {
        r[i0] = 0.0;
        for (i1 = 0; i1 < 3; i1++)
        {
            r[i0] += dv0[i0 + 3 * i1] * dv1[i1];
        }
    }

    *theta = acos(r[2]);
    *psi = psiO + atan(r[1] / r[0]);
}