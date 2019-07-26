#include "rt_nonfinite.h"
#include "rotation_theta.h"
#include "rotation_theta_initialize.h"

void rotation_theta_initialize(void)
{
  rt_InitInfAndNaN(8U);
}