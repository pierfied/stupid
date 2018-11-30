//
// Created by pierfied on 11/23/18.
//

#include "cosmology.h"

#include <math.h>

double cosmology::f(double a) {
    return 1 / sqrt((Omega_m0 + Omega_k0 * a + Omega_l0 * a * a * a) / a);
}

double cosmology::lgf(double a) {
    return 0;
}