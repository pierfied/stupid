//
// Created by pierfied on 11/23/18.
//

#include "cosmology.h"

#include <math.h>

double cosmology::f(double a) {
    return 1 / sqrt((Omega_m0 + Omega_k0 * a + Omega_l0 * a * a * a) / a);
}

double cosmology::lgf(double a) {
    // From Carroll 1992 eq 29
    return 0.5*5*a*Omega_m0/(pow(Omega_m0, 4./7.) - Omega_l0 + (1+0.5*Omega_m0)*(1+Omega_l0/70));
}

void cosmology::load_P(double *k, double *P, int size) {
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(spline, k, P, size);
}