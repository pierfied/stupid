//
// Created by pierfied on 11/23/18.
//

#include "cosmology.h"

#include <math.h>
#include <gsl/gsl_integration.h>

double cosmology::f(double a) {
    return 1 / sqrt((Omega_m0 + Omega_k0 * a + Omega_l0 * a * a * a) / a);
}

double adot(double a, double Omega_m0, double Omega_l0, double Omega_k0) {
    return a * sqrt(Omega_m0 / (a * a * a) + Omega_l0 + Omega_k0 / (a * a));
}

double D_integrand(double a_prime, void *params) {
    double *p = (double *) params;

    double Omega_m0 = p[0];
    double Omega_l0 = p[1];
    double Omega_k0 = p[2];

    return pow(adot(a_prime, Omega_m0, Omega_l0, Omega_k0), -3);
}

double cosmology::unnormalized_D(double a) {
    double eps = 1.49e-08;
    long limit = 50;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);

    double result, error;
    double params[3] = {Omega_m0, Omega_l0, Omega_k0};

    gsl_function F;
    F.function = &D_integrand;
    F.params = params;

    gsl_integration_qags(&F, 0, a, eps, eps, limit, w, &result, &error);

    gsl_integration_workspace_free(w);

    double uD = 2.5 * Omega_m0 * adot(a, Omega_m0, Omega_l0, Omega_k0) / a * result;

    return uD;
}

double cosmology::D(double a) {
    return unnormalized_D(a) / unnormalized_D(1);
}

double cosmology::Ddot(double a) {
    double unnormalized_D_dot = Omega_m0 / (2 * a * adot(a, Omega_m0, Omega_l0, Omega_k0)) *
                                (5 - 3 * unnormalized_D(a) / a - 2 * Omega_k0 * unnormalized_D(a) / Omega_m0);

    return unnormalized_D_dot / unnormalized_D(1);
}

void cosmology::load_P(double *k, double *P, long size) {
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(spline, k, P, size);
}
