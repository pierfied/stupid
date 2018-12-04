//
// Created by pierfied on 11/23/18.
//

#ifndef STUPID_COSMOLOGY_H
#define STUPID_COSMOLOGY_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


class cosmology {
public:
    const double Omega_0;
    const double Omega_m0;
    const double Omega_k0;
    const double Omega_l0;

    gsl_interp_accel *acc;
    gsl_spline *spline;


    cosmology(double Omega_m0, double Omega_k0, double Omega_l0) :
            Omega_m0(Omega_m0), Omega_k0(Omega_k0), Omega_l0(Omega_l0), Omega_0(Omega_m0 + Omega_l0) {}

    double f(double a);

    double lgf(double a);

    void load_P(double *k, double *P, int size);

    inline double P_at_k(double k) {
        return gsl_spline_eval(spline, k, acc);
    }
};


#endif //STUPID_COSMOLOGY_H
