//
// Created by pierfied on 11/23/18.
//

#ifndef STUPID_COSMOLOGY_H
#define STUPID_COSMOLOGY_H


class cosmology {
public:
    const double Omega_0;
    const double Omega_m0;
    const double Omega_k0;
    const double Omega_l0;
    const double H_0;

    cosmology(double Omega_m0, double Omega_k0, double Omega_l0, double H_0) : Omega_m0(Omega_m0), Omega_k0(Omega_k0),
                                                                               Omega_l0(Omega_l0), H_0(H_0),
                                                                               Omega_0(Omega_m0 + Omega_l0) {}

    double f(double a);
};


#endif //STUPID_COSMOLOGY_H
