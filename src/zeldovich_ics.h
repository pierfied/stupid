//
// Created by rafa on 12/1/18.
//

#ifndef STUPID_ZELDOVICH_ICS_H
#define STUPID_ZELDOVICH_ICS_H

#include "initial_conditions.h"
#include "fftw3.h"

class zeldovich_ics: public initial_conditions {
    cosmology cosmo;

    double a0;
    double delta_a;

    const int Np1;
    const int Ng1;

    array_3d<double> real_Sx, real_Sy, real_Sz;
    array_3d<fftw_complex> fourier_Sx, fourier_Sy, fourier_Sz;

    fftw_plan backward_planx, backward_plany, backward_planz;

public:
    zeldovich_ics(cosmology cosmo, particle_list &plist, grid &g, double a0, double delta_a) :
        initial_conditions(plist, g), cosmo(cosmo), a0(a0), delta_a(delta_a), Np1(int(cbrt(plist.num_particles))), Ng1(g.n),
        real_Sx(Np1, Np1, Np1), real_Sy(Np1, Np1, Np1), real_Sz(Np1, Np1, Np1),
        fourier_Sx(Np1, Np1, Np1/2+1), fourier_Sy(Np1, Np1, Np1/2+1), fourier_Sz(Np1, Np1, Np1/2+1),
        backward_planx(fftw_plan_dft_c2r_3d(Np1, Np1, Np1, fourier_Sx.data, real_Sx.data, FFTW_MEASURE)),
        backward_plany(fftw_plan_dft_c2r_3d(Np1, Np1, Np1, fourier_Sy.data, real_Sy.data, FFTW_MEASURE)),
        backward_planz(fftw_plan_dft_c2r_3d(Np1, Np1, Np1, fourier_Sz.data, real_Sz.data, FFTW_MEASURE)) {}

    ~zeldovich_ics() {
        fftw_destroy_plan(backward_planx);
        fftw_destroy_plan(backward_plany);
        fftw_destroy_plan(backward_planz);
    }

private:
    void generate_uniform_distribution();

    void fourier_displacement_vec();

    void real_displacement_vec();

    void apply_ZA();

public:
    void setup_ics() override;
};

#endif //STUPID_ZELDOVICH_ICS_H
