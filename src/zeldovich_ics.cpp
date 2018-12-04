//
// Created by rafa on 12/1/18.
//

#include "zeldovich_ics.h"
#include <random>

void zeldovich_ics::generate_uniform_distribution() {
    int cur_particle = 0;
    for (int i = 0; i < Np1; ++i) {
        for (int j = 0; j < Np1; ++j) {
            for (int k = 0; k < Np1; ++k) {
                plist->x->index(cur_particle, 0) = double(i)*Ng1/Np1;
                plist->x->index(cur_particle, 1) = double(j)*Ng1/Np1;
                plist->x->index(cur_particle, 2) = double(k)*Ng1/Np1;

                cur_particle++;
            }
        }
    }
}

void zeldovich_ics::fourier_displacement_vec() {
    double ak;
    double bk;
    double k2;
    double sigma8 = 0.835;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    for (int i=0; i<Np1; i++) {
        double kx = 2*M_PI*(i-Ng1/2+1)/Ng1;

        for (int j=0; j<Np1; j++) {
            double ky = 2*M_PI*(j-Ng1+1)/Ng1;

            for (int k=0; k<Np1/2; k++) {
                double kz = 2*M_PI*(k-Ng1+1)/Ng1;

                if (i==0 && j==0 && k==0) {
                    ak = 0;
                    bk = 0;
                } else{
                    k2 = kx*kx + ky*ky + kz*kz;
                    std::normal_distribution<> d{0, sqrt(cosmo.P_at_k(sqrt(k2)) / k2)};

                    ak = sigma8*d(gen)/2;
                    bk = sigma8*d(gen)/2;
                }

                fourier_Sx(i, j, k)[0] = kx*ak;
                fourier_Sx(i, j, k)[1] = -kx*bk;
                fourier_Sy(i, j, k)[0] = ky*ak;
                fourier_Sy(i, j, k)[1] = -ky*bk;
                fourier_Sz(i, j, k)[0] = kz*ak;
                fourier_Sz(i, j, k)[1] = -kz*bk;
            }
        }
    }
}

void zeldovich_ics::real_displacement_vec() {
    fftw_execute(backward_planx);
    fftw_execute(backward_plany);
    fftw_execute(backward_planz);
}

void zeldovich_ics::apply_ZA() {
    int cur_particle = 0;
    for (int i = 0; i < Np1; ++i) {
        for (int j = 0; j < Np1; ++j) {
            for (int k = 0; k < Np1; ++k) {
                plist->x->index(cur_particle, 0) -= cosmo.lgf(a0)*real_Sx(i, j, k)/pow(Np1, 1.5);
                plist->x->index(cur_particle, 1) -= cosmo.lgf(a0)*real_Sy(i, j, k)/pow(Np1, 1.5);
                plist->x->index(cur_particle, 2) -= cosmo.lgf(a0)*real_Sz(i, j, k)/pow(Np1, 1.5);
                printf("%lf\n", cosmo.lgf(a0)*real_Sx(i, j, k)/pow(Np1, 1.5));

                cur_particle++;
            }
        }
    }
}

void zeldovich_ics::setup_ics() {
    generate_uniform_distribution();

    fourier_displacement_vec();

    real_displacement_vec();

    apply_ZA();
}