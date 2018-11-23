//
// Created by pierfied on 11/19/18.
//

#include "grid.h"

#include <fftw3.h>
#include <math.h>

void grid::fft() {
    fftw_execute(forward_plan);
}

void grid::ifft() {
    fftw_execute(backward_plan);

    double norm = nx * ny * nz;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                real_grid(i, j, k) /= norm;
            }
        }
    }
}

void grid::apply_greens_func(double a) {
    double factor = -3 * cosmo.Omega_0 / (8 * a);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < (nz / 2 + 1); ++k) {
                double sin_x = sin(M_PI * i / nx);
                double sin_y = sin(M_PI * j / ny);
                double sin_z = sin(M_PI * k / nz);

                double G = factor / (sin_x * sin_x + sin_y * sin_y + sin_z * sin_z);

                fourier_grid(i, j, k)[0] *= G;
                fourier_grid(i, j, k)[1] *= G;
            }
        }
    }

    fourier_grid(0, 0, 0)[0] = 0;
    fourier_grid(0, 0, 0)[1] = 0;
}

void grid::compute_potential(double a) {
    populate_delta_grid();

    fft();

    apply_greens_func(a);

    ifft();
}
