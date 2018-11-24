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

double grid::grid_accel(int i, int j, int k, int dim) {
    int plus_neighbor;
    int minus_neighbor;

    switch (dim) {
        case 0:
            plus_neighbor = (i + 1) % nx;
            minus_neighbor = (i - 1) % nx;

            return -(real_grid(plus_neighbor, j, k) - real_grid(minus_neighbor, j, k)) / 2;
        case 1:
            plus_neighbor = (j + 1) % ny;
            minus_neighbor = (j - 1) % ny;

            return -(real_grid(i, plus_neighbor, k) - real_grid(i, minus_neighbor, k)) / 2;
        default:
            plus_neighbor = (k + 1) % nz;
            minus_neighbor = (k - 1) % nz;

            return -(real_grid(i, i, plus_neighbor) - real_grid(i, j, minus_neighbor)) / 2;
    }
}
