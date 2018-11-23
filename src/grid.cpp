//
// Created by pierfied on 11/19/18.
//

#include "grid.h"

#include <fftw3.h>

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