//
// Created by pierfied on 11/19/18.
//

#ifndef STUPID_GRID_H
#define STUPID_GRID_H


#include <iostream>
#include <fftw3.h>
#include "multidim_array.h"
#include "particle_list.h"
#include "cosmology.h"

class grid {
public:
    const int nx;
    const int ny;
    const int nz;

    particle_list *plist;

    cosmology cosmo;

    array_3d<double> real_grid;
    array_3d<fftw_complex> fourier_grid;

    fftw_plan forward_plan;
    fftw_plan backward_plan;

    grid(int nx, int ny, int nz, particle_list &plist, cosmology cosmo) : nx(nx), ny(ny), nz(nz), real_grid(nx, ny, nz),
                                                                          fourier_grid(nx, ny, (nz / 2 + 1)),
                                                                          plist(&plist),
                                                                          forward_plan(fftw_plan_dft_r2c_3d(nx, ny, nz,
                                                                                                            real_grid.data,
                                                                                                            fourier_grid.data,
                                                                                                            FFTW_MEASURE)),
                                                                          backward_plan(
                                                                                  fftw_plan_dft_c2r_3d(nx, ny, nz,
                                                                                                       fourier_grid.data,
                                                                                                       real_grid.data,
                                                                                                       FFTW_MEASURE)),
                                                                          cosmo(cosmo) {}

    ~grid() {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
    }

    virtual void populate_delta_grid() = 0;

    void fft();

    void ifft();

    void apply_greens_func(double a);

    void compute_potential(double a);
};


#endif //STUPID_GRID_H
