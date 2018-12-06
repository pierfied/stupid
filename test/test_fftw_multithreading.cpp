//
// Created by rafa on 11/27/18.
//

#include <gtest/gtest.h>
#include <cic_grid.h>
#include <leapfrog_integrator.h>
#include <filesystem>
#include <fstream>
#include <fftw3.h>
#include <omp.h>

TEST(fftw_multithreading, flat_density) {

    int n = 64;

    int num_particles = 128*128*128;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    srand (time(NULL));

    int cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                x(cur_particle, 0) = rand()/RAND_MAX*n;
                x(cur_particle, 1) = rand()/RAND_MAX*n;
                x(cur_particle, 2) = rand()/RAND_MAX*n;

                cur_particle++;
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmology(1, 0, 0, 0));

    grid.fft();

}