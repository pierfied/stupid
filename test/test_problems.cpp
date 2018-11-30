//
// Created by rafa on 11/29/18.
//

#include <gtest/gtest.h>
#include <cic_grid.h>
#include <leapfrog_integrator.h>
#include <filesystem>
#include <fstream>
#include <math.h>
#include <cosmology.h>


TEST(problems, one_dimensional_plane_wave) {
    cosmology cosmo(1, 0, 0);

    int n = 32;

    int npart = 16;
    int num_particles = npart*npart*npart;
    int num_dims = 3;

    double a0 = 0.5;
    double af = 0.6;
    int num_steps = 10000;
    double delta_a = (af - a0) / num_steps;

    double across = 20*a0;
    double wave_num = 2*M_PI/n;
    double amplitude = 1./wave_num/cosmo.lgf(across);

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < npart; ++i) {
        for (int j = 0; j < npart; ++j) {
            for (int k = 0; k < npart; ++k) {
                x(cur_particle, 0) = double(i)*n/npart + cosmo.lgf(a0)*amplitude*sin(wave_num*i*n/npart);
                x(cur_particle, 1) = double(j)*n/npart;
                x(cur_particle, 2) = double(k)*n/npart;

                p(cur_particle, 0) = a0*cosmo.lgf(a0-0.5*delta_a*amplitude*sin(wave_num*i*n/npart));

                cur_particle++;
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    cur_particle = 0;
    for (int i = 0; i < npart; ++i) {
        for (int j = 0; j < npart; ++j) {
            for (int k = 0; k < npart; ++k) {
                EXPECT_FLOAT_EQ(x(cur_particle, 0), double(i)*n/npart + cosmo.lgf(af)*amplitude*sin(wave_num*i*n/npart));

                cur_particle++;
            }
        }
    }
}
