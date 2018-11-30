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


TEST(leapfrog_integrator, one_dimensional_plane_wave) {
    cosmology cosmo(1, 0, 0, 0);

    int n = 32;

    int num_particles = 64*64*64;
    int num_dims = 3;

    double a0 = 0.5;
    double af = 5;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    double across = 10*a0;
    double wave_num = 2*M_PI/n;
    double amplitude = 1./wave_num/cosmo.lgf(across);

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                x(cur_particle, 0) = double(i)*n/num_particles + cosmo.lgf(a0)*amplitude*sin(wave_num*i*n/num_particles);
                x(cur_particle, 1) = double(j)*n/num_particles + cosmo.lgf(a0)*amplitude*sin(wave_num*j*n/num_particles);
                x(cur_particle, 2) = double(k)*n/num_particles + cosmo.lgf(a0)*amplitude*sin(wave_num*k*n/num_particles);

                p(cur_particle, 0) = a0*cosmo.lgf(a0-0.5*delta_a*amplitude*sin(wave_num*i*n/num_particles));
                p(cur_particle, 1) = a0*cosmo.lgf(a0-0.5*delta_a*amplitude*sin(wave_num*j*n/num_particles));
                p(cur_particle, 2) = a0*cosmo.lgf(a0-0.5*delta_a*amplitude*sin(wave_num*k*n/num_particles));

                cur_particle++;
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_NEAR(i, x(cur_particle, 0), double(i)*n/num_particles + cosmo.lgf(af)*amplitude*sin(wave_num*i*n/num_particles));
                EXPECT_NEAR(j, x(cur_particle, 1), double(j)*n/num_particles + cosmo.lgf(af)*amplitude*sin(wave_num*j*n/num_particles));
                EXPECT_NEAR(k, x(cur_particle, 2), double(k)*n/num_particles + cosmo.lgf(af)*amplitude*sin(wave_num*k*n/num_particles));

                cur_particle++;
            }
        }
    }
}
