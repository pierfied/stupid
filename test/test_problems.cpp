//
// Created by rafa on 11/29/18.
//

#include <gtest/gtest.h>
#include <cosmology.h>
#include <tsc_grid.h>
#include <leapfrog_integrator.h>
#include <zeldovich_ics.h>


TEST(problems, one_dimensional_plane_wave) {
    cosmology cosmo(1, 0, 0, 0);

    const int n = 32;

    int npart = 16;
    int num_particles = npart*npart*npart;
    int num_dims = 3;

    double a0 = 0.2;
    double af = 0.6;
    int num_steps = 1000;
    double delta_a = (af - a0) / num_steps;

    double across = 5*a0;
    double wave_num = 2*M_PI/n;
    double amplitude = 1./wave_num/cosmo.D(across);

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < npart; ++i) {
        for (int j = 0; j < npart; ++j) {
            for (int k = 0; k < npart; ++k) {
                x(cur_particle, 0) = double(i)*n/npart + cosmo.D(a0) * amplitude * sin(wave_num*i*n/npart);
                x(cur_particle, 1) = double(j)*n/npart;
                x(cur_particle, 2) = double(k)*n/npart;

                p(cur_particle, 0) = a0 * a0 * cosmo.Ddot(a0) * amplitude * sin(wave_num*i*n/npart);

                cur_particle++;
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);
    grid.rebound_positions();

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    cur_particle = 0;
    for (int i = 0; i < npart; ++i) {
        for (int j = 0; j < npart; ++j) {
            for (int k = 0; k < npart; ++k) {
                double analytic = grid.modulo(double(i)*n/npart + cosmo.D(af) * amplitude * sin(wave_num*i*n/npart), n);
                EXPECT_NEAR(x(cur_particle, 0), analytic, 1e-2*analytic);

                cur_particle++;
            }
        }
    }
}