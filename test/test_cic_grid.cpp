//
// Created by pierfied on 11/21/18.
//

#include <gtest/gtest.h>
#include <multidim_array.h>
#include <cic_grid.h>

TEST(cic_grid, flat_density) {
    int n = 10;

    int num_particles = n * n * n;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                x(cur_particle, 0) = i;
                x(cur_particle, 1) = j;
                x(cur_particle, 2) = k;

                cur_particle++;
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, n, n, plist, cosmology(0, 0, 0, 0));

    grid.populate_delta_grid();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_DOUBLE_EQ(0, grid.real_grid(i, j, k));
            }
        }
    }
}

TEST(cic_grid, alternating_density) {
    int n = 10;

    int num_particles = n * n * n;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                x(cur_particle, 0) = i - i % 2;
                x(cur_particle, 1) = j - j % 2;
                x(cur_particle, 2) = k - k % 2;

                cur_particle++;
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, n, n, plist, cosmology(0, 0, 0, 0));

    grid.populate_delta_grid();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                if (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) {
                    EXPECT_DOUBLE_EQ(7, grid.real_grid(i, j, k));
                } else {
                    EXPECT_DOUBLE_EQ(-1, grid.real_grid(i, j, k));
                }
            }
        }
    }
}

TEST(cic_grid, fft_recover) {
    int n = 10;

    int num_particles = n * n * n;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                x(cur_particle, 0) = i - i % 2;
                x(cur_particle, 1) = j - j % 2;
                x(cur_particle, 2) = k - k % 2;

                cur_particle++;
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, n, n, plist, cosmology(0, 0, 0, 0));

    grid.populate_delta_grid();

    grid.fft();
    grid.real_grid.reset_zero();
    grid.ifft();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                if (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) {
                    EXPECT_DOUBLE_EQ(7, grid.real_grid(i, j, k));
                } else {
                    EXPECT_DOUBLE_EQ(-1, grid.real_grid(i, j, k));
                }
            }
        }
    }
}

TEST(cic_grid, flat_potential) {
    int n = 10;

    int num_particles = n * n * n;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                x(cur_particle, 0) = i;
                x(cur_particle, 1) = j;
                x(cur_particle, 2) = k;

                cur_particle++;
            }
        }
    }

    cosmology cosmo(0.3, 0, 0.7, 70);

    particle_list plist(x, p);
    cic_grid grid(n, n, n, plist, cosmo);

    double a = 1;
    grid.compute_potential(a);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_EQ(0, grid.real_grid(i, j, k));
            }
        }
    }
}