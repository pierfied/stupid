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
    cic_grid grid(n, n, n, plist);

    grid.populate_delta_grid();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_DOUBLE_EQ(0, grid.delta_grid(i, j, k));
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
    cic_grid grid(n, n, n, plist);

    grid.populate_delta_grid();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                if (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) {
                    EXPECT_DOUBLE_EQ(7, grid.delta_grid(i, j, k));
                } else {
                    EXPECT_DOUBLE_EQ(-1, grid.delta_grid(i, j, k));
                }
            }
        }
    }
}