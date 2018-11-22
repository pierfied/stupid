//
// Created by pierfied on 11/21/18.
//

#include <gtest/gtest.h>
#include <multidim_array.h>
#include <cic_grid.h>

TEST(cic_grid, one_mass_per_cell) {
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
            }
        }
    }

    particle_list plist(x, p);
    cic_grid grid(n, n, n, plist);

    grid.populate_mass_grid();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_DOUBLE_EQ(1, grid.mass_grid(i, j, k));
            }
        }
    }
}