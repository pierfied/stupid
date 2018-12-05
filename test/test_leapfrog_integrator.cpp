//
// Created by pierfied on 11/25/18.
//

#include <gtest/gtest.h>
#include <cic_grid.h>
#include <leapfrog_integrator.h>
#include <filesystem>
#include <fstream>

TEST(leapfrog_integrator, motion_2body_x) {
    int n = 10;

    int num_particles = 2;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    x(1, 0) = 1;

    cosmology cosmo(1, 0, 0, 0);

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    EXPECT_NEAR(0, x(0, 1), 1e-8);
    EXPECT_NEAR(0, x(0, 2), 1e-8);
    EXPECT_NEAR(0, x(1, 1), 1e-8);
    EXPECT_NEAR(0, x(1, 2), 1e-8);

    EXPECT_LT(fabs(x(1, 0) - x(0, 0)), 1);
}

TEST(leapfrog_integrator, motion_2body_y) {
    int n = 10;

    int num_particles = 2;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    x(1, 1) = 1;

    cosmology cosmo(1, 0, 0, 0);

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    EXPECT_NEAR(0, x(0, 0), 1e-8);
    EXPECT_NEAR(0, x(0, 2), 1e-8);
    EXPECT_NEAR(0, x(1, 0), 1e-8);
    EXPECT_NEAR(0, x(1, 2), 1e-8);

    EXPECT_LT(fabs(x(1, 1) - x(0, 1)), 1);
}

TEST(leapfrog_integrator, motion_2body_z) {
    int n = 10;

    int num_particles = 2;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    x(1, 2) = 1;

    cosmology cosmo(1, 0, 0, 0);

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    EXPECT_NEAR(0, x(0, 1), 1e-8);
    EXPECT_NEAR(0, x(0, 0), 1e-8);
    EXPECT_NEAR(0, x(1, 1), 1e-8);
    EXPECT_NEAR(0, x(1, 0), 1e-8);

    EXPECT_LT(fabs(x(1, 2) - x(0, 2)), 1);
}

TEST(leapfrog_integrator, flat_density) {
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
    cic_grid grid(n, plist, cosmology(1, 0, 0, 0));

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_NEAR(i, x(cur_particle, 0), 1e-8);
                EXPECT_NEAR(j, x(cur_particle, 1), 1e-8);
                EXPECT_NEAR(k, x(cur_particle, 2), 1e-8);

                cur_particle++;
            }
        }
    }
}

TEST(leapfrog_integrator, alternating_density) {
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
    cic_grid grid(n, plist, cosmology(1, 0, 0, 0));

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    cur_particle = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_NEAR(i - i % 2, x(cur_particle, 0), 1e-8);
                EXPECT_NEAR(j - j % 2, x(cur_particle, 1), 1e-8);
                EXPECT_NEAR(k - k % 2, x(cur_particle, 2), 1e-8);

                cur_particle++;
            }
        }
    }
}

TEST(leapfrog_integrator, particle_cube) {
    int n = 10;

    int num_particles = 8;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                x(cur_particle, 0) = 2 * i;
                x(cur_particle, 1) = 2 * j;
                x(cur_particle, 2) = 2 * k;

                cur_particle++;
            }
        }
    }

    cosmology cosmo(1, 0, 0, 0);

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    cur_particle = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                EXPECT_LT(fabs(x(cur_particle, 0) - 1), 1);
                EXPECT_LT(fabs(x(cur_particle, 1) - 1), 1);
                EXPECT_LT(fabs(x(cur_particle, 2) - 1), 1);

                cur_particle++;
            }
        }
    }
}

TEST(leapfrog_integrator, pos_write) {
    int n = 10;

    int num_particles = 8;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                x(cur_particle, 0) = 2 * i;
                x(cur_particle, 1) = 2 * j;
                x(cur_particle, 2) = 2 * k;

                cur_particle++;
            }
        }
    }

    cosmology cosmo(1, 0, 0, 0);

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    std::filesystem::create_directory("test_output");

    std::string fprefix = "test_output/test";

    leapfrog_integrator li(a0, af, delta_a, grid, fprefix, num_steps);

    li.run_sim();

    std::string fname = "test_output/test_x_1.bin";

    std::ifstream file(fname, std::ios::in | std::ios::binary);

    double af_written;
    file.read((char *) &af_written, sizeof(double));

    EXPECT_DOUBLE_EQ(af, af_written);

    array_2d<double> x_written(num_particles, num_dims);
    file.read((char *) x_written.data, sizeof(double) * num_particles * num_dims);

    for (int i = 0; i < num_particles; ++i) {
        EXPECT_DOUBLE_EQ(x(i, 0), x_written(i, 0));
        EXPECT_DOUBLE_EQ(x(i, 1), x_written(i, 1));
        EXPECT_DOUBLE_EQ(x(i, 2), x_written(i, 2));
    }

    file.close();

    std::filesystem::remove_all("test_output");
}

TEST(leapfrog_integrator, mom_write) {
    int n = 10;

    int num_particles = 8;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                x(cur_particle, 0) = 2 * i;
                x(cur_particle, 1) = 2 * j;
                x(cur_particle, 2) = 2 * k;

                cur_particle++;
            }
        }
    }

    cosmology cosmo(1, 0, 0, 0);

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);

    double a0 = 0.5;
    double af = 1;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    std::filesystem::create_directory("test_output");

    std::string fprefix = "test_output/test";

    leapfrog_integrator li(a0, af, delta_a, grid, fprefix, num_steps);

    li.run_sim();

    std::string fname = "test_output/test_p_1.bin";

    std::ifstream file(fname, std::ios::in | std::ios::binary);

    double af_written;
    file.read((char *) &af_written, sizeof(double));

    EXPECT_DOUBLE_EQ(af, af_written);

    array_2d<double> p_written(num_particles, num_dims);
    file.read((char *) p_written.data, sizeof(double) * num_particles * num_dims);

    for (int i = 0; i < num_particles; ++i) {
        EXPECT_DOUBLE_EQ(p(i, 0), p_written(i, 0));
        EXPECT_DOUBLE_EQ(p(i, 1), p_written(i, 1));
        EXPECT_DOUBLE_EQ(p(i, 2), p_written(i, 2));
    }

    file.close();

    std::filesystem::remove_all("test_output");
}