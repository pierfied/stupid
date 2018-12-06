//
// Created by rafa on 11/29/18.
//

#include <gtest/gtest.h>
#include <cosmology.h>
#include <tsc_grid.h>
#include <leapfrog_integrator.h>
#include <zeldovich_ics.h>


/*TEST(problems, one_dimensional_plane_wave) {
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
*/
TEST(problems, satellite) {
    cosmology cosmo(1, 0, 0, 0);

    const int n = 100;

    int num_particles = 10010;
    int num_dims = 3;

    double a0 = 1.;
    double af = 1.001;
    int num_steps = 10;
    double delta_a = (af - a0) / num_steps;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);



    int cur_particle = 0;
    for (int i = 0; i < 10000; ++i) {
        x(cur_particle, 0) = n/2.;
        x(cur_particle, 1) = n/2.;
        x(cur_particle, 2) = n/2.;

        p(cur_particle, 0) = 0.;
        p(cur_particle, 1) = 0.;
        p(cur_particle, 2) = 0.;

        cur_particle++;
    }

    for (int i = 0; i < 10; ++i) {
        x(cur_particle, 0) = n/2.;
        x(cur_particle, 1) = n/2. + 2;
        x(cur_particle, 2) = n/2.;

        p(cur_particle, 0) = 0.;
        p(cur_particle, 1) = pow(a0, 7.5) * sqrt(10000./2./8./M_PI);
        p(cur_particle, 2) = 0.;

        cur_particle++;
    }

    particle_list plist(x, p);


    zeldovich_ics ics(cosmo, plist, )



    cic_grid grid(n, plist, cosmo);
    //grid.rebound_positions();

    cur_particle = 0;
    for (int i = 0; i < 100; ++i) {
        double analytic = n/2;
        EXPECT_NEAR(x(cur_particle, 2), analytic, 1e-8);

        cur_particle++;
    }

    leapfrog_integrator li(a0, af, delta_a, grid);

    li.run_sim();

    cur_particle = 0;
    for (int i = 0; i < 100; ++i) {
        double analytic = n/2;
        EXPECT_NEAR(x(cur_particle, 2), analytic, 1e-8);

        cur_particle++;
    }
}

/*TEST(problems, pos_write) {
    cosmology cosmo(1, 0, 0, 0);

    const int n = 100;

    int num_particles = 10010;
    int num_dims = 3;

    double a0 = 1.;
    double af = 1.007;
    int num_steps = 100;
    double delta_a = (af - a0) / num_steps;

    array_2d<double> x(num_particles, num_dims);
    array_2d<double> p(num_particles, num_dims);

    int cur_particle = 0;
    for (int i = 0; i < 10000; ++i) {
        x(cur_particle, 0) = n/2.;
        x(cur_particle, 1) = n/2.;
        x(cur_particle, 2) = n/2.;

        p(cur_particle, 0) = 0.;
        p(cur_particle, 1) = 0.;
        p(cur_particle, 2) = 0.;

        cur_particle++;
    }

    for (int i = 0; i < 10; ++i) {
        x(cur_particle, 0) = n/2. + 2;
        x(cur_particle, 1) = n/2.;
        x(cur_particle, 2) = n/2.;

        p(cur_particle, 0) = 0.;
        p(cur_particle, 1) = pow(a0, 7.5) * sqrt(10000000./2./8./M_PI);
        p(cur_particle, 2) = 0.;

        cur_particle++;
    }

    particle_list plist(x, p);
    cic_grid grid(n, plist, cosmo);
    //grid.rebound_positions();

    cur_particle = 0;
    for (int i = 0; i < 100; ++i) {
        double analytic = n/2;
        EXPECT_NEAR(x(cur_particle, 2), analytic, 1e-8);

        cur_particle++;
    }

    std::filesystem::create_directory("test_satellite");

    std::string fprefix = "test_satellite/test";

    leapfrog_integrator li(a0, af, delta_a, grid, fprefix, num_steps);

    li.run_sim();

    std::string fname = "test_satellite/test_x_1.bin";

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

    //std::filesystem::remove_all("test_output");
}*/