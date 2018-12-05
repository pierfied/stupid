#include <iostream>

#include "stupid.h"
#include "multidim_array.h"
#include "integrator.h"
#include "leapfrog_integrator.h"
#include "cic_grid.h"
#include "tsc_grid.h"
#include "zeldovich_ics.h"

#define NUM_DIMS 3

void hello() {
    std::cout << "Hello, World!" << std::endl;
}

void run_sim(stupid_args args) {
    array_2d<double> x(args.num_particles, NUM_DIMS);
    array_2d<double> p(args.num_particles, NUM_DIMS);

    for (int i = 0; i < args.num_particles * NUM_DIMS; ++i) {
        x.data[i] = args.x[i];
        p.data[i] = args.p[i];
    }

    particle_list plist(x, p);

    cosmology *cosmo;

    if (args.use_real_units) {
        cosmo = new cosmology(args.Omega_m0, args.Omega_k0, args.Omega_l0, args.sigma8);
    }else{
        cosmo = new cosmology(args.Omega_m0, args.Omega_k0, args.Omega_l0, args.sigma8);
    }

    grid *g;

    switch (args.interp_scheme) {
        case 1:
            g = new cic_grid(args.num_cells, plist, *cosmo);
        default:
            g = new tsc_grid(args.num_cells, plist, *cosmo);
    }

    integrator *i;

    std::string file_prefix(args.file_prefix);
    switch (args.integrator) {
        default:
            i = new leapfrog_integrator(args.a0, args.af, args.delta_a, *g, file_prefix, args.write_nth_step);
    }

    i->run_sim();

    delete g;
    delete i;
    delete cosmo;
}

void setup_ics(stupid_ics_args args) {
    printf("a\n");
    array_2d<double> x(args.num_particles, NUM_DIMS);
    array_2d<double> p(args.num_particles, NUM_DIMS);

    particle_list *plist;

    plist = new particle_list(x, p);

    cosmology *cosmo;

    cosmo = new cosmology(args.Omega_m0, args.Omega_k0, args.Omega_l0, args.sigma8);

    zeldovich_ics *ics;

    ics = new zeldovich_ics(*cosmo, *plist, args.num_cells, args.k, args.P, args.sizeofk, args.a0);

    ics->setup_ics();

    for (int i = 0; i < args.num_particles * NUM_DIMS; ++i) {
        args.x[i] = x.data[i];
        args.p[i] = p.data[i];
    }

    delete ics;
    delete cosmo;
}