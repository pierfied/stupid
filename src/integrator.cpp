//
// Created by pierfied on 11/24/18.
//

#include <fstream>
#include "integrator.h"

void integrator::write_positions(double a) {
    std::string fname = file_prefix + "_x_" + std::to_string(write_num) + ".bin";
    std::ofstream file(fname, std::ios::out | std::ios::binary);

    if (g->cosmo.use_real_units) {
#pragma omp parallel for
        for (long i = 0; i < g->plist->num_particles; ++i) {
            g->plist->x->index(i,0) *= a * g->cosmo.u.r0;
            g->plist->x->index(i,1) *= a * g->cosmo.u.r0;
            g->plist->x->index(i,2) *= a * g->cosmo.u.r0;
        }
    }

    file.write((char *) &a, sizeof(double));
    file.write((char *) g->plist->x->data, sizeof(double) * g->plist->num_particles * g->plist->num_dims);

    if (g->cosmo.use_real_units) {
#pragma omp parallel for
        for (long i = 0; i < g->plist->num_particles; ++i) {
            g->plist->x->index(i,0) /= a * g->cosmo.u.r0;
            g->plist->x->index(i,1) /= a * g->cosmo.u.r0;
            g->plist->x->index(i,2) /= a * g->cosmo.u.r0;
        }
    }

    file.close();
}

void integrator::write_momentum(double a) {
    std::string fname = file_prefix + "_p_" + std::to_string(write_num) + ".bin";
    std::ofstream file(fname, std::ios::out | std::ios::binary);

    double v0 = g->cosmo.u.r0 / g->cosmo.u.t0;

    if (g->cosmo.use_real_units) {
#pragma omp parallel for
        for (long i = 0; i < g->plist->num_particles; ++i) {
            g->plist->p->index(i,0) *= v0 / a;
            g->plist->p->index(i,1) *= v0 / a;
            g->plist->p->index(i,2) *= v0 / a;
        }
    }

    file.write((char *) &a, sizeof(double));
    file.write((char *) g->plist->p->data, sizeof(double) * g->plist->num_particles * g->plist->num_dims);

    if (g->cosmo.use_real_units) {
#pragma omp parallel for
        for (long i = 0; i < g->plist->num_particles; ++i) {
            g->plist->p->index(i,0) /= v0 / a;
            g->plist->p->index(i,1) /= v0 / a;
            g->plist->p->index(i,2) /= v0 / a;
        }
    }

    file.close();
}

void integrator::write_pos_and_mom(double a) {
    write_positions(a);
    write_momentum(a);

    write_num++;
}
