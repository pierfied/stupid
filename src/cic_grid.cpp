//
// Created by pierfied on 11/20/18.
//

#include "cic_grid.h"

void cic_grid::populate_delta_grid() {
    real_grid.reset_zero();

    for (int p = 0; p < plist->num_particles; ++p) {
        int i = int(plist->x->index(p, 0)) % nx;
        int j = int(plist->x->index(p, 1)) % ny;
        int k = int(plist->x->index(p, 2)) % nz;

        int i_neighbor = (i + 1) % nx;
        int j_neighbor = (j + 1) % ny;
        int k_neighbor = (k + 1) % nz;

        double dx = plist->x->index(p, 0) - i;
        double dy = plist->x->index(p, 1) - j;
        double dz = plist->x->index(p, 2) - k;

        double tx = 1 - dx;
        double ty = 1 - dy;
        double tz = 1 - dz;

        real_grid(i, j, k) += tx * ty * tz;

        real_grid(i_neighbor, j, k) += dx * ty * tz;
        real_grid(i, j_neighbor, k) += tx * dy * tz;
        real_grid(i, j, k_neighbor) += tx * ty * dz;

        real_grid(i_neighbor, j_neighbor, k) += dx * dy * tz;
        real_grid(i_neighbor, j, k_neighbor) += dx * ty * dz;
        real_grid(i, j_neighbor, k_neighbor) += tx * dy * dz;

        real_grid(i_neighbor, j_neighbor, k_neighbor) += dx * dy * dz;
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                real_grid(i, j, k) = real_grid(i, j, k) / avg_particle_density - 1;
            }
        }
    }
}

inline double cic_grid::particle_accel(int particle_ind, int dim) {
    int i = int(plist->x->index(particle_ind, 0)) % nx;
    int j = int(plist->x->index(particle_ind, 1)) % ny;
    int k = int(plist->x->index(particle_ind, 2)) % nz;

    int i_neighbor = (i + 1) % nx;
    int j_neighbor = (j + 1) % ny;
    int k_neighbor = (k + 1) % nz;

    double dx = plist->x->index(particle_ind, 0) - i;
    double dy = plist->x->index(particle_ind, 1) - j;
    double dz = plist->x->index(particle_ind, 2) - k;

    double tx = 1 - dx;
    double ty = 1 - dy;
    double tz = 1 - dz;

    double g = grid_accel(i, j, k, dim) * tx * ty * tz +
               grid_accel(i_neighbor, j, k, dim) * dx * ty * tz +
               grid_accel(i, j_neighbor, k, dim) * tx * dy * tz +
               grid_accel(i, j, k_neighbor, dim) * tx * ty * dz +
               grid_accel(i_neighbor, j_neighbor, k, dim) * dx * dy * tz +
               grid_accel(i_neighbor, j, k_neighbor, dim) * dx * ty * dz +
               grid_accel(i, j_neighbor, k_neighbor, dim) * tx * dy * dz +
               grid_accel(i_neighbor, j_neighbor, k_neighbor, dim) * dx * dy * dz;

    return g;
}
