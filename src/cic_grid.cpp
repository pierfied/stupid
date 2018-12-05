//
// Created by pierfied on 11/20/18.
//

#include "cic_grid.h"

void cic_grid::populate_delta_grid() {
    real_grid.reset_zero();

    for (int p = 0; p < plist->num_particles; ++p) {
        int i = modulo(plist->x->index(p, 0), n);
        int j = modulo(plist->x->index(p, 1), n);
        int k = modulo(plist->x->index(p, 2), n);

        int i_neighbor = (i + 1) % n;
        int j_neighbor = (j + 1) % n;
        int k_neighbor = (k + 1) % n;

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

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                real_grid(i, j, k) = real_grid(i, j, k) / avg_particle_density - 1;
            }
        }
    }
}