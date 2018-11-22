//
// Created by pierfied on 11/20/18.
//

#include "cic_grid.h"

void cic_grid::populate_mass_grid() {
    mass_grid.reset_zero();

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

        mass_grid(i, j, k) += tx * ty * tz;

        mass_grid(i_neighbor, j, k) += dx * ty * tz;
        mass_grid(i, j_neighbor, k) += tx * dy * tz;
        mass_grid(i, j, k_neighbor) += tx * ty * dz;

        mass_grid(i_neighbor, j_neighbor, k) += dx * dy * tz;
        mass_grid(i_neighbor, j, k_neighbor) += dx * ty * dz;
        mass_grid(i, j_neighbor, k_neighbor) += tx * dy * dz;

        mass_grid(i_neighbor, j_neighbor, k_neighbor) += dx * dy * dz;
    }
}