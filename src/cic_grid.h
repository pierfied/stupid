//
// Created by pierfied on 11/20/18.
//

#ifndef STUPID_CIC_GRID_H
#define STUPID_CIC_GRID_H


#include "grid.h"

class cic_grid : public grid {
public:
    cic_grid(long n, particle_list &plist, cosmology cosmo) : grid(n, plist, cosmo) {}

    void populate_delta_grid() override;

    inline double particle_accel(long particle_ind, long dim) override {
        long i = modulo(plist->x->index(particle_ind, 0), n);
        long j = modulo(plist->x->index(particle_ind, 1), n);
        long k = modulo(plist->x->index(particle_ind, 2), n);

        long i_neighbor = modulo(i + 1, n);
        long j_neighbor = modulo(j + 1, n);
        long k_neighbor = modulo(k + 1, n);

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
};


#endif //STUPID_CIC_GRID_H
