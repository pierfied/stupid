//
// Created by pierfied on 11/20/18.
//

#ifndef STUPID_CIC_GRID_H
#define STUPID_CIC_GRID_H


#include "grid.h"

class cic_grid : public grid {
public:
    cic_grid(int nx, int ny, int nz, particle_list &plist, cosmology cosmo) : grid(nx, ny, nz, plist, cosmo) {}

    void populate_delta_grid() override;

    inline double particle_accel(int particle_ind, int dim) override {
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
};


#endif //STUPID_CIC_GRID_H
