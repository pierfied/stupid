//
// Created by pierfied on 11/20/18.
//

#ifndef STUPID_TSC_GRID_H
#define STUPID_TSC_GRID_H


#include "grid.h"

class tsc_grid : public grid {
public:
    tsc_grid(int nx, int ny, int nz, particle_list &plist, cosmology cosmo) : grid(nx, ny, nz, plist, cosmo) {}

    void populate_delta_grid() override;

    inline double particle_accel(int particle_ind, int dim) override {
        int i = int(plist->x->index(particle_ind, 0)) % nx;
        int j = int(plist->x->index(particle_ind, 1)) % ny;
        int k = int(plist->x->index(particle_ind, 2)) % nz;

        int iu = (i + 1) % nx;
        int ju = (j + 1) % ny;
        int ku = (k + 1) % nz;
        int il = ((i - 1) % nx + nx) % nx;
        int jl = ((j - 1) % ny + ny) % ny;
        int kl = ((k - 1) % nz + nz) % nz;

        double dx = plist->x->index(particle_ind, 0) - i;
        double dy = plist->x->index(particle_ind, 1) - j;
        double dz = plist->x->index(particle_ind, 2) - k;

        double tx = 0.75 - dx*dx;
        double ty = 0.75 - dy*dy;
        double tz = 0.75 - dz*dz;

        double srx = 0.5*(0.5 + abs(dx))*(0.5 + abs(dx));
        double sry = 0.5*(0.5 + abs(dy))*(0.5 + abs(dy));
        double srz = 0.5*(0.5 + abs(dz))*(0.5 + abs(dz));

        double lrx = 0.5*(0.5 - abs(dx))*(0.5 - abs(dx));
        double lry = 0.5*(0.5 - abs(dy))*(0.5 - abs(dy));
        double lrz = 0.5*(0.5 - abs(dz))*(0.5 - abs(dz));

        if (abs(dx) > 1.5) lrx = 0.;
        if (abs(dy) > 1.5) lry = 0.;
        if (abs(dz) > 1.5) lrz = 0.;

        double ux, uy, uz, lx, ly, lz;

        if (dx >= 0.) {
            ux = srx;
            lx = lrx;
        }
        else if (dx < 0.) {
            ux = lrx;
            lx = srx;
        }

        if (dy >= 0.) {
            uy = sry;
            ly = lry;
        }
        else if (dy < 0.) {
            uy = lry;
            ly = sry;
        }

        if (dz >= 0.) {
            uz = srz;
            lz = lrz;
        }
        else if (dz < 0.) {
            uz = lrz;
            lz = srz;
        }

        double g = grid_accel(i, j, k, dim) * tx * ty * tz + 
                   grid_accel(iu, j, k, dim) * ux * ty * tz +
                   grid_accel(i, ju, k, dim) * tx * uy * tz +
                   grid_accel(i, j, ku, dim) * tx * ty * uz +
                   grid_accel(il, j, k, dim) * lx * ty * tz +
                   grid_accel(i, jl, k, dim) * tx * ly * tz +
                   grid_accel(i, j, kl, dim) * tx * ty * lz +

                   grid_accel(iu, ju, k, dim) * ux * uy * tz +
                   grid_accel(iu, j, ku, dim) * ux * ty * uz +
                   grid_accel(i, ju, ku, dim) * tx * uy * uz +
                   grid_accel(il, jl, k, dim) * lx * ly * tz +
                   grid_accel(il, j, kl, dim) * lx * ty * lz +
                   grid_accel(i, jl, kl, dim) * tx * ly * lz +
                   grid_accel(iu, jl, k, dim) * ux * ly * tz +
                   grid_accel(iu, j, kl, dim) * ux * ty * lz +
                   grid_accel(i, ju, kl, dim) * tx * uy * lz +
                   grid_accel(il, ju, k, dim) * lx * uy * tz +
                   grid_accel(il, j, ku, dim) * lx * ty * uz +
                   grid_accel(i, jl, ku, dim) * tx * ly * uz +

                   grid_accel(iu, ju, ku, dim) * ux * uy * uz +
                   grid_accel(iu, ju, kl, dim) * ux * uy * lz +
                   grid_accel(iu, jl, ku, dim) * ux * ly * uz +
                   grid_accel(iu, jl, kl, dim) * ux * ly * lz +
                   grid_accel(il, ju, ku, dim) * lx * uy * uz +
                   grid_accel(il, ju, kl, dim) * lx * uy * lz +
                   grid_accel(il, jl, ku, dim) * lx * ly * uz +
                   grid_accel(il, jl, kl, dim) * lx * ly * lz;

        return g;
    }
};


#endif //STUPID_TSC_GRID_H
