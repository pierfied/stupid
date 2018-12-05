//
// Created by pierfied on 11/20/18.
//

#include "tsc_grid.h"

void tsc_grid::populate_delta_grid() {
    real_grid.reset_zero();

    for (int p = 0; p < plist->num_particles; ++p) {
        int i = modulo(plist->x->index(p, 0), n);
        int j = modulo(plist->x->index(p, 1), n);
        int k = modulo(plist->x->index(p, 2), n);

        int iu = modulo((i + 1), n);
        int ju = modulo((j + 1), n);
        int ku = modulo((k + 1), n);
        int il = modulo((i - 1), n);
        int jl = modulo((j - 1), n);
        int kl = modulo((k - 1), n);

        double dx = plist->x->index(p, 0) - i;
        double dy = plist->x->index(p, 1) - j;
        double dz = plist->x->index(p, 2) - k;

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

        real_grid(i, j, k) += tx * ty * tz;

        real_grid(iu, j, k) += ux * ty * tz;
        real_grid(i, ju, k) += tx * uy * tz;
        real_grid(i, j, ku) += tx * ty * uz;
        real_grid(il, j, k) += lx * ty * tz;
        real_grid(i, jl, k) += tx * ly * tz;
        real_grid(i, j, kl) += tx * ty * lz;

        real_grid(iu, ju, k) += ux * uy * tz;
        real_grid(iu, j, ku) += ux * ty * uz;
        real_grid(i, ju, ku) += tx * uy * uz;
        real_grid(il, jl, k) += lx * ly * tz;
        real_grid(il, j, kl) += lx * ty * lz;
        real_grid(i, jl, kl) += tx * ly * lz;
        real_grid(iu, jl, k) += ux * ly * tz;
        real_grid(iu, j, kl) += ux * ty * lz;
        real_grid(i, ju, kl) += tx * uy * lz;
        real_grid(il, ju, k) += lx * uy * tz;
        real_grid(il, j, ku) += lx * ty * uz;
        real_grid(i, jl, ku) += tx * ly * uz;

        real_grid(iu, ju, ku) += ux * uy * uz;
        real_grid(iu, ju, kl) += ux * uy * lz;
        real_grid(iu, jl, ku) += ux * ly * uz;
        real_grid(iu, jl, kl) += ux * ly * lz;
        real_grid(il, ju, ku) += lx * uy * uz;
        real_grid(il, ju, kl) += lx * uy * lz;
        real_grid(il, jl, ku) += lx * ly * uz;
        real_grid(il, jl, kl) += lx * ly * lz;
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