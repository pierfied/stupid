//
// Created by pierfied on 11/20/18.
//

#include "tsc_grid.h"

void tsc_grid::populate_delta_grid() {
    real_grid.reset_zero();

    for (int p = 0; p < plist->num_particles; ++p) {
        int i = int(plist->x->index(p, 0)) % nx;
        int j = int(plist->x->index(p, 1)) % ny;
        int k = int(plist->x->index(p, 2)) % nz;

        /*
        int i_neighbor = (i + 1) % nx;
        int j_neighbor = (j + 1) % ny;
        int k_neighbor = (k + 1) % nz;
        */

        double dx = plist->x->index(p, 0) - i;
        double dy = plist->x->index(p, 1) - j;
        double dz = plist->x->index(p, 2) - k;

        double srx = 0.75 - dx*dx;
        double sry = 0.75 - dy*dy;
        double srz = 0.75 - dz*dz;

        double lrx = 0.5*(1.5 - abs(dx))*(1.5 - abs(dx));
        double lry = 0.5*(1.5 - abs(dy))*(1.5 - abs(dy));
        double lrz = 0.5*(1.5 - abs(dz))*(1.5 - abs(dz));

        real_grid(i, j, k) += srx * sry * srz;

        double ux, uy, uz, lx, ly, lz;

        if (abs(dx) < 0.5) {
            ux = srx;
            lx = lrx;
        }
        else if (abs(dx) < 1.5) {
            ux = lrx;
            lx = srx;
        }
        else {
            ux = 0.;
            lx = 0.;
        }

        if (abs(dy) < 0.5) {
            uy = sry;
            ly = lry;
        }
        else if (abs(dy) < 1.5) {
            uy = lry;
            ly = sry;
        }
        else {
            uy = 0.;
            ly = 0.;
        }

        if (abs(dz) < 0.5) {
            uz = srz;
            lz = lrz;
        }
        else if (abs(dz) < 1.5) {
            uz = lrz;
            lz = srz;
        }
        else {
            uz = 0.;
            lz = 0.;
        }

        real_grid(i+1, j, k) += ux * sry * srz;
        real_grid(i, j+1, k) += srx * uy * srz;
        real_grid(i, j, k+1) += srx * sry * uz;
        real_grid(i+1, j+1, k) += ux * uy * srz;
        real_grid(i+1, j, k+1) += ux * sry * uz;
        real_grid(i, j+1, k+1) += srx * uy * uz;
        real_grid(i+1, j+1, k+1) += ux * uy * uz;


        real_grid(i-1, j, k) += lx * sry * srz;
        real_grid(i-1, j+1, k) += lx * uy * srz;
        real_grid(i-1, j, k+1) += lx * sry * uz;
        real_grid(i-1, j+1, k+1) += dx * dy * dz;


        real_grid(i, j-1, k) += srx * ly * srz;
        real_grid(i+1, j-1, k) += ux * ly * srz;
        real_grid(i, j-1, k+1) += srx * ly * uz;
        real_grid(i+1, j-1, k+1) += ux * ly * uz;


        real_grid(i, j, k-1) += srx * sry * lz;
        real_grid(i+1, j, k-1) += ux * sry * lz;
        real_grid(i, j+1, k-1) += srx * uy * lz;
        real_grid(i+1, j+1, k-1) += ux * uy * lz;


        real_grid(i-1, j-1, k) += lx * ly * srz;
        real_grid(i-1, j-1, k+1) += lx * ly * uz;


        real_grid(i-1, j, k-1) += lx * sry * lz;
        real_grid(i-1, j+1, k-1) += lx * uy * lz;


        real_grid(i, j-1, k-1) += srx * ly * lz;
        real_grid(i+1, j-1, k-1) += ux * ly * lz;


        real_grid(i-1, j-1, k-1) += lx * ly * lz;
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                real_grid(i, j, k) = real_grid(i, j, k) / avg_particle_density - 1;
            }
        }
    }
}