//
// Created by pierfied on 11/19/18.
//

#ifndef STUPID_GRID_H
#define STUPID_GRID_H


#include "multidim_array.h"

class grid {
public:
    const int nx;
    const int ny;
    const int nz;

    array_3d<double> mass_grid;

    grid(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz), mass_grid(nx, ny, nz) {}

    virtual void populate_mass_grid() = 0;
};


#endif //STUPID_GRID_H
