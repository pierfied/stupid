//
// Created by pierfied on 11/19/18.
//

#ifndef STUPID_GRID_H
#define STUPID_GRID_H


#include <iostream>
#include "multidim_array.h"
#include "particle_list.h"

class grid {
public:
    const int nx;
    const int ny;
    const int nz;

    particle_list *plist;

    array_3d<double> delta_grid;

    grid(int nx, int ny, int nz, particle_list &plist) : nx(nx), ny(ny), nz(nz), delta_grid(nx, ny, nz), plist(&plist) {}

    virtual void populate_delta_grid() = 0;
};


#endif //STUPID_GRID_H
