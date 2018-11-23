//
// Created by pierfied on 11/20/18.
//

#ifndef STUPID_TSC_GRID_H
#define STUPID_TSC_GRID_H


#include "grid.h"

class tsc_grid : public grid {
public:
    tsc_grid(int nx, int ny, int nz, particle_list &plist) : grid(nx, ny, nz, plist) {}

    void populate_delta_grid() override;
};


#endif //STUPID_TSC_GRID_H
