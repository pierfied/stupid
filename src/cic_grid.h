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

    inline double particle_accel(int particle_ind, int dim) override;
};


#endif //STUPID_CIC_GRID_H
