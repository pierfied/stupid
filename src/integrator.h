//
// Created by pierfied on 11/24/18.
//

#ifndef STUPID_INTEGRATOR_H
#define STUPID_INTEGRATOR_H


#include "grid.h"

class integrator {
    grid *g;

public:
    integrator(grid &g) : g(&g) {}

    virtual void run_sim() = 0;
};


#endif //STUPID_INTEGRATOR_H
