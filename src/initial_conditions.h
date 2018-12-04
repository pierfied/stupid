//
// Created by rafa on 12/1/18.
//

#ifndef STUPID_INITIAL_CONDITIONS_H
#define STUPID_INITIAL_CONDITIONS_H

#include "cosmology.h"
#include "particle_list.h"
#include "grid.h"

class initial_conditions {
public:
    particle_list *plist;

    grid *g;

    initial_conditions(particle_list &plist, grid &g) : plist(&plist), g(&g){}

    virtual void setup_ics() = 0;

};

#endif //STUPID_INITIAL_CONDITIONS_H
