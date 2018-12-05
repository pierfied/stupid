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

    initial_conditions(particle_list &plist) : plist(&plist) {}

    virtual void setup_ics() = 0;

};

#endif //STUPID_INITIAL_CONDITIONS_H
