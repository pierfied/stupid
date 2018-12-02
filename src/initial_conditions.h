//
// Created by rafa on 12/1/18.
//

#ifndef STUPID_INITIAL_CONDITIONS_H
#define STUPID_INITIAL_CONDITIONS_H

#include "cosmology.h"
#include "particle_list.h"

class initial_conditions {

public:
    cosmology *cosmo;
    particle_list *plist;

    initial_conditions(cosmology &cosmo, particle_list &plist) : cosmo(&cosmo), plist(&plist){}

    virtual void setup_ics() = 0;

};

#endif //STUPID_INITIAL_CONDITIONS_H
