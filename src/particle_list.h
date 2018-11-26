//
// Created by pierfied on 11/18/18.
//

#ifndef STUPID_PARTICLE_LIST_H
#define STUPID_PARTICLE_LIST_H

#include "multidim_array.h"

class particle_list {
public:
    const int num_particles;
    array_2d<double> *x;
    array_2d<double> *p;

    particle_list(array_2d<double> &x, array_2d<double> &p) : x(&x), p(&p), num_particles(x.nx) {}
};

#endif //STUPID_PARTICLE_LIST_H
