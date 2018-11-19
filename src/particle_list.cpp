//
// Created by pierfied on 11/18/18.
//

#include "particle_list.h"

template<int _num_particles, int _num_dims>
particle_list<_num_particles, _num_dims>::particle_list(double (&x)[_num_particles][_num_dims], double (&p)[_num_particles][_num_dims]) {
    particle_list::x = &x;
    particle_list::p = &p;
}