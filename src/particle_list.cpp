//
// Created by pierfied on 11/18/18.
//

#include "particle_list.h"

particle_list::particle_list(int num_particles, double x[][NUM_DIMS], double p[][NUM_DIMS]) : num_particles(num_particles), x(x),
                                                                               p(p) {}