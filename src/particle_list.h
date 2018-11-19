//
// Created by pierfied on 11/18/18.
//

#ifndef STUPID_PARTICLE_LIST_H
#define STUPID_PARTICLE_LIST_H

#define NUM_DIMS 3

class particle_list {
public:
    const int num_particles;
    double (*x)[NUM_DIMS];
    double (*p)[NUM_DIMS];

    particle_list(int num_particles, double x[][NUM_DIMS], double p[][NUM_DIMS]);
};

#endif //STUPID_PARTICLE_LIST_H
