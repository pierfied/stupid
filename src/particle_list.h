//
// Created by pierfied on 11/18/18.
//

#ifndef STUPID_PARTICLE_LIST_H
#define STUPID_PARTICLE_LIST_H

template<int _num_particles, int _num_dims>
class particle_list {
public:
    static const int num_particles = _num_particles;
    static const int num_dims = _num_dims;
    double (*x)[num_particles][num_dims];
    double (*p)[num_particles][num_dims];

    particle_list(double (&x)[_num_particles][_num_dims], double (&p)[_num_particles][_num_dims]);
};

#endif //STUPID_PARTICLE_LIST_H
