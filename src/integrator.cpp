//
// Created by pierfied on 11/24/18.
//

#include <fstream>
#include "integrator.h"

void integrator::write_positions(double a) {
    std::string fname = file_prefix + "_x_" + std::to_string(write_num) + ".bin";
    std::ofstream file(fname, std::ios::out | std::ios::binary);

    file.write((char *) &a, sizeof(double));
    file.write((char *) g->plist->x->data, sizeof(double) * g->plist->num_particles * g->plist->num_dims);

    file.close();
}

void integrator::write_momentum(double a) {
    std::string fname = file_prefix + "_p_" + std::to_string(write_num) + ".bin";
    std::ofstream file(fname, std::ios::out | std::ios::binary);

    file.write((char *) &a, sizeof(double));
    file.write((char *) g->plist->p->data, sizeof(double) * g->plist->num_particles * g->plist->num_dims);

    file.close();
}

void integrator::write_pos_and_mom(double a) {
    write_positions(a);
    write_momentum(a);

    write_num++;
}
