//
// Created by pierfied on 11/24/18.
//

#include <fstream>
#include "integrator.h"

void integrator::write_positions(double a) {
    std::string fname = file_prefix + "_" + std::to_string(write_num++) + ".bin";
    std::ofstream file(fname, std::ios::out | std::ios::binary);

    file.write((char *) &a, sizeof(double));
    file.write((char *) g->plist->x->data, sizeof(double) * g->plist->num_particles * g->plist->num_dims);

    file.close();
}
