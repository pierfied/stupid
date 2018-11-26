//
// Created by pierfied on 11/24/18.
//

#include <fstream>
#include "integrator.h"

#define NUMDIMS 3

void integrator::write_positions() {
    std::string fname = file_prefix + std::to_string(write_num++);
    std::ofstream file(fname, std::ios::out | std::ios::binary);

    file.write((char *) g->plist->x->data, sizeof(double) * g->plist->num_particles * NUMDIMS);

    file.close();
}
