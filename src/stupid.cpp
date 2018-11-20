#include "stupid.h"
#include "multidim_array.h"
#include "particle_list.h"
#include "cic_grid.h"

#include <iostream>

void hello() {
    std::cout << "Hello, World!" << std::endl;
}

int main() {
    hello();

    int num_particles = 10;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    x(0, 0) = 1;
    x(1, 1) = 2;

    array_2d<double> p(num_particles, num_dims);

    particle_list plist(x, p);

    std::cout << plist.x(0, 0) << "\t\t" << plist.x(1, 1) << std::endl;

    cic_grid a(2, 2, 2, plist);

    x(0, 0) = -1;

    std::cout << a.plist.x(0,0) << std::endl;
}