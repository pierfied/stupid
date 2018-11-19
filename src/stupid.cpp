#include "stupid.h"
#include "particle_list.h"

#include <iostream>

void hello() {
    std::cout << "Hello, World!" << std::endl;
}

int main() {
    hello();

    int num_particles = 1000;

    double x[num_particles][NUM_DIMS];
    double p[num_particles][NUM_DIMS];

    x[0][2] = 1;

    particle_list plist(num_particles, x, p);

    std::cout << plist.x[0][2] << "\t\t" << plist.x[0][0] << std::endl;
}