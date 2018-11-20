#include "stupid.h"
#include "multidim_array.h"

#include <iostream>

void hello() {
    std::cout << "Hello, World!" << std::endl;
}

int main() {
    hello();

    int num_particles = 10;
    int num_dims = 3;

    array_2d<double> x(num_particles, num_dims);
    x(0,0) = 1;
    x(1,1) = 2;

    std::cout << x(0,0) << "\t\t" << x(1,1) << "\t\t" << x(0,1) << std::endl;
}