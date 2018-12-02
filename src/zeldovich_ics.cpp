//
// Created by rafa on 12/1/18.
//

#include "zeldovich_ics.h"

void zeldovich_ics::generate_uniform_distribution() {}

void zeldovich_ics::fourier_displacement_vec() {}

void zeldovich_ics::real_displacement_vec() {}

void zeldovich_ics::apply_ZA() {}

void zeldovich_ics::setup_ics() {
    generate_uniform_distribution();

    fourier_displacement_vec();

    real_displacement_vec();

    apply_ZA();
}