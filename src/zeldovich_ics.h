//
// Created by rafa on 12/1/18.
//

#ifndef STUPID_ZELDOVICH_ICS_H
#define STUPID_ZELDOVICH_ICS_H

#include "initial_conditions.h"

class zeldovich_ics: public initial_conditions {

private:
    void generate_uniform_distribution();

    void fourier_displacement_vec();

    void real_displacement_vec();

    void apply_ZA();

public:
    void setup_ics() override;
};

#endif //STUPID_ZELDOVICH_ICS_H
