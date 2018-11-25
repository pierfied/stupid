//
// Created by pierfied on 11/24/18.
//

#ifndef STUPID_LEAPFROG_INTEGRATOR_H
#define STUPID_LEAPFROG_INTEGRATOR_H


#include "integrator.h"
#include <math.h>

class leapfrog_integrator : public integrator {
    const int num_steps;
    const double delta_a;
    const double a0;

public:
    leapfrog_integrator(double a0, int num_steps, double delta_a, grid &g) : integrator(g), num_steps(num_steps),
                                                                             delta_a(delta_a), a0(a0) {}

    leapfrog_integrator(double a0, double a_final, double delta_a, grid &g) : integrator(g), delta_a(delta_a), a0(a0),
                                                                              num_steps(
                                                                                      (int) ceil((a_final - a0) /
                                                                                                 delta_a)) {}

private:
    void forward_half_step_p(double a);

    void backward_half_step_p(double a);

    void full_step(double a);

public:
    void run_sim() override;
};


#endif //STUPID_LEAPFROG_INTEGRATOR_H
