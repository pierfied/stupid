//
// Created by pierfied on 11/24/18.
//

#ifndef STUPID_LEAPFROG_INTEGRATOR_H
#define STUPID_LEAPFROG_INTEGRATOR_H


#include "integrator.h"
#include <math.h>

class leapfrog_integrator : public integrator {
    const long num_steps;
    const double delta_a;
    const double a0;
    const bool write_pos;
    const long write_nth_step;

public:
    leapfrog_integrator(double a0, long num_steps, double delta_a, grid &g) :
            integrator(g, std::string()), num_steps(num_steps), delta_a(delta_a), a0(a0), write_pos(false),
            write_nth_step(0) {}

    leapfrog_integrator(double a0, long num_steps, double delta_a, grid &g, std::string file_prefix, long write_nth_step)
            : integrator(g, file_prefix), num_steps(num_steps), delta_a(delta_a), a0(a0), write_pos(true),
              write_nth_step(write_nth_step) {}

    leapfrog_integrator(double a0, double a_final, double delta_a, grid &g) :
            integrator(g, std::string()), delta_a(delta_a), a0(a0), num_steps((int) ceil((a_final - a0) / delta_a)),
            write_pos(false), write_nth_step(0) {}

    leapfrog_integrator(double a0, double a_final, double delta_a, grid &g, std::string file_prefix, long write_nth_step)
            : integrator(g, file_prefix), delta_a(delta_a), a0(a0), num_steps((int) ceil((a_final - a0) / delta_a)),
              write_pos(true), write_nth_step(write_nth_step) {}

private:
    void forward_half_step_p(double a);

    void backward_half_step_p(double a);

    void full_step(double a);

public:
    void run_sim() override;
};


#endif //STUPID_LEAPFROG_INTEGRATOR_H
