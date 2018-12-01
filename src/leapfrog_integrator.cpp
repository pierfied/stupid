//
// Created by pierfied on 11/24/18.
//

#include "leapfrog_integrator.h"

void leapfrog_integrator::forward_half_step_p(double a) {
    g->compute_potential(a);

    for (int i = 0; i < g->plist->num_particles; ++i) {
        g->plist->p->index(i, 0) += g->cosmo.f(a) * g->particle_accel(i, 0) * delta_a / 2;
        g->plist->p->index(i, 1) += g->cosmo.f(a) * g->particle_accel(i, 1) * delta_a / 2;
        g->plist->p->index(i, 2) += g->cosmo.f(a) * g->particle_accel(i, 2) * delta_a / 2;
    }
}

void leapfrog_integrator::backward_half_step_p(double a) {
    g->compute_potential(a);

    for (int i = 0; i < g->plist->num_particles; ++i) {
        g->plist->p->index(i, 0) -= g->cosmo.f(a) * g->particle_accel(i, 0) * delta_a / 2;
        g->plist->p->index(i, 1) -= g->cosmo.f(a) * g->particle_accel(i, 1) * delta_a / 2;
        g->plist->p->index(i, 2) -= g->cosmo.f(a) * g->particle_accel(i, 2) * delta_a / 2;
    }
}

void leapfrog_integrator::full_step(double a) {
    g->compute_potential(a);

    for (int i = 0; i < g->plist->num_particles; ++i) {
        g->plist->p->index(i, 0) += g->cosmo.f(a) * g->particle_accel(i, 0) * delta_a;
        g->plist->p->index(i, 1) += g->cosmo.f(a) * g->particle_accel(i, 1) * delta_a;
        g->plist->p->index(i, 2) += g->cosmo.f(a) * g->particle_accel(i, 2) * delta_a;
    }

    a += 0.5 * delta_a;
    for (int i = 0; i < g->plist->num_particles; ++i) {
        g->plist->x->index(i, 0) += g->cosmo.f(a) * g->plist->p->index(i, 0) * delta_a / (a * a);
        g->plist->x->index(i, 1) += g->cosmo.f(a) * g->plist->p->index(i, 1) * delta_a / (a * a);
        g->plist->x->index(i, 2) += g->cosmo.f(a) * g->plist->p->index(i, 2) * delta_a / (a * a);
    }

    g->rebound_positions();
}

void leapfrog_integrator::run_sim() {
    double a = a0;

    if (write_pos) {
        write_pos_and_mom(a);
    }

    backward_half_step_p(a);

    for (int i = 0; i < num_steps; ++i) {
        full_step(a);

        a += delta_a;

        if (write_pos && (i + 1) % write_nth_step == 0) {
            write_pos_and_mom(a);
        }
    }

    forward_half_step_p(a);

    if (write_pos && num_steps % write_nth_step != 0) {
        write_pos_and_mom(a);
    }
}
