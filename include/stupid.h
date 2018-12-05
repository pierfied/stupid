#ifndef STUPID_LIBRARY_H
#define STUPID_LIBRARY_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int num_particles;
    double *x;
    double *p;
    int num_cells;
    double a0;
    double af;
    double delta_a;
    double Omega_m0;
    double Omega_k0;
    double Omega_l0;
    double sigma8;
    char *file_prefix;
    int interp_scheme;
    int integrator;
    int write_nth_step;
    bool use_real_units;
    double r0;
    double t0;
} stupid_args;

void run_sim(stupid_args args);

void hello();

#ifdef __cplusplus
}
#endif

#endif