from ctypes import *
import os
import numpy as np


class _Args_Struct(Structure):
    _fields_ = [
        ('num_particles', c_int),
        ('x', POINTER(c_double)),
        ('p', POINTER(c_double)),
        ('num_cells', c_int),
        ('a0', c_double),
        ('af', c_double),
        ('delta_a', c_double),
        ('Omega_m0', c_double),
        ('Omega_k0', c_double),
        ('Omega_l0', c_double),
        ('sigma8', c_double),
        ('file_prefix', c_char_p),
        ('interp_scheme', c_int),
        ('integrator', c_int),
        ('write_nth_step', c_int),
        ('use_real_units', c_bool),
        ('r0', c_double),
        ('t0', c_double)
    ]


class Cosmology:
    def __init__(self, Omega_m0, Omega_k0, Omega_l0, sigma8=0, H0=0):
        self.Omega_m0 = Omega_m0
        self.Omega_k0 = Omega_k0
        self.Omega_l0 = Omega_l0
        self.sigma8 = sigma8
        self.H0 = H0


class STUPID:
    def __init__(self, r, v, a0, af, delta_a, num_cells, cosmo, file_prefix, interp_scheme, integrator,
                 write_nth_step, use_real_units=True, box_len=None):
        if use_real_units:
            self.r0 = box_len / num_cells
            self.t0 = 1 / cosmo.H0
            v0 = self.r0 / self.t0

            x = r / (self.r0 * a0)
            p = v * a0 / v0
        else:
            self.r0 = 0
            self.t0 = 0

            x = r
            p = v
            self.r0 = 0
            self.t0 = 0

        lib_path = os.path.join(os.path.dirname(__file__), 'libstupid.so')
        self.stupid_lib = cdll.LoadLibrary(lib_path)

        self.stupid_run = self.stupid_lib.run_sim
        self.stupid_run.argtypes = [_Args_Struct]
        self.stupid_run.restype = None

        self.x = np.ascontiguousarray(x, dtype=np.double)
        self.p = np.ascontiguousarray(p, dtype=np.double)
        self.a0 = a0
        self.af = af
        self.delta_a = delta_a
        self.num_cells = num_cells
        self.cosmo = cosmo
        self.file_prefix = file_prefix
        self.write_nth_step = write_nth_step
        self.use_real_units = use_real_units

        if interp_scheme == 'cic':
            self.interp_scheme = 1
        elif interp_scheme == 'tsc':
            self.interp_scheme = 2
        else:
            self.interp_scheme = 0

        if integrator == 'leapfrog':
            self.integrator = 1
        else:
            self.integrator = 0

    def run_sim(self):
        args = _Args_Struct()
        args.num_particles = self.x.shape[0]
        args.x = self.x.ctypes.data_as(POINTER(c_double))
        args.p = self.p.ctypes.data_as(POINTER(c_double))
        args.num_cells = self.num_cells
        args.a0 = self.a0
        args.af = self.af
        args.delta_a = self.delta_a
        args.Omega_m0 = self.cosmo.Omega_m0
        args.Omega_k0 = self.cosmo.Omega_k0
        args.Omega_l0 = self.cosmo.Omega_l0
        args.sigma8 = self.cosmo.sigma8
        args.file_prefix = c_char_p(self.file_prefix.encode('utf-8'))
        args.interp_scheme = self.interp_scheme
        args.integrator = self.integrator
        args.write_nth_step = self.write_nth_step
        args.use_real_units = self.use_real_units
        args.r0 = self.r0
        args.t0 = self.t0

        print('Starting Run')

        self.stupid_run(args)

        print('Run Finished')
