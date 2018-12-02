from ctypes import *
import os
import numpy as np


class _Args_Struct(Structure):
    _fields_ = [
        ('num_particles', c_int),
        ('x', POINTER(c_double)),
        ('p', POINTER(c_double)),
        ('num_cells', c_double),
        ('a0', c_double),
        ('af', c_double),
        ('delta_a', c_double),
        ('Omega_m0', c_double),
        ('Omega_k0', c_double),
        ('Omega_l0', c_double),
        ('file_prefix', c_char_p),
        ('interp_scheme', c_int),
        ('integrator', c_int),
        ('write_nth_step', c_int)
    ]


class Cosmology:
    def __init__(self, Omega_m0, Omega_k0, Omega_l0):
        self.Omega_m0 = Omega_m0
        self.Omega_k0 = Omega_k0
        self.Omega_l0 = Omega_l0


class STUPID:
    def __init__(self, x, p, a0, af, delta_a, num_cells, cosmo, file_prefix, interp_scheme, integrator,
                 write_nth_step):
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
        self.interp_scheme = interp_scheme
        self.integrator = integrator
        self.write_nth_step = write_nth_step

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
        args.file_prefix = c_char_p(self.file_prefix.encode('utf-8'))
        args.interp_scheme = self.interp_scheme
        args.integrator = self.integrator
        args.write_nth_step = self.write_nth_step

        print('Starting Run')

        self.stupid_run(args)

        print('Run Finished')
