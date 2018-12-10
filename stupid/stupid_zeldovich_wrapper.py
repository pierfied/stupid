from ctypes import *
import os
import numpy as np
from ctypes.util import find_library


class _Args_ZAICs_Struct(Structure):
    _fields_ = [
        ('x', POINTER(c_double)),
        ('p', POINTER(c_double)),
        ('num_particles', c_long),
        ('num_cells', c_long),
        ('sizeofk', c_long),
        ('k', POINTER(c_double)),
        ('P', POINTER(c_double)),
        ('a0', c_double),
        ('Omega_m0', c_double),
        ('Omega_k0', c_double),
        ('Omega_l0', c_double),
        ('H0', c_double),
    ]


class STUPID_ICS:
    def __init__(self, x, p, k, P, a0, delta_a, num_cells, cosmo):
        lib_path = os.path.join(os.path.dirname(__file__), 'libstupid.so')
        self.stupid_lib = cdll.LoadLibrary(lib_path)

        self.stupid_ics = self.stupid_lib.setup_ics
        self.stupid_ics.argtypes = [_Args_ZAICs_Struct]
        self.stupid_ics.restype = None

        self.x = np.ascontiguousarray(x, dtype=np.double)
        self.p = np.ascontiguousarray(p, dtype=np.double)
        self.k = np.ascontiguousarray(k, dtype=np.double)
        self.P = np.ascontiguousarray(P, dtype=np.double)
        self.a0 = a0
        self.delta_a = delta_a
        self.num_cells = num_cells
        self.cosmo = cosmo

    def setup_ics(self):
        args = _Args_ZAICs_Struct()
        args.num_particles = self.x.shape[0]
        args.x = self.x.ctypes.data_as(POINTER(c_double))
        args.p = self.p.ctypes.data_as(POINTER(c_double))
        args.num_cells = self.num_cells
        args.sizeofk = self.k.size
        args.k = self.k.ctypes.data_as(POINTER(c_double))
        args.P = self.P.ctypes.data_as(POINTER(c_double))
        args.a0 = self.a0
        args.delta_a = self.delta_a
        args.Omega_m0 = self.cosmo.Omega_m0
        args.Omega_k0 = self.cosmo.Omega_k0
        args.Omega_l0 = self.cosmo.Omega_l0
        args.H0 = self.cosmo.H0
        print(self.cosmo.H0)

        print('Set up initial conditions')

        self.stupid_ics(args)

        SIZE = self.x.shape[0]*3
        libc = CDLL(find_library('c'))
        libc.malloc.restype = c_void_p
        x_pointer = cast(args.x,POINTER(c_double))
        p_pointer = cast(args.p,POINTER(c_double))
        self.new_x = np.ctypeslib.as_array(x_pointer,shape=(SIZE,))
        self.new_p = np.ctypeslib.as_array(p_pointer,shape=(SIZE,))

        self.new_x = np.reshape(self.new_x, self.x.shape)
        self.new_p = np.reshape(self.new_p, self.p.shape)

        print('Finished setting up ICs')
