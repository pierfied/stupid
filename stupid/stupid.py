import ctypes
import os

import numpy

class STUPID():
    def __init__(self):
        lib_path = os.path.join(os.path.dirname(__file__), 'libstupid.so')
        self.stupid_lib = ctypes.cdll.LoadLibrary(lib_path)