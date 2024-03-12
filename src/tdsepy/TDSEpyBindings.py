import ctypes
# Load shared library
lib = ctypes.CDLL('./libTDSEpy.so')


lib.test.argtypes = [ctypes.c_int]
lib.test.restype = None

def test(x):
    x = ctypes.c_int(x)
    return lib.test(x)

lib.test2.argtypes = [ctypes.c_int, ctypes.c_double]
lib.test2.restype = ctypes.c_double

def test2(x, y):
    x = ctypes.c_int(x)
    y = ctypes.c_double(y)
    return lib.test2(x, y)

lib.nottest.argtypes = []
lib.nottest.restype = ctypes.c_int

def nottest():
    return lib.nottest()
