import ctypes
# Load shared library
lib = ctypes.CDLL('./libTDSEpy.so')

#inst = ctypes.c_void_p(SimulationManager()) #NEED TO ADD PARAMETERS

lib.createSimulation.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
lib.createSimulation.restype = ctypes.c_void_p

def createSimulation(xmin:float, xmax:float, dx:float, dt:float, maxT:float):
    """
    Creates a TDSE Simulation object.

    Parameters
    ----------
    xmin : float
        Left boundary position.
    xmax : float
        Right boundary position.
    dx : float
        Spatial step size.
    dt : float
        Temporal step size.
    maxT : float
        Maximum time in simulation (end time).

    Returns
    -------
    int
        Pointer to simulation instance.
    """
    xmin = ctypes.c_double(xmin)
    xmax = ctypes.c_double(xmax)
    dx = ctypes.c_double(dx)
    dt = ctypes.c_double(dt)
    maxT = ctypes.c_double(maxT)
    return lib.createSimulation(xmin, xmax, dx, dt, maxT)

lib.deleteSimulation.argtypes = [ctypes.c_void_p]
lib.deleteSimulation.restype = ctypes.c_int

def deleteSimulation(ptr:ctypes.c_void_p):
    """
    Deletes a TDSE simulation instance.

    Parameters
    ----------
    ptr : int
        Pointer to simulation instance.

    Returns
    -------
    int
        0 if successful.
    """
    ptr = ctypes.c_void_p(ptr)
    return lib.deleteSimulation(ptr)

lib.addPot_FilePotential.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_wchar_p, ctypes.c_int]
lib.addPot_FilePotential.restype = ctypes.c_void_p

def addPot_FilePotential(ptr:ctypes.c_void_p, offset:float, fil:str, refPoint:int):
    """
    Adds potential stored in file.

    Parameters
    ----------
    ptr : int
        Pointer to simulation instance.
    offset : float
        Positional offset.
    fil : str
        File containing potential.
    refPoint : int
        Maximum time in simulation (end time).

    Returns
    -------
    int
        Pointer to simulation instance.
    """
    ptr = ctypes.c_void_p(ptr)
    offset = ctypes.c_double(offset)
    fil = ctypes.c_wchar_p(fil)
    refPoint = ctypes.c_int(refPoint)
    return lib.addPot_FilePotential(ptr, offset, fil, refPoint)
