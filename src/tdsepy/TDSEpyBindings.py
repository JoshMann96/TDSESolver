import ctypes
lib = ctypes.CDLL('./libTDSEpy.so')

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
    ctypes.c_void_p
        Pointer to simulation instance.
    """
    xmin = ctypes.c_double(xmin)
    xmax = ctypes.c_double(xmax)
    dx = ctypes.c_double(dx)
    dt = ctypes.c_double(dt)
    maxT = ctypes.c_double(maxT)
    return ctypes.c_void_p(lib.createSimulation(xmin, xmax, dx, dt, maxT))

lib.deleteSimulation.argtypes = [ctypes.c_void_p]
lib.deleteSimulation.restype = ctypes.c_int

def deleteSimulation(ptr:ctypes.c_void_p):
    """
    Deletes a TDSE simulation instance.

    Parameters
    ----------
    ptr : c_void_p
        Pointer to simulation instance.

    Returns
    -------
    int
        0 if successful.
    """
    return lib.deleteSimulation(ptr)

lib.addPot_FilePotential.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_wchar_p, ctypes.c_int]
lib.addPot_FilePotential.restype = None

def addPot_FilePotential(ptr:ctypes.c_void_p, offset:float, fil:str, refPoint:int):
    """
    Adds potential stored in file.

    Parameters
    ----------
    ptr : c_void_p
        Pointer to simulation instance.
    offset : float
        Positional offset.
    fil : str
        File containing potential.
    refPoint : int
        Maximum time in simulation (end time).
    """
    offset = ctypes.c_double(offset)
    fil = ctypes.c_wchar_p(fil)
    refPoint = ctypes.c_int(refPoint)
    return lib.addPot_FilePotential(ptr, offset, fil, refPoint)

lib.addPot_JelliumPotentialBacked.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int]
lib.addPot_JelliumPotentialBacked.restype = None

def addPot_JelliumPotentialBacked(ptr:ctypes.c_void_p, center:float, ef:float, w:float, backStart:float, backWidth:float, refPoint:int):
    """
    Adds Jellium potential.

    Parameters
    ----------
    ptr : c_void_p
        Pointer to simulation instance.
    center : float
        Jellium center (sigmoid center).
    ef : float
        Fermi energy.
    w : float
        Work function.
    backStart : float
        Internal start position of polynomial smooth spline backing.
    backWidth : float
        Width of polynomial smooth spline backing.
    refPoint : int
        Maximum time in simulation (end time).
    """
    center = ctypes.c_double(center)
    ef = ctypes.c_double(ef)
    w = ctypes.c_double(w)
    backStart = ctypes.c_double(backStart)
    backWidth = ctypes.c_double(backWidth)
    refPoint = ctypes.c_int(refPoint)
    return lib.addPot_JelliumPotentialBacked(ptr, center, ef, w, backStart, backWidth, refPoint)
