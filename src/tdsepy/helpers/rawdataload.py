from io import BufferedReader
import json
import numpy as np
from typing import Literal
import os

_INT32_SIZE = np.dtype(np.int32).itemsize
_DOUBLE_SIZE = np.dtype(np.float64).itemsize
_C_DTYPES = Literal["int", "int32", "double", "char"]
_CONSTANT_NAMES = Literal["dx", "dt", "emax", "lam", "tau", "rad", "ef", "wf", "nElec", "nPts", "nSteps", "abs_rate", "abs_width"]
_CONSTANT_DTYPES = {
    "intSize" : "int32",
    "dx" : "double",
    "dt" : "double",
    "emax" : "double",
    "lam" : "double",
    "tau" : "double",
    "rad" : "double",
    "ef" : "double",
    "wf" : "double",
    "nElec" : "int",
    "nPts" : "int",
    "nSteps" : "int",
    "abs_rate" : "double",
    "abs_width" : "double"
}

def readData(fil:BufferedReader, dtype:_C_DTYPES, shape:int|tuple=1, INT_SIZE:int=None) -> np.ndarray|int|float|str:
    if shape is not tuple:
        shape = (shape)
        
    match dtype:
        case "int32":
            dat = np.array(np.fromfile(fil, np.int32, np.prod(shape)))
        case "double":
            dat = np.array(np.fromfile(fil, np.float64, np.prod(shape)))
        case "char":
            dat = np.array(np.fromfile(fil, np.ubyte, np.prod(shape)))
            dat = ''.join([chr(it) for it in dat])
        case "complex":
            dat = np.array(np.fromfile(fil, np.double, np.prod(shape)*2))
            dat = dat[0::2] + 1.0j*dat[1::2]
        case "int":
            match INT_SIZE:
                case 4:
                    dat = np.array(np.fromfile(fil, np.int32, np.prod(shape)))
                case 8:
                    dat = np.array(np.fromfile(fil, np.int64, np.prod(shape)))
                case _:
                    raise ValueError("INT_SIZE must be 4 or 8 for dtype 'int'")
        case _:
            raise ValueError("dtype must be 'int32', 'double', 'char', 'complex' or 'int'")
    
    if np.prod(shape) == 1:
        dat = dat[0]
    elif dtype != "char":
        if len(dat) == np.prod(shape):
            dat = np.reshape(dat, shape)
        elif len(dat) == 0:
            raise ValueError("No data read.")
        else:
            numel = (len(dat) // np.prod(shape[1:])) * np.prod(shape[1:])
            dat = np.reshape(dat[:numel], (-1,) + (shape[1:]))
        
    return dat

def getConstant(name:_CONSTANT_NAMES, fol:str, dtype:_C_DTYPES = None) -> tuple[int|float, int]:
    """
    Reads a constant (just a single value) from a .dat file.
    Args:
        name (_CONSTANT_NAMES): Name of the constant to read.
        fol (str): Folder where the .dat file is located.
        dtype (_C_DTYPES, optional): Data type to read. If None, uses default for the constant. Defaults to None.
    Returns:
        value: The value of the constant.
        typ (int): Index identifier of the measurer type.
    """
    with open(os.path.join(fol, name+".dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        dat = readData(fil, _CONSTANT_DTYPES[name] if dtype is None else dtype, INT_SIZE=INT_SIZE)
    return dat, typ

def getPsi2t(fol:str) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Reads the probability density vs time data from a psi2t.dat file.
    Args:
        fol (str): Folder where the psi2t.dat file is located.
    Returns:
        dat (np.ndarray): Probability density data with shape (nElec, nt, nx).
        xs (np.ndarray): Spatial grid points.
        ts (np.ndarray): Temporal grid points.
        typ (int): Index identifier of the measurer type.
    """
    nElec,_ = getConstant("nElec", fol)
    with open(os.path.join(fol, "psi2t.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        nx = readData(fil, "int", INT_SIZE=INT_SIZE)
        nt = readData(fil, "int", INT_SIZE=INT_SIZE)
        dat = readData(fil, "double", (nt,nElec,nx)).swapaxes(0,1)
        try:
            xs = readData(fil, "double", nx)
            ts = readData(fil, "double", nt)
        except ValueError:
            xs = np.linspace(0,dat.shape[-1] / nx, dat.shape[-1])
            ts = np.linspace(0,dat.shape[-2] / nt, dat.shape[-2])
    return dat, xs, ts, typ

def getVfunct(fol:str, index:int = -1) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Reads the potential vs time data from a Vfunct.dat file.
    Args:
        fol (str): Folder where the Vfunct.dat file is located.
        index (int, optional): Index of the potential to read. If -1, reads Vfunct.dat. Defaults to -1.
    Returns:
        dat (np.ndarray): Potential data with shape (nt, nx).
        xs (np.ndarray): Spatial grid points.
        ts (np.ndarray): Temporal grid points.
        typ (int): Index identifier of the measurer type.
    """
    with open(os.path.join(fol, "Vfunct.dat" if index < 0 else f"{index:d}Vfunct.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        nx = readData(fil, "int", INT_SIZE=INT_SIZE)
        nt = readData(fil, "int", INT_SIZE=INT_SIZE)
        dat = readData(fil, "double", (nt, nx))
        try:
            xs = readData(fil, "double", nx)
            ts = readData(fil, "double", nt)
        except ValueError:
            xs = np.linspace(0,dat.shape[-1] / nx, dat.shape[-1])
            ts = np.linspace(0,dat.shape[-2] / nt, dat.shape[-2])
    return dat, xs, ts, typ

def getWghts(fol:str) -> tuple[np.ndarray, int]:
    """
    Reads the orbital weights from a wghts.dat file.
    Args:
        fol (str): Folder where the wghts.dat file is located.
    Returns:
        wghts (np.ndarray): Weights of the orbitals with shape (nElec,).
        typ (int): Index identifier of the measurer type.
    """
    with open(os.path.join(fol, "wghts.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        nElec = readData(fil, "int", INT_SIZE=INT_SIZE)
        wghts = readData(fil, "double", nElec)
    return wghts, typ

def getFluxSpecVD(fol:str, vdNum:int=0) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, str, int]:
    """
    Reads the bidirectional flux spectrum data from a fluxspecvd.dat file.
    Args:
        fol (str): Folder where the fluxspecvd.dat file is located.
        vdNum (int, optional): Index of the virtual detector. Defaults to 0.
    Returns:
        energies (np.ndarray): Signed energy (negative for left-moving, positive for right-moving) grid points with shape (nSamp*2-1).
        momenta (np.ndarray): Momentum (wavenumber) grid points with shape (nSamp*2-1).
        psik (np.ndarray): Wavefunction in signed energy space with shape (nElec, nSamp*2-1). 
            The normalization is chosen such that the integral of |psik|^2 over the energies provides the total probability flux.
            Therefore, |psik|^2 = dP/dE.
        posIdx (int): Position index of the virtual detector.
        name (str): Name of the virtual detector.
        typ (int): Index identifier of the measurer type.
    """
    nElec,_ = getConstant("nElec", fol)
    with open(os.path.join(fol, f"{vdNum:d}" + "fluxspecvd.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        readData(fil, "int32") #skip VD index
        name = readData(fil, "char", 4)
        posIdx = readData(fil, "int", INT_SIZE=INT_SIZE)
        nSamp = readData(fil, "int", INT_SIZE=INT_SIZE)
        energies = readData(fil, "double", (nSamp*2-1))
        momenta = readData(fil, "double", (nSamp*2-1))
        psik = readData(fil, "complex", (nElec, nSamp*2-1))
    return energies, momenta, psik, posIdx, name, typ

def getUnidirectionalFluxSpecVD(fol:str, vdNum:int=0) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, str, int]:
    """
    Reads the unidirectional flux spectrum data from a unifluxspecvd.dat file.
    Args:
        fol (str): Folder where the unifluxspecvd.dat file is located.
        vdNum (int, optional): Index of the virtual detector. Defaults to 0.
    Returns:
        momenta (np.ndarray): Momentum (wavenumber) grid points with shape (nSamp,).
        energies (np.ndarray): Energy grid points with shape (nSamp,).
        psift (np.ndarray): Wavefunction in energy space with shape (nElec, nSamp). 
            The normalization is chosen such that the integral of |psift|^2 over the energies provides the total probability flux.
            Therefore, |psift|^2 = dP/dE.
        posIdx (int): Position index of the virtual detector.
        name (str): Name of the virtual detector.
        typ (int): Index identifier of the measurer type.
    """
    nElec,_ = getConstant("nElec", fol)
    with open(os.path.join(fol, f"{vdNum:d}" + "unifluxspecvd.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        readData(fil, "int32") #skip VD index
        name = readData(fil, "char", 4)
        posIdx = readData(fil, "int", INT_SIZE=INT_SIZE)
        nSamp = readData(fil, "int", INT_SIZE=INT_SIZE)
        energies = readData(fil, "double", (nSamp))
        momenta = readData(fil, "double", (nSamp))
        psift = readData(fil, "complex", (nElec, nSamp))
    return energies, momenta, psift, posIdx, name, typ

def getClassicalSpecVD(fol:str, vdNum:int=0) -> tuple[np.ndarray, np.ndarray, int, str, int]:
    """
    Reads the classical flux spectrum data from a classicalFlux.dat file.
    Args:
        fol (str): Folder where the classicalFlux.dat file is located.
        vdNum (int, optional): Index of the virtual detector. Defaults to 0.
    Returns:
        momenta (np.ndarray): Momentum (wavenumber) grid points with shape (nSamp*2-1).
        yields (np.ndarray): Yield data with shape (nElec, nSamp*2-1). 
            The normalization is chosen such that the integral of yields over the momenta provides the total probability flux.
            Thus, ylds = dP/dk.
        posIdx (int): Position index of the virtual detector.
        name (str): Name of the virtual detector.
        typ (int): Index identifier of the measurer type.
    """
    nElec,_ = getConstant("nElec", fol)
    with open(os.path.join(fol, f"{vdNum:d}" + "classicalfluxspecvd.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        readData(fil, "int32") #skip VD index
        name = readData(fil, "char", 4)
        posIdx = readData(fil, "int", INT_SIZE=INT_SIZE)
        nSamp = readData(fil, "int", INT_SIZE=INT_SIZE)
        momenta = readData(fil, "double", (nSamp*2-1))
        yields = readData(fil, "double", (nElec, nSamp*2-1))
    return momenta, yields, posIdx, name, typ

def getExpectE0(fol:str) -> tuple[np.ndarray, int]:
    """
    Reads the initial expectation values of energy from an expectE0.dat file.
    Args:
        fol (str): Folder where the expectE0.dat file is located.
    Returns:
        e0 (np.ndarray): Initial expectation values of energy with shape (nElec,).
        typ (int): Index identifier of the measurer type.
    """
    nElec,_ = getConstant("nElec", fol)
    with open(os.path.join(fol, "expectE0.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        e0 = readData(fil, "double", nElec)
    return e0, typ

def getTs(fol:str) -> tuple[np.ndarray, int]:
    """
    Reads the time value for each time step from a ts.dat file.
    Args:
        fol (str): Folder where the ts.dat file is located.
    Returns:
        ts (np.ndarray): Temporal grid points with shape (nt,).
        typ (int): Index identifier of the measurer type.
    """
    nt,_ = getConstant("nSteps", fol)
    try:
        with open(os.path.join(fol, "ts.dat"), 'rb') as fil:
            INT_SIZE = readData(fil, 'int32')
            typ = readData(fil, 'int32')
            ts = readData(fil, "double", nt)
    except FileNotFoundError:
        dt,_ = getConstant("dt", fol)
        ts = np.arange(nt)*dt
        typ = 0
    return ts, typ

SIMPLE_ARRAY_PERWF_QUANTITIES = Literal["expectA", "expectE", "expectX", "expectP", "totProb"]

def getSimpleArrayPerWFData(fol:str, quantity:SIMPLE_ARRAY_PERWF_QUANTITIES) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Reads simple array-per-wavefunction data from a .dat file.
    Intended for internal use, see getExpect* and getTotProb.
    Args:
        fol (str): Folder where the .dat file is located.
        quantity (SIMPLE_ARRAY_PERWF_QUANTITIES): Name of the quantity to read.
    Returns:
        ts (np.ndarray): Temporal grid points with shape (nt,).
        dat (np.ndarray): Data array with shape (nt, nElec).
        typ (int): Index identifier of the measurer type.
    """
    nElec,_ = getConstant("nElec", fol)
    ts,_ = getTs(fol)
    with open(os.path.join(fol, f"{quantity}.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        dat = readData(fil, "double", (len(ts), nElec))
    return ts, dat, typ

def getExpectA(fol:str):
    """
    Reads the expectation values of acceleration from an expectA.dat file.
    Args:
        fol (str): Folder where the expectA.dat file is located.
    Returns:
        ts (np.ndarray): Temporal grid points with shape (nt,).
        dat (np.ndarray): Expectation values of acceleration with shape (nt, nElec).
        typ (int): Index identifier of the measurer type.
    """
    return getSimpleArrayPerWFData(fol, "expectA")
def getExpectX(fol:str):
    """
    Reads the expectation values of position from an expectX.dat file.
    Args:
        fol (str): Folder where the expectX.dat file is located.
    Returns:
        ts (np.ndarray): Temporal grid points with shape (nt,).
        dat (np.ndarray): Expectation values of position with shape (nt, nElec).
        typ (int): Index identifier of the measurer type.
    """
    return getSimpleArrayPerWFData(fol, "expectX")
def getExpectE(fol:str):
    """
    Reads the expectation values of energy from an expectE.dat file.
    Args:
        fol (str): Folder where the expectE.dat file is located.
    Returns:
        ts (np.ndarray): Temporal grid points with shape (nt,).
        dat (np.ndarray): Expectation values of energy with shape (nt, nElec).
        typ (int): Index identifier of the measurer type.
    """
    return getSimpleArrayPerWFData(fol, "expectE")
def getExpectP(fol:str):
    """
    Reads the expectation values of momentum from an expectP.dat file.
    Args:
        fol (str): Folder where the expectP.dat file is located.
    Returns:
        ts (np.ndarray): Temporal grid points with shape (nt,).
        dat (np.ndarray): Expectation values of momentum with shape (nt, nElec).
        typ (int): Index identifier of the measurer type.
    """
    return getSimpleArrayPerWFData(fol, "expectP")
def getTotProb(fol:str):
    """
    Reads the total probability from a totProb.dat file.
    Args:
        fol (str): Folder where the totProb.dat file is located.
    Returns:
        ts (np.ndarray): Temporal grid points with shape (nt,).
        dat (np.ndarray): Total probability with shape (nt, nElec).
        typ (int): Index identifier of the measurer type.
    """
    return getSimpleArrayPerWFData(fol, "totProb")