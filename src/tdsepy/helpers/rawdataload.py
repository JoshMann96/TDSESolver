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

def readData(fil:BufferedReader, dtype:_C_DTYPES, shape:int|tuple=1, INT_SIZE:int=None):
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

def combinePath(fol:str, fil:str):
    """Combines path components, ensuring '/' between fol and fil.

    Args:
        fol (str): Folder, with or without final '/'.
        fil (str): File name without slashes.

    Returns:
        str: Final path.
    """    
    if '/' in fil:
        raise ValueError("fil must not be a path (contains '/')")
    if fol[-1] != '/':
        fol += '/'
    return fol + fil

def getConstant(name:_CONSTANT_NAMES, fol:str, dtype:_C_DTYPES = None):
    with open(combinePath(fol, name + ".dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        dat = readData(fil, _CONSTANT_DTYPES[name] if dtype is None else dtype, INT_SIZE=INT_SIZE)
    return dat, typ

def getPsi2t(fol:str):
    nElec,_ = getConstant("nElec", fol)
    with open(combinePath(fol, "psi2t.dat"), 'rb') as fil:
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

def getVfunct(fol:str, index:int = -1):
    with open(combinePath(fol, "Vfunct.dat" if index < 0 else f"{index:d}Vfunct.dat"), 'rb') as fil:
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

def getWghts(fol:str):
    with open(combinePath(fol, "wghts.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        nElec = readData(fil, "int", INT_SIZE=INT_SIZE)
        wghts = readData(fil, "double", nElec)
    return wghts, typ

def getFluxSpecVD(fol:str, vdNum:int=0):
    nElec,_ = getConstant("nElec", fol)
    with open(combinePath(fol, f"{vdNum:d}" + "fluxspecvd.dat"), 'rb') as fil:
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

def getExpectE0(fol:str):
    nElec,_ = getConstant("nElec", fol)
    with open(combinePath(fol, "expectE0.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        e0 = readData(fil, "double", nElec)
    return e0, typ

def getTs(fol:str):
    nt,_ = getConstant("nSteps", fol)
    try:
        with open(combinePath(fol, "ts.dat"), 'rb') as fil:
            INT_SIZE = readData(fil, 'int32')
            typ = readData(fil, 'int32')
            ts = readData(fil, "double", nt)
    except FileNotFoundError:
        dt,_ = getConstant("dt", fol)
        ts = np.arange(nt)*dt
        typ = 0
    return ts, typ

SIMPLE_ARRAY_PERWF_QUANTITIES = Literal["expectA", "expectE", "expectX", "expectP", "totProb"]

def getSimpleArrayPerWFData(fol:str, quantity:SIMPLE_ARRAY_PERWF_QUANTITIES):
    nElec,_ = getConstant("nElec", fol)
    ts,_ = getTs(fol)
    with open(combinePath(fol, f"{quantity}.dat"), 'rb') as fil:
        INT_SIZE = readData(fil, 'int32')
        typ = readData(fil, 'int32')
        dat = readData(fil, "double", (len(ts), nElec))
    return ts, dat, typ

def getExpectA(fol:str):
    return getSimpleArrayPerWFData(fol, "expectA")
def getExpectX(fol:str):
    return getSimpleArrayPerWFData(fol, "expectX")
def getExpectE(fol:str):
    return getSimpleArrayPerWFData(fol, "expectE")
def getExpectP(fol:str):
    return getSimpleArrayPerWFData(fol, "expectP")
def getTotProb(fol:str):
    return getSimpleArrayPerWFData(fol, "totProb")