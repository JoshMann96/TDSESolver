from io import BufferedReader
import json
import numpy as np
from typing import Literal
import os

_INT_SIZE = np.dtype(np.int32).itemsize
_DOUBLE_SIZE = np.dtype(np.float64).itemsize
_CONSTANT_NAMES = Literal["dx", "dt", "emax", "lam", "tau", "rad", "ef", "wf", "nElec", "nPts", "nSteps"]
_CONSTANT_DTYPES = {
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
    "nSteps" : "int"
}

def readData(fil:BufferedReader, dtype:Literal["int", "double", "char"], shape=1):
    if shape is not tuple:
        shape = (shape)
        
    match dtype:
        case "int":
            dat = np.array(np.fromfile(fil, np.int32, np.prod(shape)))
        case "double":
            dat = np.array(np.fromfile(fil, np.float64, np.prod(shape)))
        case "char":
            dat = np.array(np.fromfile(fil, np.ubyte, np.prod(shape)))
            dat = ''.join([chr(it) for it in dat])
        case "complex":
            dat = np.array(np.fromfile(fil, np.double, np.prod(shape)*2))
            dat = dat[0::2] + 1.0j*dat[1::2]
    
    if np.prod(shape) == 1:
        dat = dat[0]
    elif dtype != "char":
        dat = np.reshape(dat, shape)
    
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

def getConstant(name:_CONSTANT_NAMES, fol:str):
    with open(combinePath(fol, name + ".dat"), 'rb') as fil:
        typ = readData(fil, 'int')
        dat = readData(fil, _CONSTANT_DTYPES[name])
    return dat, typ

def getPsi2t(fol:str):
    nElec,_ = getConstant("nElec", fol)
    with open(combinePath(fol, "psi2t.dat"), 'rb') as fil:
        typ = readData(fil, "int")
        nx = readData(fil, "int")
        nt = readData(fil, "int")
        dat = readData(fil, "double", (nt,nElec,nx)).swapaxes(0,1)
        xs = readData(fil, "double", nx)
        ts = readData(fil, "double", nt)
    return dat, xs, ts, typ

def getWghts(fol:str):
    with open(combinePath(fol, "wghts.dat"), 'rb') as fil:
        typ = readData(fil, "int")
        nElec = readData(fil, "int")
        wghts = readData(fil, "double", nElec)
    return wghts, typ

def getFluxSpecVD(fol:str, vdNum:int=0):
    nElec,_ = getConstant("nElec", fol)
    with open(combinePath(fol, f"{vdNum:d}" + "fluxspecvd.dat"), 'rb') as fil:
        typ = readData(fil, "int")
        readData(fil, "int") #skip VD index
        name = readData(fil, "char", 4)
        posIdx = readData(fil, "int", 1)
        nSamp = readData(fil, "int", 1)
        maxE = readData(fil, "double", 1)
        
        dftl = readData(fil, "complex", (nElec, nSamp))
        dftr = readData(fil, "complex", (nElec, nSamp))
    return dftl, dftr, maxE, posIdx, name, typ

def getExpectE0(fol:str):
    nElec,_ = getConstant("nElec", fol)
    with open(combinePath(fol, "expectE0.dat"), 'rb') as fil:
        typ = readData(fil, "int")
        e0 = readData(fil, "double", nElec)
    return e0, typ

def getSubFols(basefol:str):
    return [os.path.join(basefol, o) + "/" for o in os.listdir(basefol) if os.path.isdir(os.path.join(basefol,o))]