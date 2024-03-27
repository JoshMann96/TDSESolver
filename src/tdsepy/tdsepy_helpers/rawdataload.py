import json
import numpy as np
from typing import Literal

_INT_SIZE = np.dtype(np.int32).itemsize
_DOUBLE_SIZE = np.dtype(np.float64).itemsize
_CONSTANT_NAMES = Literal["dx", "dt", "emax", "lam", "tau", "rad", "ef", "wf", "nElec", "nPts", "nSteps"]
_CONSTANT_DTYPES = {
    "dx" : np.float64,
    "dt" : np.float64,
    "emax" : np.float64,
    "lam" : np.float64,
    "tau" : np.float64,
    "rad" : np.float64,
    "ef" : np.float64,
    "wf" : np.float64,
    "nElec" : np.int32,
    "nPts" : np.int32,
    "nSteps" : np.int32
}

def getConstant(name:_CONSTANT_NAMES, fol:str):
    if fol[-1] != '/':
        fol += '/'
    with open(fol+name+".dat", 'rb') as fil:
        typ, = np.frombuffer(fil.read(_INT_SIZE), np.int32)
        emax, = np.frombuffer(fil.read(np.dtype(_CONSTANT_DTYPES[name]).itemsize),_CONSTANT_DTYPES[name])
    return emax, typ