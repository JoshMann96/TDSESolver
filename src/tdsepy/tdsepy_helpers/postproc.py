import json
from .rawdataload import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def plot1DElectronDensity(fol:str, elecNum:int = -1, vmin:float=20, vmax:float=27, cmap="magma"):
    dat, xs, ts, _ = getPsi2t(fol)
    wghts, _ = getWghts(fol)
    
    print(dat.shape)
    print(wghts.shape)
    
    #retain weight only for desired states
    if elecNum != -1:
        desWght = wghts[elecNum]
        wghts *= 0
        wghts[elecNum] = desWght
    
    s = plt.pcolormesh(xs, ts, np.log10(np.tensordot(wghts, dat, (0,0))), cmap=cmap, vmin=vmin, vmax=vmax)
    plt.show()
    
def get1DFluxSpectrum(fol:str, vdNum:int = 0, minE:float = 0, maxE:float = 500):
    """Gets the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k))

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        minE (float, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500.
    Returns:
        es [eV]: Signed kinetic energy.
        yld [ 1 / m^2 eV ]: Weighed bidirectional flux spectrum, shape (nElec, nSamp).
    """
    dftl, dftr, maxE_ = getFluxSpecVD(fol, vdNum)[0:3]
    maxE_ = maxE_ / 1.602e-19
    dx,_ = getConstant("dx", fol)
    dt,_ = getConstant("dt", fol)
    wghts,_ = getWghts(fol)
    
    nSamp = dftl.shape[1]
    
    ces = maxE_/nSamp * np.arange(1, nSamp)
    ces = ces[ces <= max(abs(minE), abs(maxE))]
    ks = np.sqrt(2*9.11e-31*ces*1.602e-19/1.11212e-68)
    ces = ces[ks < np.pi/dx]
    ks = ks[ks < np.pi/dx]
    
    dftl = dftl[:,1:len(ces)+1]*dt
    dftr = dftr[:,1:len(ces)+1]*dt
    
    phi = np.exp(0.5j*dx*ks)
    phist = np.conj(phi)
    muldiv = 1.0/(-2.0j*np.sin(dx*ks))
    
    r = (phist*dftl - phi*dftr  )*muldiv
    l = (-phi*dftl  + phist*dftr)*muldiv
    
    nes = np.concatenate((np.flip(-ces), np.array([0]), ces))
    ces_exp = np.broadcast_to(ces[None,:]**(0.25), (len(wghts), len(ces)))
    yld = (1.602e-19)**(3.0/2) / (np.pi*1.0546e-34*np.sqrt(2*9.11e-31)) * np.abs(np.concatenate((\
        np.flip(l*ces_exp  , axis=1), np.zeros((len(wghts),1)), r*ces_exp),\
        axis=1))
    
    es = nes[(nes < maxE) & (nes > minE)]
    yld = interp1d(nes, yld, axis=-1)(es)*np.broadcast_to(wghts[:,None], (len(wghts), len(es)))
    
    return es, yld
    
def plot1DFluxSpectrum(fol:str, vdNum:int = 0, elecNum = -1, minE:float = 0, maxE:float = 500):
    """Plots the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k))

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        elecNum (int | list, optional): Selected electron states. -1 to include all, or a list to include selected states. Defaults to -1.
        minE (float, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500.
    """
    es, yld = get1DFluxSpectrum(fol, vdNum, minE, maxE)
    
    if type(elecNum) is list:
        yld = np.sum(yld[np.array(elecNum), :], axis=0)
    elif elecNum == -1:
        yld = np.sum(yld, axis=0)
    else:
        yld = yld[elecNum, :]
    
    plt.semilogy(es, yld)
    
    plt.show()
    
def get1DStateYield(fol:str, vdNum:int = 0, minE:float = 0, maxE:float=500):
    """Gets yield for each 1-D state within energy range.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        minE (float, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500.

    Returns:
        yld [ 1 / m^2 ]: Weighed yield, shape (nElec).
    """    
    es, spc = get1DFluxSpectrum(fol, vdNum, minE, maxE)
    return np.trapz(x = es, y = spc, axis=1)

def get1DTotalYield(fol:str, vdNum:int = 0, elecNum = -1, minE:float = 0, maxE:float = 500):
    """Gets total yield for selected 1-D states within energy range.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        elecNum (int | list, optional): Selected electron states. -1 to include all, or a list to include selected states. Defaults to -1.
        minE (float, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500.
        
    Returns:
        yld [ 1 / m^2 ]: Weighed yield.
    """    
    ylds = get1DStateYield(fol, vdNum, minE, maxE)
    
    if type(elecNum) is list:
        yld = np.sum(ylds[elecNum])
    elif elecNum == -1:
        yld = np.sum(ylds)
    else:
        yld = ylds[elecNum]
        
    return yld

def getIntrinsicMTEs(fol:str):
    """Gets intrinsic mean transverse energy (MTE) due to transverse crystal momentum in material, assuming the free electron gas model.

    Args:
        fol (str): Folder containing data.

    Returns:
        mtes [J]: Intrinsic MTE of each state, shape (nElec).
    """    
    e0s = getExpectE0(fol)
    ef = getConstant("ef", fol)
    wf = getConstant("wf", fol)
    e0s += ef + wf
    
    mtes = 0.5*(ef - e0s)
    return mtes

