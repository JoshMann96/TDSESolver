import json
from .rawdataload import *
import matplotlib.pyplot as plt

def plotElectronDensity(fol, elecNum = -1, vmin=20, vmax=27, cmap="magma"):
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
    
def getFluxSpectrum(fol:str, vdNum = 0, elecNum = -1, minE = 0, maxE:float = 500):
    """Gets the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k))

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        elecNum (int | list, optional): Selected electron states. -1 to include all, or a list to include selected states. Defaults to -1.
        minE (int, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500 eV.
    Returns:
        es [eV]: Signed kinetic energy.
        yld [ 1 / m^2 eV ]: Bidirectional flux spectrum.
    """
    dftl, dftr, maxE_ = getFluxSpecVD(fol, vdNum)[0:3]
    dx,_ = getConstant("dx", fol)
    dt,_ = getConstant("dt", fol)
    wghts,_ = getWghts(fol)
    
    nSamp = dftl.shape[1]
    
    ces = maxE_/nSamp * np.arange(1, nSamp)
    ces = ces[ces <= max(abs(minE), abs(maxE))]
    ks = np.sqrt(2*9.11e-31*ces*1.602e-19/1.11212e-68)
    ces = ces[ks < np.pi/dx]
    ks = ks[ks < np.pi/dx]
    
    dftl = dftl[:,1:nSamp]*dt
    dftr = dftr[:,1:nSamp]*dt
    
    phi = np.exp(0.5j*dx*ks)
    phist = np.conj(phi)
    muldiv = 1.0/(-2.0j*np.sin(dx*ks))
    
    r = (phist*dftl - phi*dftr  )*muldiv
    l = (-phi*dftl  + phist*dftr)*muldiv
    
    es = np.concatenate((np.flip(ces), np.array([0]), ces))
    
    return es, yld
    
def plotFluxSpectrum(fol:str, vdNum = 0, elecNum = -1, minE = 0, maxE:float = 500):
    pass