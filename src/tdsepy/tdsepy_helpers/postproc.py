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
    
def getFluxSpectrum(fol, vdNum = 0, elecNum = -1, minE = 0, maxE = 500*1.602e-19):
    dftl, dftr, maxE_, _ = getFluxSpecVD(fol, vdNum)
    dx = getConstant("dx", fol)
    dt = getConstant("dt", fol)
    wghts = getWghts(fol)
    
    nSamp = len(dftl)/len(wghts)
    
    ces = maxE_/nSamp * np.arange(1, nSamp)
    ces = ces[ces <= np.max(abs(minE), abs(maxE))]
    ks = np.sqrt(2*9.11e-31*ces*1.602e-19/1.11212e-68)
    ces = ces[ks < np.pi/dx]
    ks = ks[ks < np.pi/dx]
    
    dftl = dftl[:,2:len(ces)+1]*dt
    dftr = dftr[:,2:len(ces)+1]*dt
    
    phi = np.exp(0.5j*dx*ks)
    phist = np.conj(phi)
    muldiv = 1.0/(-2.0j*np.sin(dx*ks))
    
    r = (phist*dftl - phi*dftr  )*muldiv
    l = (-phi*dftl  + phist*dftr)*muldiv
    
    nes = np.concatenate(np.flip(ces), 0, ces)