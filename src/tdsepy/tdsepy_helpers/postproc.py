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