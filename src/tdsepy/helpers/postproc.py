from .rawdataload import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def plot1DElectronDensity(fol:str, elecNum:int = -1, vmin:float=20, vmax:float=27, cmap="magma", plot:bool = True):
    """Plots the 1-D collective electron density as a function of time for a selection of states.

    Args:
        fol (str): Folder containing data.
        elecNum (int | list , optional): Selected 1-D state. 
                                    Integer for a selected single state.
                                    List for multiple selected states.
                                    -1 for all states. Defaults to -1.
        vmin (float, optional): Log-scale minimum value. Defaults to 20.
        vmax (float, optional): Log-scale maximum value. Defaults to 27.
        cmap (str, optional): Colormap. Defaults to "magma".
        plot (bool, optional): Run plt.show(). Defaults to True.
        
    Returns:
        fig: Matplotlib figure.
        ax: Matplotlib axis.
    """
    
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
    
    if plot:
        plt.show()
    
    return plt.gcf(), plt.gca()

def plotPotential(fol:str, plot:bool=True):
    """Plots the potential as a function of time.

    Args:
        fol (str): Folder containing data.
        plot (bool, optional): Run plt.show(). Defaults to True.
    """
    dat, xs, ts, _ = getVfunct(fol)
    plt.pcolormesh(xs, ts, dat)
    
    if plot:
        plt.show()
    
    return plt.gcf(), plt.gca()
    
    
def get1DStateFluxSpectrum(fol:str, vdNum:int = 0, minE:float = 0, maxE:float = 500):
    """Gets the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k)) for each state.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        minE (float, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500.
    Returns:
        es [eV]: Signed kinetic energy, shape (nSamp).
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
   
def get1DTotalFluxSpectrum(fol:str, vdNum:int = 0, elecNum = -1, minE:float = 0, maxE:float = 500):
    """Gets the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k)) summed over all states.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        minE (float, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500.
    Returns:
        es [eV]: Signed kinetic energy, shape (nSamp).
        yld [ 1 / m^2 eV ]: Weighed bidirectional flux spectrum, shape (nSamp).
    """
    
    es, spc = get1DStateFluxSpectrum(fol, vdNum, minE, maxE)
    
    if type(elecNum) is list:
        spc = np.sum(spc[np.array(elecNum), :], axis=0)
    elif elecNum == -1:
        spc = np.sum(spc, axis=0)
    else:
        spc = spc[elecNum, :]
    
    return es, spc

def plot1DFluxSpectrum(fol:str, vdNum:int = 0, elecNum = -1, minE:float = 0, maxE:float = 500, plot:bool = True):
    """Plots the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k))

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        elecNum (int | list, optional): Selected electron states. -1 to include all, or a list to include selected states. Defaults to -1.
        minE (float, optional) [eV]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [eV]: Maximum signed kinetic energy. Defaults to 500.
    """
    es, spc = get1DTotalFluxSpectrum(fol, vdNum, elecNum, minE, maxE)
    
    plt.semilogy(es, spc)
    
    if plot:
        plt.show()
    else:
        return plt.gcf(), plt.gca()
    
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
    es, spc = get1DStateFluxSpectrum(fol, vdNum, minE, maxE)
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
        Vacuum must be potential reference level.

    Args:
        fol (str): Folder containing data.

    Returns:
        mtes [J]: Intrinsic MTE of each state, shape (nElec).
    """    
    e0s,_ = getExpectE0(fol)
    ef,_ = getConstant("ef", fol)
    wf,_ = getConstant("wf", fol)
    e0s += ef + wf
    
    mtes = 0.5*(ef - e0s)
    return mtes

