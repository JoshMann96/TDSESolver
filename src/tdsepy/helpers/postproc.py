from .rawdataload import *
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.interpolate import interp1d
from scipy.signal import windows
from scipy import constants as cons
from scipy.fft import fft
from typing import Any

def plot1DElectronDensity(fol:str, elecNum:int = -1, vmin:float=-11, vmax:float=-4, cmap="magma", ax = None, difference=False) -> tuple[plt.Figure, plt.Axes]|None:
    """Plots the 1-D collective electron density as a function of time for a selection of states.
    AXIS | VAR | UNIT
       x |  x  | nm
       y |  t  | fs
       z |  n  | 1/a0^3

    Args:
        fol (str): Folder containing data.
        elecNum (int | list , optional): Selected 1-D state. 
                                    Integer for a selected single state.
                                    List for multiple selected states.
                                    -1 for all states. Defaults to -1.
        vmin (float, optional): Log-scale minimum value. Defaults to 20.
        vmax (float, optional): Log-scale maximum value. Defaults to 27.
        cmap (str, optional): Colormap. Defaults to "magma".
        ax (axis, optional): Axis to plot on. Defaults to None (create own fig, ax and return).
        difference (bool, optional): Plot difference between from initial state. Defaults to False.
        
    Returns:
        if plot is None:
            fig: Matplotlib figure.
            ax: Matplotlib axis.
    """
    
    dat, xs, ts, _ = getPsi2t(fol)
    wghts, _ = getWghts(fol)
    
    #retain weight only for desired states
    if elecNum != -1:
        desWght = wghts[elecNum]
        wghts *= 0
        wghts[elecNum] = desWght
    
    fig = None
    if ax is None:
        fig, ax = plt.subplots()
    
    if difference:
        im = ax.pcolormesh(xs*1e9, ts*1e15, (np.tensordot(wghts, dat, (0,0)) - np.tensordot(wghts, dat[:,0,:], (0,0)))*(cons.physical_constants["atomic unit of length"][0]**3), cmap=cmap, norm = colors.CenteredNorm())
    else:
        im = ax.pcolormesh(xs*1e9, ts*1e15, np.log10(np.tensordot(wghts, dat, (0,0))*(cons.physical_constants["atomic unit of length"][0]**3)), cmap=cmap, vmin=vmin, vmax=vmax)
    
    ax.set_xlabel(r"$x$ (nm)")
    ax.set_ylabel(r"$t$ (fs)")
    ax.set_title("Electron Density")
    cb = plt.colorbar(im, ax=ax)
    if difference:
        cb.set_label(r"$\Delta n$ (1/a$_0^3$)")
    else:
        cb.set_label(r"$\log n$ (1/a$_0^3$)")

    if fig is not None:
        return im, fig, ax
    else:
        return im

def plotPotential(fol:str, ax = None, potIndex = -1) -> Any|tuple[Any, plt.Figure, plt.Axes]:
    """Plots the potential as a function of time.
    AXIS | VAR | UNIT
       x |  x  | nm
       y |  t  | fs
       z |  U  | eV

    Args:
        fol (str): Folder containing data.
        ax (axis, optional): Axis to plot on. Defaults to None (create own fig, ax and return).
        potIndex (int, optional): Index of potential file. Defaults to -1 (no index).
    
    Returns:
        if plot is None:
            im: Matplotlib pcolormesh object.
            fig: Matplotlib figure.
            ax: Matplotlib axis.
        else:
            im: Matplotlib pcolormesh object.
    """
    dat, xs, ts, _ = getVfunct(fol, potIndex)
    
    fig = None
    if ax is None:
        fig, ax = plt.subplots()
        
    im = ax.pcolormesh(xs*1e9, ts*1e15, dat/cons.eV)
    ax.set_xlabel(r"$x$ (nm)")
    ax.set_ylabel(r"$t$ (fs)")
    ax.set_title("Potential")
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(r"$V$ (eV)")
    
    if fig is not None:
        return im, fig, ax
    else:
        return im
    
    
def get1DStateFluxSpectrum(fol:str, vdNum:int = 0, minE:float = 0, maxE:float = 500*cons.e) -> tuple[np.ndarray, np.ndarray]:
    """Gets the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k)) for each state.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        minE (float, optional) [J]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [J]: Maximum signed kinetic energy. Defaults to 500 eV.
    Returns:
        es [eV]: Signed kinetic energy, shape (nSamp).
        yld [ 1 / m^2 eV ]: Weighed bidirectional flux spectrum, shape (nElec, nSamp).
    """
    nes, momenta, psik = getFluxSpecVD(fol, vdNum)[:3]
    wghts,_ = getWghts(fol)
    
    yld = abs(psik)**2
    yld[:, nes == 0.0] = 0  # zero flux for zero momentum
    
    es = nes[(nes < maxE) & (nes > minE)]
    yld = interp1d(nes, yld, axis=-1)(es)*np.broadcast_to(wghts[:,None], (len(wghts), len(es)))

    return es/cons.e, yld
   
def get1DTotalFluxSpectrum(fol:str, vdNum:int = 0, elecNum = -1, minE:float = 0, maxE:float = 500*cons.e) -> tuple[np.ndarray, np.ndarray]:
    """Gets the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k)) summed over all states.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        elecNum (int | list, optional): Selected electron states. -1 to include all, or a list to include selected states. Defaults to -1.
        minE (float, optional) [J]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [J]: Maximum signed kinetic energy. Defaults to 500 eV.
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

def plot1DFluxSpectrum(fol:str, vdNum:int = 0, elecNum = -1, minE:float = 0, maxE:float = 500*cons.e, ax = None) -> tuple[Any, plt.Figure, plt.Axes]|Any:
    """Plots the bidirectional density flux spectrum with respect to the signed kinetic energy (sgn(E) = sgn(k))
    AXIS | VAR | UNIT
       x |  E  | eV
       y | flx | 1 / m^2 eV
       
    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        elecNum (int | list, optional): Selected electron states. -1 to include all, or a list to include selected states. Defaults to -1.
        minE (float, optional) [J]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [J]: Maximum signed kinetic energy. Defaults to 500 eV.
        ax (axis, optional): Axis to plot on. Defaults to None (create own fig, ax and return).
    
    Returns:
        if plot is None:
            im: Matplotlib semilogy object.
            fig: Matplotlib figure.
            ax: Matplotlib axis.
        else:
            im: Matplotlib semilogy object.
    """
    es, spc = get1DTotalFluxSpectrum(fol, vdNum, elecNum, minE, maxE)

    fig = None
    if ax is None:
        fig, ax = plt.subplots()
    
    im = ax.semilogy(es, spc)
    ax.set_xlabel(r"$E$ (eV)")
    ax.set_ylabel(r"$n(E)$ (1/m$^2$ eV)")
    
    if fig is not None:
        return im, fig, ax
    else:
        return im
    
def get1DStateYield(fol:str, vdNum:int = 0, minE:float = 0, maxE:float=500*cons.e) -> np.ndarray:
    """Gets yield for each 1-D state within energy range.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        minE (float, optional) [J]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [J]: Maximum signed kinetic energy. Defaults to 500 eV.

    Returns:
        yld [ 1 / m^2 ]: Weighed yield, shape (nElec).
    """    
    es, spc = get1DStateFluxSpectrum(fol, vdNum, minE, maxE)
    return np.trapezoid(x = es, y = spc, axis=1)

def get1DTotalYield(fol:str, vdNum:int = 0, elecNum = -1, minE:float = 0, maxE:float = 500*cons.e) -> float:
    """Gets total yield for selected 1-D states within energy range.

    Args:
        fol (str): Folder containing data.
        vdNum (int, optional): Virtual detector index. Defaults to 0.
        elecNum (int | list, optional): Selected electron states. -1 to include all, or a list to include selected states. Defaults to -1.
        minE (float, optional) [J]: Minimum signed kinetic energy. Defaults to 0.
        maxE (float, optional) [J]: Maximum signed kinetic energy. Defaults to 500 eV.
        
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

def getIntrinsicMTEs(fol:str) -> np.ndarray:
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

def getDipoleRadiationSpectrum(fol:str, tmin:float = None, tmax:float = None, tukeyAlpha:float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """Gets the dipole radiation spectrum.

    Args:
        fol (str): Folder containing data.
        tmin (float, optional) [s]: Minimum time. Defaults to None.
        tmax (float, optional) [s]: Maximum time. Defaults to None.

    Returns:
        es [eV]: Energy.
        spc [#/eV /nm^4]: Dipole radiation spectrum (photon count).
    """
    ts, expectA = getExpectA(fol)[0:2]
    wghts,_ = getWghts(fol)
    dt,_ = getConstant("dt", fol)

    timeMask = np.logical_and(ts >= tmin if tmin is not None else np.ones(ts.shape),
                              ts <= tmax if tmax is not None else np.ones(ts.shape))
    expectA = expectA.T[:, timeMask]
    ts = ts[timeMask]

    nt = len(ts)

    dip = cons.e * (wghts @ expectA)
    dip = dip - np.mean(dip)

    es = np.arange(0, nt//2-1) * (2*cons.pi/(ts[-1]-ts[0])) * cons.hbar/cons.eV

    spc = fft(dip * windows.tukey(nt, alpha=tukeyAlpha, sym=False))*dt
    spc = cons.mu_0/(6*cons.pi*cons.h*cons.c) * np.abs(spc)**2
    spc = spc[1:nt//2] / (es/cons.e) / cons.e * 1e-36 #convert to #/eV/nm^4
    spc[0] = 0

    return es, spc