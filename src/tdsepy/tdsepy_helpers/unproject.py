""" Maps 1-D TDSE outputs to 2-D results (i.e., for nanoblades).
    Essentially the inverse of a projection.
    Data should be structured as:
    superfol/
        data1/
        data2/
        ...
    With each data subfolder containing the results for a single field.
"""

from .postproc import *
from .rawdataload import *
import os
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.constants as cons

def getSubFols(basefol:str):
    return [os.path.join(basefol, o) + "/" for o in os.listdir(basefol) if os.path.isdir(os.path.join(basefol,o))]

def loadRadialSpectrumData(superfol:str, vdNum:int=0, minE:float=0, maxE:float=1200, calcYldMTE:bool=False):
    fols = getSubFols(superfol)
    nFols = len(fols)
    
    #get total flux spectrum, summed over states, take log10
    ess, spcs = zip(*[((res := get1DTotalFluxSpectrum(fol, vdNum, -1, minE, maxE))[0], np.log10(res[1])) for fol in fols])
    eMaxs = [getConstant("emax", fol)[0] for fol in fols]
    
    if calcYldMTE:
        stateYlds = [get1DStateYield(fol, vdNum, minE, maxE) for fol in fols]
        imtes = [np.sum(yld * getIntrinsicMTEs(fol)/cons.e / np.sum(yld))  for (yld, fol) in zip(stateYlds, fols)]
    
    return ess, spcs, eMaxs, imtes, nFols, maxE
    

def getRadialSpectrum(radialSpectrumData, outputMaxE:float=600, 
                      fieldMax:float=80e9, fieldProfile=lambda theta:np.cos(theta), nTheta=400, nRadial=400,
                      calcYldMTE:bool=False, structureDim:int=1):
    
    ess, spcs, eMaxs, imtes, nFols, maxE = radialSpectrumData
    
    if outputMaxE > maxE:
        raise ValueError("outputMaxE must be less than or equal to maxE")
    
    minSpc = np.min(spcs)
    
    #map each spectrum individually to ponderomotive units
    #makes interpolation across theta/field strength more smooth
    #number of radial points is max of output shape
    nR = max([len(spc) for spc in spcs])
    spc_pu = np.zeros((nFols, nR))
    eval_ponds = np.linspace(0, maxE / fieldMax**2, nR)
    for i, (eMax, es, spc) in enumerate(zip(eMaxs, ess, spcs)):
        if eMax == 0:
            spc_pu[i,:] = minSpc
        else:
            f = interpolate.interp1d(es / eMax**2, spc, kind='linear', fill_value=minSpc, bounds_error=False)
            spc_pu[i,:] = f(eval_ponds)
    
    #sort by field
    sort_idx = np.argsort(eMaxs)
    spc_pu = np.array(spc_pu[sort_idx,:])
    eMaxs = np.array(eMaxs)[sort_idx]
    if calcYldMTE:
        imtes = np.array(imtes)[sort_idx]
    
    #interpolate between spectra
    thetas = np.linspace(np.pi/2, 0, nTheta)
    eval_fields = fieldMax*np.array([fieldProfile(theta) for theta in thetas])
    f = interpolate.RectBivariateSpline(eMaxs, eval_ponds, spc_pu)
    spc_pu = f(eval_fields, eval_ponds)
    
    #map spectra back to eV on uniform grid
    es = np.linspace(0.0, outputMaxE, nRadial)
    spc = np.zeros((nTheta, nRadial))
    for i in range(nTheta):
        f = interpolate.interp1d(eval_ponds * eval_fields[i]**2, spc_pu[i,:], kind='linear', fill_value=minSpc, bounds_error=False)
        spc[i,:] = f(es)
        
    #calculate brightness, yield, mte
    if calcYldMTE:
        ylds = np.trapz(10**spc, es, axis=1)
        te_curv = np.trapz(10**spc * es, es , axis=1) * np.sin(thetas)**2 / ylds
        
        f = interpolate.interp1d(eMaxs, imtes, kind='linear', fill_value=min(imtes), bounds_error=False)
        te_int = f(eval_fields)
        
        #tip (Jacobian factor is sin(theta) for tip)
        if structureDim == 0:
            yld = -np.trapz(ylds * np.sin(thetas), thetas, axis=0)
            mte_curv = -np.trapz(te_curv * np.sin(thetas) * ylds, thetas, axis=0) / yld
            mte_int = -np.trapz(te_int * np.sin(thetas) * (0.5+0.5*np.cos(thetas)**2) * ylds, thetas, axis=0) / yld #cos is from half of iMTE becoming longitudinal
            rms_x = np.sqrt(-np.trapz(ylds * np.sin(thetas) * thetas**2, thetas, axis=0)/yld)
        #blade (Jacobain factor is 1 for blade)
        elif structureDim == 1:
            yld = -np.trapz(ylds, thetas, axis=0)
            mte_curv = -np.trapz(te_curv * ylds, thetas, axis=0) / yld
            mte_int = -np.trapz(te_int * ylds * (0.5+0.5*np.cos(thetas)**2), thetas, axis=0) / yld
            rms_x = np.sqrt(-np.trapz(ylds * thetas**2, thetas, axis=0)/yld)
        else:
            raise ValueError("structureDim must be 0 (tip) or 1 (blade)")
        
    if calcYldMTE:
        return thetas, eval_fields, es, spc, yld, mte_curv, mte_int, rms_x
    else:
        return thetas, eval_fields, es, spc
    
def plotRadialSpectrum(superfol:str, vdNum:int=0, minE:float=0, maxE:float=1200, outputMaxE:float=600, 
                      fieldMax:float=80e9, fieldProfile=lambda theta:np.cos(theta), nTheta=400, nRadial=400,
                      vmin:float=20, vmax:float=27):
    
    thetas, _, es, spc = getRadialSpectrum(superfol=superfol, vdNum=vdNum, minE=minE, maxE=maxE, outputMaxE=outputMaxE, 
                      fieldMax=fieldMax, fieldProfile=fieldProfile, nTheta=nTheta, nRadial=nRadial)
    
    #duplicate for full 180 degrees
    spc = np.concatenate((spc, np.flip(spc[:-1,:], axis=0)), axis=0)
    thetas = np.concatenate((thetas, -np.flip(thetas[:-1])))
    
    #plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set(thetamin=-90, thetamax=90, theta_zero_location='N')
    im=ax.pcolormesh(thetas, es, spc.T, cmap = 'inferno', vmin=vmin, vmax=vmax)
    #colorbar
    cax = fig.add_axes([0.85, 0.27, 0.03, 0.5])
    fig.colorbar(im, cax=cax, orientation='vertical')
    cax.set_title(r"log$_{10}$ $I$ ($e$ m$^{-2}$ eV$^{-1}$)")
    
    plt.show()

def calcBrightness(radialSpectrumData, outputMaxE:float=600, 
                      fieldMax:float=80e9, fieldProfile=lambda theta:np.cos(theta), nTheta=400, nRadial=400,
                      structureDim:int=1, structureRadius:float=20e-9, illuminationLength:float=1e-6):
    
    yld, mte_curv, mte_int, rms_x = getRadialSpectrum(radialSpectrumData=radialSpectrumData, outputMaxE=outputMaxE, 
                      fieldMax=fieldMax, fieldProfile=fieldProfile, nTheta=nTheta, nRadial=nRadial, calcYldMTE=True, structureDim=structureDim)[-4:]
    
    if structureDim == 0:
        eps_x = rms_x * structureRadius * np.sqrt((mte_int + mte_curv)/(cons.m_e*cons.c**2/cons.e))
        eps_y = eps_x
        b4 = yld * cons.e * structureRadius**2 / eps_x / eps_y
    elif structureDim == 1:
        eps_x = rms_x * structureRadius * np.sqrt((mte_int + mte_curv)/(cons.m_e*cons.c**2/cons.e))
        eps_y = illuminationLength*np.sqrt(mte_int/(cons.m_e*cons.c**2/cons.e))
        b4 = yld * cons.e * structureRadius*illuminationLength / eps_x / eps_y
    else:
        raise ValueError("structureDim must be 0 (tip) or 1 (blade)")
    
    return b4, yld, eps_x, eps_y, mte_curv, mte_int