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
import matplotlib.pyplot as pls

def getSubFols(basefol:str):
    return [os.path.join(basefol, o) + "/" for o in os.listdir(basefol) if os.path.isdir(os.path.join(basefol,o))]

def getRadialSpectrum(superfol:str, vdNum:int=0, minE:float=0, maxE:float=600, outputMaxE:float=600, 
                      fieldMax:float=80e9, fieldProfile=lambda theta:np.cos(theta), nTheta=100, nRadial=200):
    
    if outputMaxE > maxE:
        raise ValueError("outputMaxE must be less than or equal to maxE")
    
    #get list of folders, log-spectra, peak fields
    fols = getSubFols(superfol)
    nFols = len(fols)
    
    #get total flux spectrum, summed over states, take log10
    (ess, spcs) = zip(*[((res := get1DStateFluxSpectrum(fol, vdNum, minE, maxE))[0], np.log10(np.sum(res[1], axis=0).squeeze())) for fol in fols])
    minSpc = min([min(spc) for spc in spcs])
    eMaxs = [getConstant("emax", fol) for fol in fols]
    
    #other stuff probably for MTE calculation, probably for other function
    #stateYlds = [get1DStateYield(fol, vdNum, minE, maxE) for fol in fols]
    #imtes = [yld * getIntrinsicMTEs(fol) / np.sum(yld)  for (yld, fol) in zip(ylds, fols)]
    
    #map each spectrum individually to ponderomotive units
    #makes interpolation across theta/field strength more smooth
    #number of radial points is max of output shape
    nR = max([len(spc) for spc in spcs])
    spc_pu = np.zeros((nFols, nR))
    eval_ponds = np.linspace(0, maxE / fieldMax**2, nR)
    for (i,(eMax, es, spc)) in enumerate(zip(eMax, ess, spcs)):
        if eMax == 0:
            spc_pu[i,:] = minSpc
        else:
            f = interpolate.interp1d(es / eMax**2, spc, kind='linear', fill_value=minSpc, bounds_error=False)
            spc_pu[i,:] = f(eval_ponds)
    
    #sort by field
    sort_idx = np.argsort(eMax)
    spc_pu = np.array(spc_pu[sort_idx,:])
    eMax = np.array(eMax[sort_idx])
    
    #interpolate between spectra
    thetas = np.linspace(np.pi/2, 0, nTheta)
    eval_fields = fieldMax*np.array([fieldProfile(theta) for theta in thetas])
    f = interpolate.RectBivariateSpline(eMax, eval_ponds, spc_pu)
    spc_pu = f(eval_fields, eval_ponds)
    
    
    ##test
    plt.pcolor(eval_fields, eval_ponds, spc_pu)
    plt.show()
    ## NEED TO TEST HERE, SEE IF PLOTS LOOK REASONABLE
    
    #map spectra back to eV on uniform grid
    es = np.linspace(0.0, outputMaxE, nRadial)
    spc = np.zeros((nTheta, nRadial))
    for i in range(nTheta):
        f = interpolate.interp1d(eval_ponds * eval_fields[i]**2, spc_pu[i,:], kind='linear', fill_value=minSpc, bounds_error=False)
        spc[i,:] = f(es)
    
    return eval_fields, es, spc
    
    