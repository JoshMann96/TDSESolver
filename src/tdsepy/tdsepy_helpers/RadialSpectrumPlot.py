import matlab.engine
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from scipy import interpolate
from scipy.optimize import minimize
from scipy.io import savemat

eng = matlab.engine.start_matlab()
mpl.rcParams.update({'font.size':14})
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'mathtext.fontset':'stix'})

def getIntrinsicMTECont(fol):
    e0s = np.array(eng.getExpectE0(fol, nargout=1))
    ef = eng.getEf(fol, nargout=1)
    wf = eng.getW(fol, nargout=1)
    e0s += ef + wf
    
    mtes = 0.5*(ef - e0s)
    return mtes

spec_res = {}

#returns the energies (eV) and spectral intensity (electrons per eV m^2) and intrinsic MTE (eV)
def getSpectrum(fol, maxE = 600):
    txt = fol + "__maxE__" + str(maxE)
    if txt not in spec_res:
        (es, int, posIndex, name) = eng.getFluxSpectrum(fol, 0, 0, maxE, nargout=4)
        es = np.concatenate(([[0]], es), axis=1)
        int = np.concatenate((np.array(int)[:,0][:,np.newaxis], int), axis=1)
        
        mte_conts = getIntrinsicMTECont(fol)
        mte_avg = np.sum(mte_conts*int, axis=0) / np.sum(int, axis=0)
        spec_res[txt] = ( np.array(np.transpose(es)), np.array(np.transpose(np.sum(int, axis=0))), mte_avg/1.602e-19 )
    
    return spec_res[txt]

yld_res = {}

def getYield(fol, maxE = 600):
    txt = fol + "__maxE__" + str(maxE)
    if txt not in yld_res:
        yld_res[txt] = eng.getYieldByFluxSpec(fol, 0, 0, maxE)
    return yld_res[txt]

def getEmax(fol):
    return np.array( eng.getEmax(fol) )
    
def getSubFols(basefol):
    return [os.path.join(basefol, o) + "/" for o in os.listdir(basefol) if os.path.isdir(os.path.join(basefol,o))]

def getRadialSpectrum(basefol, maxE=600, nTPts=100, nRPts=200, fieldMax=80e9, fieldFunc=lambda theta: 80e9*np.cos(theta),
                        pondInterpWindowSize=10, plotSpectraCheck=False):
                        
    #array initialziation
    thetas = np.linspace(np.pi/2, 0, nTPts)
    Rs = np.linspace(0, maxE, nRPts)
    eval_ponds = np.linspace(0, maxE / fieldMax**2 * pondInterpWindowSize, nRPts*pondInterpWindowSize)
    eval_fields = np.array([fieldFunc(theta) for theta in thetas])


    #get list of folders, log-spectra, max fields
    fols = getSubFols(basefol)
    (es, int, mtes) = zip(*[getSpectrum(fol) for fol in fols])
    es = np.squeeze(np.array(es))
    mtes = np.squeeze(np.array(mtes))
    int0 = np.log10(np.squeeze(np.array(int)))
    emax = np.array([getEmax(fol) for fol in fols])
    ylds = np.array([getYield(fol) for fol in fols])

    minInt = np.amin(int0)

    #map each spectrum individually to ponderomotive units (by interpolation)
    int_pu = np.zeros((len(fols), nRPts*pondInterpWindowSize))
    mte_pu = np.zeros((len(fols), nRPts*pondInterpWindowSize))
    for i in range(len(fols)):
        if emax[i] == 0:
            int_pu[i,:] = minInt
            mte_pu[i,:] = 0.0
        else:
            f = interpolate.interp1d(es[i,:] / emax[i]**2, int0[i,:], kind='linear', fill_value=minInt, bounds_error=False)
            int_pu[i,:] = f(eval_ponds)
            f = interpolate.interp1d(es[i,:] / emax[i]**2, mtes[i,:], kind='linear', fill_value=0, bounds_error=False)
            mte_pu[i,:] = f(eval_ponds)
    del es

    #check spectra if desired
    if plotSpectraCheck:
        plt.plot(eval_ponds, int_pu.T)
        plt.legend(emax)
        plt.show()

    #sort by field
    sort_idx = np.argsort(emax)
    int_pu = np.array(int_pu[sort_idx, :])
    mte_pu = np.array(mte_pu[sort_idx, :])
    emax = np.array(emax[sort_idx])
    ylds = np.array(ylds[sort_idx])

    #interpolate between spectra in ponderomotive units
    interp = interpolate.RectBivariateSpline(emax, eval_ponds, int_pu)
    int_pu = interp(eval_fields, eval_ponds)
    #interpolate iMTE
    interp = interpolate.RectBivariateSpline(emax, eval_ponds, mte_pu)
    mte_pu = interp(eval_fields, eval_ponds)
    #interpolate yield
    interp = interpolate.interp1d(emax, ylds)
    ylds = interp(eval_fields)

    #map each spectrum back to eV by interpolation
    int = np.zeros((nTPts, nRPts))
    mte = np.zeros((nTPts, nRPts))
    for i in range(nTPts):
        f = interpolate.interp1d(eval_ponds * eval_fields[i]**2, int_pu[i,:], kind='linear', fill_value=minInt, bounds_error=False)
        int[i,:] = f(Rs)
        f = interpolate.interp1d(eval_ponds * eval_fields[i]**2, mte_pu[i,:], kind='linear', fill_value=0, bounds_error=False)
        mte[i,:] = f(Rs)
        
    #duplicate for full 180 degrees
    int = np.concatenate((int, np.flip(int[:-1,:], axis=0)), axis=0)
    mte = np.concatenate((mte, np.flip(mte[:-1,:], axis=0)), axis=0)
    thetas = np.concatenate((thetas, -np.flip(thetas[:-1])))
    ylds = np.concatenate((ylds, np.flip(ylds[:-1])))

    #meshgrid
    Rs, thetas = np.meshgrid(Rs, thetas)
    
    return Rs, thetas, int, ylds, mte

def plotRadialSpectrum(basefol, maxE=600, radialGrids=np.array([0,1,2,3,4,5,6])*100, 
                        nTPts=100, nRPts=200, cmin=5, fieldMax=80e9,
                        fieldFunc=lambda theta: 80e9*np.cos(theta),
                        pondInterpWindowSize=10,
                        plotTitle=r"Angle-Resolved Spectrum $E_{max}=80$ V/nm"):
    
    Rs, thetas, int, ylds, imte = getRadialSpectrum(basefol, maxE=maxE, nTPts=nTPts, nRPts=nRPts,
                    fieldMax=fieldMax, fieldFunc=fieldFunc,
                    pondInterpWindowSize=pondInterpWindowSize)

    #plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set(thetamin=-90, thetamax=90, theta_zero_location='N')
    im=ax.pcolormesh(thetas, Rs, int, cmap = 'inferno', vmin=cmin)
    #colorbar
    cax = fig.add_axes([0.85, 0.27, 0.03, 0.5])
    fig.colorbar(im, cax=cax, orientation='vertical')
    cax.set_title(r"log$_{10}$ $I$ ($e$ m$^{-2}$ eV$^{-1}$)")
    #10U_p cutoff
    ax.plot(thetas[:,0], 10*(1.602e-19)/(4*9.11e-31*(2*np.pi*3e8)**2)*(800e-9)**2 * (fieldFunc(thetas[:,0]))**2, 'r-', label=r"$10U_p$ Cutoff")
    #Yield
    ax.plot(thetas[:,0], ylds * maxE / np.max(ylds)*0.9, 'c-', label="Yield (arb. u.)")
    ax.set_title(plotTitle)
    ax.legend(loc = 'upper left', fancybox=False)

    #radial grids
    for rad in radialGrids:
        ax.plot(np.linspace(-np.pi/2, np.pi/2, 100), np.full(100, rad), 'w--', linewidth=0.5)
    ax.set_rgrids(radialGrids)
    plt.gcf().canvas.draw()
    labels = []
    for label in ax.get_yticklabels():
        x,y = label.get_position()
        lab = ax.text(x,y, label.get_text(), transform=label.get_transform(),
                      ha=label.get_ha(), va=label.get_va())
        lab.set_rotation(30)
        labels.append(lab)
    ax.set_yticklabels([])

    label_position=ax.get_rlabel_position()
    ax.text(-15*np.pi/24,2*ax.get_rmax()/3.,r"$E$ (eV)",
            rotation=0,ha='center',va='center')
            
    ax.set_position([0.083, 0.11, 0.9-0.083, 0.88-0.11])

    return fig, ax

#nRPts is used both for the energy spectra resolution and the resolution for energy in sweeps
#intrinsicMTE is the MTE due to transverse momentum in the material, typically ~500 meV
#sourceDim is the dimensionality of the source (tip = 0, blade = 1, surface = 2 [not implemented])
#NEED TO CALCULATE INTRINSIC MTE INTERNALLY... MAY BE MUCH MORE COMPLICATED
def plotBeamVariables_OneField(basefol, maxE=600, nTPts=100, nRPts=200, fieldMax=80e9, fieldFunc=lambda theta: 80e9*np.cos(theta),
                        pondInterpWindowSize=10, intrinsicMTE=0, sourceDim=1, mtemax=-1, pulse_len=8e-15):
                        
    Rs, thetas, int, ylds, imte = getRadialSpectrum(basefol, maxE=maxE, nTPts=nTPts, nRPts=nRPts,
                    fieldMax=fieldMax, fieldFunc=fieldFunc,
                    pondInterpWindowSize=pondInterpWindowSize)
    
    
    fig1, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    
    #Yield
    col = 'k'
    yld = np.log10(np.sum(np.power(10,int), axis=0) / int.shape[0] * np.pi)
    ax1.plot(Rs[0,:], yld, color = col, label='Yield')
    ax1.set_ylabel(r'log$_{10}$ Yield ($e$ m$^{-2}$ eV$^{-1}$)', color=col)
    ax1.tick_params(axis='y', labelcolor=col)
    
    #Curvature MTE
    col = 'tab:red'
    te = Rs * np.sin(thetas)**2
    if sourceDim == 0:
        #weighting is sin(theta)
        mte = np.sum(te*np.power(10,int)*np.sin(thetas), axis=0)/np.sum(np.power(10,int)*np.sin(thetas), axis=0)
    elif sourceDim == 1:
        #weighting is constant
        mte = np.sum(te*np.power(10,int), axis=0)/np.sum(np.power(10,int), axis=0)
    else:
        #no geometric MTE
        mte = np.full(Rs[0,:].shape, intrinsicMTE)
    ax2.plot(Rs[0,:], mte, color=col, label='Curvature MTE')
    print("Maximum Curvature MTE: ", np.max(mte))
    ax1.set_xlabel(r'$E$ (eV)')
    ax2.set_ylabel(r'MTE (eV)', color=col)
    ax2.tick_params(axis='y', labelcolor=col)
    if mtemax != -1:
        ax2.set_ylim([0, mtemax]) 
    ax1.set_xlim([0, maxE])
    
    #Intrinsic MTE
    col = 'tab:purple'
    if sourceDim == 0:
        #weighting is sin(theta)
        imte = np.sum(imte*np.power(10, int)*np.abs(np.sin(thetas)), axis=0)/np.sum(np.power(10, int)*np.abs(np.sin(thetas)), axis=0)
    else:
        #weighting is constant
        imte = np.sum(imte*np.power(10, int), axis=0)/np.sum(np.power(10, int), axis=0)
    ax2.plot(Rs[0,:], imte, color=col, label='Intrinsic MTE')
    
    #Transverse beam size
    col = 'tab:blue'
    if sourceDim == 0:
        sigx = np.sqrt(np.sum(np.abs(np.sin(thetas)**3)*np.power(10, int), axis=0)/np.sum(np.abs(np.sin(thetas))*np.power(10,int), axis=0))
    else:
        sigx = np.sqrt(np.sum(np.sin(thetas)**2*np.power(10, int), axis=0)/np.sum(np.power(10,int), axis=0))
    ax2.plot(Rs[0,:], sigx*mtemax/np.max(sigx), color=col, label='RMS Size (arb. u.)')
    
    fig1.legend(loc = (0.4,0.67), fancybox=False)
    fig1.tight_layout()
    
    #Find local maxima for brightness
    startPts = np.array([maxE/3, maxE*9/10])*fieldMax**2/80e9**2
    if sourceDim == 0:
        totEmit = sigx**2*np.sqrt(mte**2+imte**2)/5.11e5
    else:
        totEmit = sigx*np.sqrt(np.sqrt(mte**2+imte**2)*imte)/5.11e5
        
    emitx = sigx * np.sqrt(np.sqrt(mte**2+imte**2)/5.11e5) #across blade
    emity = np.sqrt(imte/5.11e5) #along blade
    
    def rms_ener_gaus(inp):
        cent = inp[0]
        width = np.abs(inp[1])
        gaus = np.exp(-(Rs[0,:]-cent)**2/(2*width**2))
        mn_E = np.sum(np.power(10, yld) * Rs[0,:] * gaus)/np.sum(np.power(10, yld) * gaus)
        rms_E = np.sqrt(np.sum(np.power(10, yld) * (Rs[0,:]-mn_E)**2 * gaus)/np.sum(np.power(10, yld) * gaus))
        return mn_E, rms_E
    
    def rms_ener_full():
        mn_E = np.sum(np.power(10, yld) * Rs[0,:])/np.sum(np.power(10, yld))
        rms_E = np.sqrt(np.sum(np.power(10, yld) * (Rs[0,:]-mn_E)**2)/np.sum(np.power(10, yld)))
        return mn_E, rms_E
    
    def neg_brightness_gaus(inp):
        cent = inp[0]
        width = np.abs(inp[1])
        gaus = np.exp(-(Rs[0,:]-cent)**2/(2*width**2))
        (mn_E, rms_E) = rms_ener_gaus(inp)
        return -1.602e-19/pulse_len*np.sum(np.power(10, yld) * gaus)**2 / np.sum(np.power(10, yld) * totEmit * gaus) * (Rs[0,-1] - Rs[0,0]) / len(Rs[0,:]) / np.sqrt(rms_E/5.11e5)
        
    def neg_brightness_5_gaus(inp):
        cent = inp[0]
        width = np.abs(inp[1])
        gaus = np.exp(-(Rs[0,:]-cent)**2/(2*width**2))
        (mn_E, rms_E) = rms_ener_gaus(inp)
        return -1.602e-19/pulse_len*np.sum(np.power(10, yld) * gaus)**2 / np.sum(np.power(10, yld) * totEmit * gaus) * (Rs[0,-1] - Rs[0,0]) / len(Rs[0,:])
        
    def emittance_gaus(inp):
        cent = inp[0]
        width = np.abs(inp[1])
        gaus = np.exp(-(Rs[0,:]-cent)**2/(2*width**2))
        #full, x, y
        return np.sum(np.power(10, yld) * totEmit * gaus) / np.sum(np.power(10, yld) * gaus), np.sum(np.power(10, yld) * emitx * gaus) / np.sum(np.power(10, yld) * gaus), np.sum(np.power(10, yld) * emity * gaus) / np.sum(np.power(10, yld) * gaus)
    
    def neg_brightness_full():
        (mn_E, rms_E) = rms_ener_full()
        return -1.602e-19/pulse_len*np.sum(np.power(10, yld))**2 / np.sum(np.power(10, yld) * totEmit) * (Rs[0,-1] - Rs[0,0]) / len(Rs[0,:]) / np.sqrt(rms_E/5.11e5)
    
    def neg_brightness_5_full():
        (mn_E, rms_E) = rms_ener_full()
        return -1.602e-19/pulse_len*np.sum(np.power(10, yld))**2 / np.sum(np.power(10, yld) * totEmit) * (Rs[0,-1] - Rs[0,0]) / len(Rs[0,:])
    
    def emittance_full():
        #full, x, y
        return np.sum(np.power(10, yld) * totEmit)/np.sum(np.power(10, yld)), np.sum(np.power(10, yld) * emitx)/np.sum(np.power(10, yld)), np.sum(np.power(10, yld) * emity)/np.sum(np.power(10, yld))
    
    print("B6 calculated as pulse-averaged current / (x-y emittance * sqrt(RMS energy spread / average normal energy))")
    
    B6 = np.zeros(len(startPts)+1)
    B5 = np.zeros(len(startPts)+1)
    E4 = np.zeros(len(startPts)+1)
    Ex = np.zeros(len(startPts)+1)
    Ey = np.zeros(len(startPts)+1)
    
    
    for i, strt in enumerate(startPts):
        res = minimize(neg_brightness_gaus, [strt, maxE/20])
        #assumes 1-D structure for units
        print("Maximum brightness of \t\tB6 = " + str(-neg_brightness_gaus(res.x)*1e-18) + " A/um^2 mrad^2")
        print("\t\t\t\tB5 = " + str(-neg_brightness_5_gaus(res.x)*1e-18) + " A/um^2 mrad^2")
        emits = emittance_gaus(res.x)
        print("\t with an emittance of \texy = " + str(emits[0]*1e6) + " LR mrad^2")
        print("\t\t\t\tex = " + str(emits[1]*1e3) + " R mrad")
        print("\t\t\t\tey = " + str(emits[2]*1e3) + " L mrad")
        print("\t with a Gaussian at \tE = " + str(res.x[0]) + " +/- " + str(res.x[1]))
        
        B6[i] = -neg_brightness_gaus(res.x)
        B5[i] = -neg_brightness_5_gaus(res.x)
        E4[i] = emits[0]
        Ex[i] = emits[1]
        Ey[i] = emits[2]
        
        ax2.errorbar(res.x[0], mtemax/2, None, res.x[1], 'k.')
    emits = emittance_full()
    print("Full brightness of \t\tB6 = " + str(-neg_brightness_full()*1e-18) + " A/um^2 mrad^2")
    print("\t\t\t\tB5 = " + str(-neg_brightness_5_full()*1e-18) + " A/um^2 mrad^2")
    print("\t with an emittance of \texy = " + str(emits[0]*1e6) + " LR mrad^2")
    print("\t\t\t\tex = " + str(emits[1]*1e3) + " R mrad")
    print("\t\t\t\tey = " + str(emits[2]*1e3) + " L mrad")
    print("")
    
    B6[-1] = -neg_brightness_full()
    B5[-1] = -neg_brightness_5_full()
    E4[-1] = emits[0]
    Ex[-1] = emits[1]
    Ey[-1] = emits[2]
        
        
    return fig1, ax1, ax2, B6, B5, E4, Ex, Ey

#Plot radial spectrum
#parameters
basefol = "E:/Joshua/Desktop/School Stuff (Unsynced)/PBPL HHG/data/bigsim_ic/lam800/"
maxE = 100 #eV
radialGrids = np.array([0,1,2,3,4,5,6])*100
nTPts = 100
nRPts = 200
cmin = 5
#fieldMax V/m:
#fieldMax = 80e9
fieldMax = 31.74e9 #80 eV
#fieldMax = 79.3e9 #500 eV
#fieldMax = 25.1e9 #50 eV
#fieldMax = 7.94e9 #5 eV
fieldFunc = lambda theta: fieldMax*np.cos(theta)
pondInterpWindowSize = 10 #extra space for ponderomotive interpolation
#plotTitle = r"Angle-Resolved Spectrum $E_{max}=80$ V/nm"
plotTitle = r"Image Charge"

# plt.rcParams["figure.figsize"]=[8.5,6]

# (fig, ax) = plotRadialSpectrum(basefol, maxE=maxE, radialGrids=radialGrids, nTPts=nTPts, nRPts=nRPts,
                    # cmin=cmin, fieldMax=fieldMax, fieldFunc=fieldFunc,
                    # pondInterpWindowSize=pondInterpWindowSize, plotTitle=plotTitle)

# plt.subplots_adjust(left=0.08, bottom=0.08, right=0.886, top=0.84)
# #plt.show()
# plt.savefig(fname="radspec_ic.png", dpi=600)


#Sweep through fields, save result
# fieldMaxs = np.linspace(1, 80, 100)*1e9
# B6s = np.zeros((len(fieldMaxs), 3))
# B5s = np.zeros((len(fieldMaxs), 3))
# E4s = np.zeros((len(fieldMaxs), 3))
# Exs = np.zeros((len(fieldMaxs), 3))
# Eys = np.zeros((len(fieldMaxs), 3))

# for i in range(len(fieldMaxs)):
    # print(str(i/len(fieldMaxs)*100)+"%")
    # fieldFunc = lambda theta: fieldMaxs[i]*np.cos(theta)
    # (fig1, ax1, ax2, B6, B5, E4, Ex, Ey) = plotBeamVariables_OneField(basefol, maxE=maxE, nTPts=nTPts, nRPts=nRPts, fieldMax=fieldMaxs[i], fieldFunc=fieldFunc, pondInterpWindowSize=pondInterpWindowSize, mtemax=8)
    # B6s[i,:] = B6
    # B5s[i,:] = B5
    # E4s[i,:] = E4
    # Exs[i,:] = Ex
    # Eys[i,:] = Ey
    # plt.close('all')

# mdic = {"B6": B6s, "B5": B5s, "E4": E4s, "Ex": Exs, "Ey": Eys, "fields": fieldMaxs}
# savemat("brightness_emittance_nc.mat", mdic)


(fig1, ax1, ax2, B6, B5, E4, Ex, Ey) = plotBeamVariables_OneField(basefol, maxE=maxE, nTPts=nTPts, nRPts=nRPts, fieldMax=fieldMax, fieldFunc=fieldFunc, pondInterpWindowSize=pondInterpWindowSize, mtemax=-1)
print(B6, B5, E4, Ex, Ey)
#max MTE at ~50 eV
#plt.show()
plt.savefig(fname="yldmte_ic_b_atc.svg")