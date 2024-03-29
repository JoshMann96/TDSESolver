import sys, os
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)

sys.path.insert(0, "/home/jmann/TDSESolveLinux/build/lib")
from tdsepy import *
import numpy as np
import scipy.constants as cons
import json

def pond_U(emax, lam):
    return cons.e**2 * emax**2 * lam**2 / (4* cons.m_e * (2*cons.pi*cons.c)**2)

def pond_a(emax, lam):
    return cons.e * emax * lam**2 / (cons.m_e * (2*cons.pi*cons.c)**2)

def runSimSweepFields(emaxs:list, lam:float=800e-9, rad:float=20e-9, ef:float=5.51*1.602e-19, wf:float=5.1*1.602e-19, tau:float=8e-15, data_fol:str="data/", callback=None, 
                  target_total_truncation_error:float = 0.01, min_emitted_energy:float=1.602e-19, target_elec_num:float = 50, abs_width:float=20e-9, min_timesteps:int=2000, measure_density:bool=True):
    """Runs a series of rescattering simulations within a range of peak field strengths.
        Current quantities being output:
        nPts, nSteps, dx, dt, <a>, nElec, VDFluxSpec, Weights
        Optional: Psi2t, Vfunct (via parameter measure_density)
        
        Parameters 'emax' through 'tau' are also saved as constants for each simulation.
        
        Parameters 'emaxs' through 'tau' are also saved as a dict in data_fol/params.json.
        
        All parameters are in SI units.

        PARAMETERS
        ----------
        emaxs : list
            Peak electric fields to sweep over.
        lam : float 
            Laser wavelength. Defaults to 800e-9.
        rad : float
            Apex radius of curvature (for collective image fields). Defaults to 20e-9.
        ef : float
            Fermi energy. Defaults to 5.51*1.602e-19.
        wf : float
            Work function. Defaults to 5.1*1.602e-19.
        tau : float
             Full-width half-max power. Defaults to 8e-15.
        data_fol : str 
            Folder to store output data in.
            May include final '/', or not.
            Defaults to "data/".
        callback 
            Callback function which must take an integer between 0 and 100, inclusive.
            callback is called during time-stepping with the current percentage complete, as an integer. 
            Defaults to None.
        target_total_truncation_error : float
            Upper bound for total relative truncation error of the wavefunction. Defaults to 0.01.
        min_emitted_energy : float
            Minumum emittable energy. 
            Sets temporal length of simulation such that a particle emitted with this energy at the peak of the laser pulse reaches the rightmost boundary before timestepping stops. 
            Defaults to 1.602e-19.
        target_elec_num : float
            Target number of states. 
            Sets the width of the Jellium slab such that there are about this many 1-D states below the Fermi level. Actual number of states may deviate.
            Defaults to 50.
        abs_width : float
            Width of the absorptive boundary. Defaults to 20e-9.
        min_timesteps : int
            Minimum number of timesteps. 
            Overrides target_total_truncation_error if the number of timesteps would be too few. 
            Defaults to 2000.
        measure_density : bool
            Whether to use the Psi2t and Vfunct measurers. Defaults to True.
    """    
    os.makedirs(data_fol, exist_ok=True)
    if data_fol[-1] != '/':
        data_fol += '/'
    
    emaxs = list(emaxs)
    
    paramDict = locals()
    del paramDict["callback"]
    
    with open(f"{data_fol}params.json", 'w') as fil:
        json.dump(paramDict, fil)
    for (i, emax) in enumerate(emaxs):
        runSingleSimulation(emax, lam, rad, ef, wf, tau, f"{data_fol}{i}/", callback, target_total_truncation_error, min_emitted_energy, target_elec_num, abs_width, min_timesteps, measure_density)
    
    
def runSingleSimulation(emax:float=20e9, lam:float=800e-9, rad:float=20e-9, ef:float=5.51*1.602e-19, wf:float=5.1*1.602e-19, tau:float=8e-15, data_fol:str="data/", callback=None, 
                  target_total_truncation_error:float = 0.01, min_emitted_energy:float=1.602e-19, target_elec_num:float = 50, abs_width:float=20e-9, min_timesteps:int=2000, measure_density:bool=True):
    """Runs a rescattering simulation while recording various quantites. 
        Current quantities being output:
        nPts, nSteps, dx, dt, <a>, nElec, VDFluxSpec, Weights
        Optional: Psi2t, Vfunct (via parameter measure_density)
        
        Parameters 'emax' through 'tau' are also saved as constants.
        
        All parameters are in SI units.

        PARAMETERS
        ----------
        emax : float
            Peak electric field. Defaults to 20e9.
        lam : float 
            Laser wavelength. Defaults to 800e-9.
        rad : float
            Apex radius of curvature (for collective image fields). Defaults to 20e-9.
        ef : float
            Fermi energy. Defaults to 5.51*1.602e-19.
        wf : float
            Work function. Defaults to 5.1*1.602e-19.
        tau : float
             Full-width half-max power. Defaults to 8e-15.
        data_fol : str 
            Folder to store output data in.
            May include final '/', or not.
            Defaults to "data/".
        callback 
            Callback function which must take an integer between 0 and 100, inclusive.
            callback is called during time-stepping with the current percentage complete, as an integer. 
            Defaults to None.
        target_total_truncation_error : float
            Upper bound for total relative truncation error of the wavefunction. Defaults to 0.01.
        min_emitted_energy : float
            Minumum emittable energy. 
            Sets temporal length of simulation such that a particle emitted with this energy at the peak of the laser pulse reaches the rightmost boundary before timestepping stops. 
            Defaults to 1.602e-19.
        target_elec_num : float
            Target number of states. 
            Sets the width of the Jellium slab such that there are about this many 1-D states below the Fermi level. Actual number of states may deviate.
            Defaults to 50.
        abs_width : float
            Width of the absorptive boundary. Defaults to 20e-9.
        min_timesteps : int
            Minimum number of timesteps. 
            Overrides target_total_truncation_error if the number of timesteps would be too few. 
            Defaults to 2000.
        measure_density : bool
            Whether to use the Psi2t and Vfunct measurers. Defaults to True.
    """    
    
    
    os.makedirs(data_fol, exist_ok=True)
    if data_fol[-1] != '/':
        data_fol += '/'
    ### GET SIMULATION PARAMETERS ###
    
    #simulation must hold 2 ponderomotive amplitudes, and must prevent 10Up electrons from emitting for 2*tau from their emission
    #for low-field limit, must also have 5 penetration depths at Fermi level (assuming square well)
    #   xmax = 2*ap + 2*tau*sqrt(2/m*10Up) + 5/sqrt(2*m*W/hbar^2)
    #must have enough time after 5*tau (center) for min_emitted_energy to escape
    #   duration = 5*tau + xmax / sqrt(2/m * min_emitted_energy)
    
    peak_t = 5*tau
    xmax = 2*pond_a(emax, lam) + 2*tau*np.sqrt(2.0/cons.m_e * 10.0 * pond_U(emax, lam)) + 5/np.sqrt(2.0*cons.m_e*wf)*cons.hbar
    duration = peak_t + xmax / np.sqrt(2/cons.m_e * min_emitted_energy)
    
    #well must be large enough to store target_elec_num electrons
    # nelec = sqrt(ef*8*m/(h^2)*l^2))
    # -> l = h nelec / sqrt(8 ef m)
    
    well_width = cons.h * target_elec_num / np.sqrt(8 * ef * cons.m_e)
    jell_back = well_width/4.0
    xmin = -well_width - jell_back
    
    #dx must be able to resolve maximum energy scale:
    #   20 U_p, E_f + W, 
    #dt must be chosen to get small total error
    # TE = dt^2 t / (24 hbar m) * V'^2 < target_total_truncation_error
    #   t is length of simu
    #   V' is maximum potential gradient
    #       max for Jellium (Hartrees):
    #       (E_f + W)^2/4 / (2 (E_f+W)/kf - 1)
    #       and combine with max for field
    
    efw_h = ef + wf / cons.physical_constants["Hartree energy"][0]
    kf_h = np.sqrt(ef/cons.physical_constants["Hartree energy"][0])
    max_jel_grad = efw_h**2/(4*(2*efw_h/kf_h-1)) \
        * cons.physical_constants["Hartree energy"][0] \
            * cons.physical_constants["hartree-inverse meter relationship"][0]
        
    max_energy = 10*pond_U(emax, lam)
    energy_resolution = 2*max_energy + ef + wf
    gradient_resolution = max_jel_grad + cons.e*emax

    dx = np.sqrt( cons.hbar**2/(2*cons.m_e*energy_resolution) )
    dt = np.sqrt( target_total_truncation_error * 24.0 * cons.hbar * cons.m_e / \
        (duration * (max_jel_grad + cons.e*emax)**2) )
    dt = min(duration/min_timesteps, dt)

    ### INITIALIZE SIMULATION ###

    sim = Simulation(xmin=xmin-abs_width, xmax=xmax+abs_width, dx=dx, dt=dt, maxT=duration, callback=callback)
    
    sim.setDens(DirectDensity())
    sim.setWght(FermiGasDistro(ef))

    ### INITIALIZE POTENTIALS ###

    jellPot = Potentials.JelliumPotential(sim, 0, ef, wf, -well_width, jell_back, xmax)
    fieldPot = Potentials.PulsePotential(
        sim, 
        Potentials.FileFieldProfile(sim, 0.0, xmax - 2*abs_width, xmin, abs_width, emax, "au35_cr5_si_800nm.field"),
        Potentials.CosSquaredEnvelope(tau, peak_t),
        np.pi/2, peak_t, lam, xmax)
    imagPot = Potentials.CylindricalImagePotential(
        sim, ef, wf, rad, xmin, xmax, 0.0, xmax)

    sim.addPot(jellPot)
    sim.addPot(fieldPot)
    sim.addPot(imagPot)

    ### INITIALIZE MEASURERS ###

    sim.addMeas(Measurers.Basic(sim, data_fol))
    sim.addMeas(Measurers.ExpectA(sim, data_fol))
    sim.addMeas(Measurers.NElec(sim, data_fol))
    if measure_density:
        sim.addMeas(Measurers.Psi2t(sim, 800, 800, data_fol))
        sim.addMeas(Measurers.Vfunct(sim, 800, 800, data_fol))
    sim.addMeas(Measurers.VDFluxSpec(sim, xmax, 0, 10000, 1000*1.602e-19, "vacc", data_fol))
    sim.addMeas(Measurers.Weights(sim, data_fol))
    
    sim.addMeas(Measurers.Constant(emax, "emax", data_fol))
    sim.addMeas(Measurers.Constant(lam, "lam", data_fol))
    sim.addMeas(Measurers.Constant(rad, "rad", data_fol))
    sim.addMeas(Measurers.Constant(ef, "ef", data_fol))
    sim.addMeas(Measurers.Constant(wf, "wf", data_fol))
    sim.addMeas(Measurers.Constant(tau, "tau", data_fol))

    ### SET KINETIC OPERATOR, ADD ABSORPTIVE BOUNDARIES ###

    sim.setKin(PSM_FreeElec(sim, 1.0))

    ### RUN SIMULATION ###

    print("\tFinishing initialization...")
    sim.finishInit()

    print("\tEigensolving...")
    sim.eigenSolve(-wf-ef, -wf)

    print("\tNegating self-consistent potentials...")
    imagPot.negatePotential(sim)
    
    print("\tSimulating...")
    sim.runOS_UW2TUW()

    print("\tDone!")

# /runSimulation/