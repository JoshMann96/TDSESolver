import sys, os
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)

from tdsepy.tdsepy import *
import numpy as np
import scipy.constants as cons
import json
import time
#from multiprocess import Pool

def pond_U(emax, lam):
    return cons.e**2 * emax**2 * lam**2 / (4* cons.m_e * (2*cons.pi*cons.c)**2)

def pond_a(emax, lam):
    return cons.e * emax * lam**2 / (cons.m_e * (2*cons.pi*cons.c)**2)

"""
def runSimSweepFieldsMPI(emaxs:list, lam:float=800e-9, rad:float=20e-9, ef:float=5.51*1.602e-19, wf:float=5.1*1.602e-19, tau:float=8e-15, data_fol:str="data/", callback=None, 
                  target_total_truncation_error:float = 0.01, min_emitted_energy:float=1.602e-19, target_elec_num:float = 50, abs_width:float=20e-9, abs_rate:float=5e17, min_timesteps:int=2000, measure_density:bool=True):
""""""Runs a series of rescattering simulations within a range of peak field strengths.
        Runs SLURM_NTASKS simulations at a time.
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
        abs_rate : float
            Rate of absorption. Defaults to 5e17.
        min_timesteps : int
            Minimum number of timesteps. 
            Overrides target_total_truncation_error if the number of timesteps would be too few. 
            Defaults to 2000.
        measure_density : bool
            Whether to use the Psi2t and Vfunct measurers. Defaults to True.
""""""   
    
    os.makedirs(data_fol, exist_ok=True)
    if data_fol[-1] != '/':
        data_fol += '/'
    
    emaxs = list(emaxs)
    
    paramDict = locals()
    del paramDict["callback"]
    if "runSim" in paramDict.keys():
        del paramDict["runSim"]
    
    with open(f"{data_fol}params.json", 'w') as fil:
        json.dump(paramDict, fil)
    
    runSim = lambda i, emax : runSimTimed(i, emax, lam, rad, ef, wf, tau, data_fol, callback, target_total_truncation_error, min_emitted_energy, target_elec_num, abs_width, abs_rate, min_timesteps, measure_density)
    with Pool(int(os.environ['SLURM_NTASKS'])) as pool:
        pool.starmap(runSim, enumerate(emaxs))
"""

def runSimTimed(i, emax, lam, rad, ef, wf, tau, data_fol, callback, target_total_truncation_error, min_emitted_energy, target_elec_num, abs_width, abs_rate, min_timesteps, measure_density):
        strt = time.time()
        runSingleSimulation(emax, lam, rad, ef, wf, tau, f"{data_fol}{i}/", callback, target_total_truncation_error, min_emitted_energy, target_elec_num, abs_width, abs_rate, min_timesteps, measure_density)
        stop = time.time()
        dur = stop-strt
        
        if dur > 3600:
            print(f"{i:3d}:{dur/3600:.2f} h")
        else:
            print(f'{i:3d}{dur:.1f} s')
        print() 

def runSimSweepFields(emaxs:list, lam:float=800e-9, rad:float=20e-9, ef:float=5.51*1.602e-19, wf:float=5.1*1.602e-19, tau:float=8e-15, data_fol:str="data/", callback=None, 
                  target_total_truncation_error:float = 0.01, min_emitted_energy:float=1.602e-19, target_elec_num:float = 50, abs_width:float=20e-9, abs_rate:float=5e17, min_timesteps:int=2000, measure_density:bool=True):
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
        abs_rate : float
            Rate of absorption. Defaults to 5e17.
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
        print(f"{i:3d}/{len(emaxs):3d}")
        strt = time.time()
        runSingleSimulation(emax, lam, rad, ef, wf, tau, f"{data_fol}{i}/", callback, target_total_truncation_error, min_emitted_energy, target_elec_num, abs_width, abs_rate, min_timesteps, measure_density)
        stop = time.time()
        if callback is not None:
            print()
        
        dur = stop-strt
        if dur > 3600:
            print(f"{dur/3600:.2f} h")
        else:
            print(f'{dur:.2f} s')
        print()
    
    
def runSingleSimulation(emax:float=20e9, lam:float=800e-9, rad:float=20e-9, ef:float=5.51*1.602e-19, wf:float=5.1*1.602e-19, tau:float=8e-15, cep:float=np.pi/2, data_fol:str="data/", callback=None, 
                  target_total_truncation_error:float = 0.01, min_emitted_energy:float=1.602e-19, target_elec_num:float = 50, abs_width:float=20e-9, abs_rate:float=5e17, min_timesteps:int=2000,
                  exchange_correlation=True, bulk_hartree=True, measure_density:bool=True, verbose:bool=False):
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
        cep : float
            Carrier envelope phase. Defaults to pi/2.
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
        abs_rate : float
            Rate of absorption. Defaults to 5e17.
        exchange_correlation : bool
            Whether to use exchange-correlation potentials. Defaults to True.
        bulk_hartree : bool
            Whether to use bulk Hartree potentials. Defaults to True.
        min_timesteps : int
            Minimum number of timesteps. 
            Overrides target_total_truncation_error if the number of timesteps would be too few. 
            Defaults to 2000.
        measure_density : bool
            Whether to use the Psi2t and Vfunct measurers. Defaults to True.
    """    
    
    
    if verbose:
        message = print
    else:
        message = lambda *args: None
    
    os.makedirs(data_fol, exist_ok=True)
    if data_fol[-1] != '/':
        data_fol += '/'

    ### GET SIMULATION PARAMETERS ###
    peak_t = 5*tau
    xmin, xmax, duration, dx, dt, max_energy, well_width, jell_back = getSimulationParameters(emax, lam, ef, wf, tau, peak_t, min_emitted_energy, target_elec_num, target_total_truncation_error, min_timesteps)

    ### INITIALIZE SIMULATION ###

    sim = Simulation(xmin=xmin-abs_width, xmax=xmax+abs_width, dx=dx, dt=dt, maxT=duration, callback=callback)

    dens = CylindricalDensity(-rad, rad, xmin-abs_width)
    dens.setBaseDensity(DirectDensity())
    sim.setDens(dens)
    sim.setWght(FermiGasDistro(ef))

    ### INITIALIZE POTENTIALS ###

    jellPot = Potentials.JelliumPotential(sim, 0, ef, wf, -well_width, jell_back, xmax)
    fieldPot = Potentials.PulsePotential(
        sim, 
        Potentials.FileFieldProfile(sim, 0.0, xmax, xmin, abs_width, emax, "au35_cr5_si_800nm.field"),
        Potentials.CosSquaredEnvelope(tau, peak_t),
        cep, peak_t, lam, xmax)
    
    sim.addPot(jellPot)
    sim.addPot(fieldPot)
    
    if exchange_correlation:
        xPot = Potentials.LDAFunctional(sim, Potentials.X_SLATER, xmax)
        cPot = Potentials.LDAFunctional(sim, Potentials.C_PW, xmax)
        sim.addPot(xPot)
        sim.addPot(cPot)
    
    if bulk_hartree:
        selfPot = Potentials.PlanarToCylindricalHartreePotential(sim, rad, xmin, xmax, 0.0, xmax)
    else:
        selfPot = Potentials.CylindricalImagePotential(sim, ef, wf, rad, xmin, xmax, 0.0, xmax)
        
    sim.addPot(selfPot)

    ### INITIALIZE MEASURERS ###

    sim.addMeas(Measurers.Basic(sim, data_fol))
    sim.addMeas(Measurers.ExpectA(sim, data_fol))
    sim.addMeas(Measurers.NElec(sim, data_fol))
    sim.addMeas(Measurers.ExpectE0(sim, data_fol))
    if measure_density:
        sim.addMeas(Measurers.Psi2t(sim, 800, 800, data_fol))
        sim.addMeas(Measurers.Vfunct(sim, 800, 800, data_fol))
    sim.addMeas(Measurers.VDFluxSpec(sim, xmax, 0, 10000, max_energy*1.5 + (ef+wf), "vacc", data_fol))
    sim.addMeas(Measurers.Weights(sim, data_fol))
    
    sim.addMeas(Measurers.Constant(emax, "emax", data_fol))
    sim.addMeas(Measurers.Constant(lam, "lam", data_fol))
    sim.addMeas(Measurers.Constant(rad, "rad", data_fol))
    sim.addMeas(Measurers.Constant(ef, "ef", data_fol))
    sim.addMeas(Measurers.Constant(wf, "wf", data_fol))
    sim.addMeas(Measurers.Constant(tau, "tau", data_fol))
    sim.addMeas(Measurers.Constant(cep, "cep", data_fol))
    sim.addMeas(Measurers.Constant(abs_width, "abs_width", data_fol))
    sim.addMeas(Measurers.Constant(abs_rate, "abs_rate", data_fol))

    ### SET KINETIC OPERATOR, ADD ABSORPTIVE BOUNDARIES ###

    sim.setKin(PSM_FreeElec(sim, 1.0))
    
    sim.addLeftAbsBdy(abs_rate, abs_width)
    sim.addRightAbsBdy(abs_rate, abs_width)

    ### RUN SIMULATION ###

    message("\tFinishing initialization...")
    sim.finishInit()

    message("\tEigensolving...")
    sim.eigenSolve(-wf-ef, -wf)

    message("\tNegating self-consistent potentials...")
    selfPot.negatePotential(sim)
    if exchange_correlation:
        xPot.negatePotential(sim)
        cPot.negatePotential(sim)
    
    message("\tSimulating...")
    sim.runOS_UW2TUW()

    message("\tDone!")

# /runSimulation/

def getSimulationParameters(emax:float, lam:float, ef:float, wf:float, tau:float, peak_t:float, min_emitted_energy:float=1.602e-19, target_elec_num:int=50, target_total_truncation_error:float=0.01, min_timesteps:int=2000):
    """Calculates desirable simulation parameters for a rescattering simulation.

    Args:
        emax (float): Electric field amplitude.
        lam (float): Wavelength.
        ef (float): Fermi energy.
        wf (float): Work function.
        tau (float): Full-width half-max power.
        peak_t (float): Time of envelope maximum.
        min_emitted_energy (float, optional): Minimum emitted energy to wait for. Defaults to 1.602e-19.
        target_elec_num (int, optional): Target number of eigenstates. Defaults to 50.
        target_total_truncation_error (float, optional): Target total truncation error. Defaults to 0.01.
        min_timesteps (int, optional): Minimum timesteps. Defaults to 2000.

    Returns:
        xmin : Minimum x position.
        xmax : Maximum x position.
        duration : Duration of simulation.
        dx : x grid spacing.
        dt : Time step.
        max_energy : Maximum emitted energy (10 Up).
        well_width : Width of Jellium slab.
        jell_back : Width of Jellium backing.
    """
    
    ### GET SIMULATION PARAMETERS ###
    
    #simulation must hold 2 ponderomotive amplitudes, and must prevent 10Up electrons from emitting for 2*tau from their emission
    #for low-field limit, must also have 5 penetration depths at Fermi level (assuming square well)
    #   xmax = 2*ap + 2*tau*sqrt(2/m*10Up) + 5/sqrt(2*m*W/hbar^2)
    #must have enough time after 5*tau (center) for min_emitted_energy to escape
    #   duration = 5*tau + xmax / sqrt(2/m * min_emitted_energy)
    
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
        (duration * (gradient_resolution)**2) )
    dt = min(duration/min_timesteps, dt)
    
    return xmin, xmax, duration, dx, dt, max_energy, well_width, jell_back