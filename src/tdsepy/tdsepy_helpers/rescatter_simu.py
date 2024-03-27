import sys, os
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)

try :
    from tdsepy import *
except ImportError:
    from ..lib.tdsepy import *
import numpy as np
import scipy.constants as cons

def pond_U(emax, lam):
    return cons.e**2 * emax**2 * lam**2 / (4* cons.m_e * (2*cons.pi*cons.c)**2)

def pond_a(emax, lam):
    return cons.e * emax * lam**2 / (cons.m_e * (2*cons.pi*cons.c)**2)

def runSimulation(emax=20e9, lam=800e-9, rad=20e-9, ef=5.51*1.602e-19, wf=5.1*1.602e-19, tau=8*1e-15, data_fol="data/", callback=None, 
                  target_total_truncation_error = 0.01, min_emitted_energy=1.602e-19, target_elec_num=50, abs_width=20e-9, min_timesteps=2000, measure_density=True):
    os.makedirs(data_fol, exist_ok=True)
    
    ### GET SIMULATION PARAMETERS ###
    
    #simulation must hold 2 ponderomotive amplitudes, and must prevent 10Up electrons from emitting for 2*tau from their emission
    #   xmax = 2*ap + 2*tau*sqrt(2/m*10Up)
    #must have enough time after 5*tau (center) for min_emitted_energy to escape
    #   duration = 5*tau + xmax / sqrt(2/m * min_emitted_energy)
    
    peak_t = 5*tau
    xmax = 2*pond_a(emax, lam) + 2*tau*np.sqrt(2.0/cons.m_e * 10.0 * pond_U(emax, lam))
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
    # error = dt^2 t / (24 hbar m) * V'^2 < 0.01
    #   t is length of simu
    #   V' is maximum potential gradient
    #       max for Jellium (Hartrees):
    #       (E_f + W)^2/4 / (2 (E_f+W)/kf - 1)
    #       and combine with max for field
    
    efw_h = ef + wf / cons.physical_constants["Hartree energy"][0]
    kf_h = np.sqrt(ef/cons.physical_constants["Hartree energy"][0])
    max_jel_grad = efw_h**2/(4*(2*efw_h/kf_h-1)) \
        * cons.physical_constants["Hartree energy"][0] * cons.physical_constants["hartree-inverse meter relationship"][0]
        
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
