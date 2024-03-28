import sys, os
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)

from tdsepy import *
import numpy as np

for field in range(5,80+1,5):

    MAX_X = 150e-9
    MIN_X = -60e-9
    ABS_S = 15e-9
    DUR   = 100e-15
    MAX_E = field*1e9
    PEAK_T= 50e-15
    DATA_FOL = f"data/{field:02d}/"
    WORKF = 5.1*1.602e-19
    FERMIE= 5.53*1.602e-19

    os.makedirs(DATA_FOL, exist_ok=True)

    ### INITIALIZE SIMULATION ###

    print(f"E = {field:02d} V/nm : ")
    print("\tBuilding simulation...")

    sim = Simulation(xmin=MIN_X, xmax=MAX_X, dx=20e-12, dt=10e-19, maxT=DUR, callback=lambda x : print(f"{x:02d} ", end=''))
    
    sim.setDens(DirectDensity())
    sim.setWght(FermiGasDistro(5.53*1.602e-19))

    ### INITIALIZE POTENTIALS ###

    jellPot = Potentials.JelliumPotential(sim, 0, FERMIE, WORKF, MIN_X/2, 10e-9, MAX_X-ABS_S)
    fieldPot = Potentials.PulsePotential(
        sim, 
        Potentials.FileFieldProfile(sim, 0.0, MAX_X-ABS_S, MIN_X, 10e-9, MAX_E, "au35_cr5_si_800nm.field"),
        Potentials.CosSquaredEnvelope(8e-15, PEAK_T),
        np.pi/2, PEAK_T, 800e-9, MAX_X-ABS_S)
    imagPot = Potentials.CylindricalImagePotential(
        sim, FERMIE, WORKF, 20e-9, MIN_X+ABS_S, MAX_X-ABS_S, 0.0, MAX_X-ABS_S)

    sim.addPot(jellPot)
    sim.addPot(fieldPot)
    sim.addPot(imagPot)

    ### INITIALIZE MEASURERS ###

    sim.addMeas(Measurers.Basic(sim, DATA_FOL))
    sim.addMeas(Measurers.ExpectA(sim, DATA_FOL))
    sim.addMeas(Measurers.NElec(sim, DATA_FOL))
    #sim.addMeas(Measurers.Psi2t(sim, 800, 800, DATA_FOL))
    #sim.addMeas(Measurers.Vfunct(sim, 800, 800, DATA_FOL))
    sim.addMeas(Measurers.VDFluxSpec(sim, MAX_X-ABS_S, 0, 10000, 1000*1.602e-19, "vacc", DATA_FOL))
    sim.addMeas(Measurers.Weights(sim, DATA_FOL))

    sim.setKin(PSM_FreeElec(sim, 1.0))

    ### RUN SIMULATION ###

    print("\tFinishing initialization...")
    sim.finishInit()

    print("\tEigensolving...")
    sim.eigenSolve(-WORKF-FERMIE, -WORKF)

    print("\tNegating self-consistent potentials...")
    imagPot.negatePotential(sim)
    
    print("\tSimulating...")
    sim.runOS_UW2TUW()

    print("\tDone!")

# /for

print("Finished Python")

#print(sim)

#print(pot.getV())