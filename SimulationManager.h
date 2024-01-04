#pragma once

class SimulationManager {
public:
	virtual void addMeasurer(Measurers::Measurer * nMeas) = 0;
	// Adds potential to simulation.
	virtual void addPotential(Potentials::Potential * nPot) = 0;
	// Adds absorptive region to potential (OBSOLETE).
	virtual void addAbsorptiveRegion(AbsorptiveRegions::AbsorptiveRegion * nReg) = 0;
	// Adds imaginary time rotation to kinetic energy portion of hamiltonian (absorptive).
	virtual void addTimeRotation(std::complex<double> * arr) = 0;
	// Adds spatial absorptive region.
	virtual void addSpatialDamp(double* arr) = 0;
	// Runs simulation.
	virtual void run() = 0;
	virtual void runLogProgress(int idx) = 0;
	// Times individual components of simulation.
	virtual void runTimed() = 0;
	// Runs simulation and prints progress.
	virtual void runPrintProgress() = 0;
	// Runs the various parts of the simulation in parallel, and prints progress. Keep in mind that with psi-dependent potentials that a 2-step old psi is used to do the calculation.
	virtual void runParallel() = 0;
	virtual void runParallelLogProgress(int idx) = 0;
	virtual void runSCFLogProgress(double thresh) = 0;
	virtual void runPertLogProgress(int idx, int nit) = 0;
	virtual void runOS_U2TULogProgress(int idx) = 0;
	virtual void runOS_UW2TUWLogProgress(int idx) = 0;
	// OBSOLETE
	virtual void runParallelOMPOld() = 0;
	// Finishes initialization of manager (REQUIRED BEFORE RUNNING).
	virtual void finishInitialization() = 0;
	// Attemps to find ground state.
	virtual void findGroundState(int nSteps, double prec) = 0;
	virtual void findGroundStatePrintProgress(int nSteps, double prec) = 0;
	// Attempts to find ground state including a wave function dependent potential.
	virtual void findGroundStateWithWaveFuncDepPot(int nSteps, double prec) = 0;
	virtual void findEigenState(double energ, double tolerance) = 0;
	virtual void findMetallicInitialState(double fermie, double w, double maxT, double rate) = 0;
	virtual void findMetallicInitialState_PSM(double fermie, double w, double maxT, double rate) = 0;
	virtual void findMetallicInitialState_HOD(double fermie, double w, double maxT, double rate, int order) = 0;
	// Sets the wave function of the simulation.
	virtual void setPsi(std::complex<double>* npsi) = 0;

	// Returns the number of points in the simulation.
	virtual int getNumPoints() = 0;
	// Returns the dx or dt spacing.
	virtual double getDX() = 0;
	virtual double getDT() = 0;
	virtual double getMaxT() = 0;
	// Returns a pointer to the current psi.
	virtual std::complex<double> * getPsi() = 0;
	virtual int getNElec() = 0;
	// Returns a pointer to the potential manager being used.
	virtual Potentials::Potential * getPotPointer() = 0;
};