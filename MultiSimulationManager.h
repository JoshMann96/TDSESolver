#pragma once
// Manages simulation by controling potentials, measurements, and TDSE iterator(s) for multiple electrons at the same time.

class MultiSimulationManager :
	public SimulationManager {
private:
	Potentials::PotentialManager * pot;
	Measurers::MeasurementManager * meas;
	AbsorptiveRegions::AbsorptiveRegionManager * reg;
	TDSEIterator1D * it;
	double *ts, maxT, dt, dx;
	double **vs;
	int nPts, index, nelec;
	std::complex<double> **psis, *scratch1, *scratch2, *tpsi;
	double getTotalEnergy(std::complex<double> * psi, double * v);

	int getVPAR(int idx, int idxPsi);
	int stepItPAR(int idx0, int idx1);
	int measPAR(int idx);

public:
	MultiSimulationManager(int nPts, double dx, double dt, double max);
	~MultiSimulationManager();
	// Adds a measurer to the simulation.
	void addMeasurer(Measurers::Measurer * nMeas);
	// Adds potential to simulation.
	void addPotential(Potentials::Potential * nPot);
	// Adds absorptive region to potential (OBSOLETE).
	void addAbsorptiveRegion(AbsorptiveRegions::AbsorptiveRegion * nReg);
	// Adds imaginary time rotation to kinetic energy portion of hamiltonian (absorptive).
	void addTimeRotation(std::complex<double> * arr);
	// Adds spatial absorptive region.
	void addSpatialDamp(double* arr);
	// Runs simulation.
	void run();
	void runLogProgress(ProgressTracker * prg, int idx);
	// Times individual components of simulation.
	void runTimed();
	// Runs simulation and prints progress.
	void runPrintProgress();
	// Runs the various parts of the simulation in parallel, and prints progress. Keep in mind that with psi-dependent potentials that a 2-step old psi is used to do the calculation.
	void runParallel();
	void runParallelLogProgress(ProgressTracker * prg, int idx);
	void runSCFLogProgress(ProgressTracker * prg, int idx, double thresh);
	void runPertLogProgress(ProgressTracker* prg, int idx, int nit);
	void runOS_U2TULogProgress(ProgressTracker* prg, int idx);
	void runOS_UW2TUWLogProgress(ProgressTracker* prg, int idx);
	// OBSOLETE
	void runParallelOMPOld();
	// Finishes initialization of manager (REQUIRED BEFORE RUNNING).
	void finishInitialization();
	// Attemps to find ground state.
	void findGroundState(int nSteps, double prec);
	void findGroundStatePrintProgress(int nSteps, double prec);
	// Attempts to find single electronground state including a wave function dependent potential.
	void findGroundStateWithWaveFuncDepPot(int nSteps, double prec);
	void findEigenState(double energ, double tolerance);
	// Finds all states in the metal, sets initial wave functions to these
	void findMetallicInitialState(double fermie, double w, double maxT, double rate);
	// Same as above, using pseudospectral method (may take up a lot of memory)
	void findMetallicInitialState_PSM(double fermie, double w, double maxT, double rate);
	// Same as above, using higher order representation of second derivative (order 1 is 3-point, 2 is 5-point, etc)
	void findMetallicInitialState_HOD(double fermie, double w, double maxT, double rate, int order);
	// Sets the wave function of the simulation.
	void setPsi(std::complex<double>* npsi);

	void iterateIndex();
	int nextIndex();
	int prevIndex();
	int prevPrevIndex();

	// Returns the number of points in the simulation.
	int getNumPoints();
	// Returns the dx or dt spacing.
	double getDX();
	double getDT();
	double getMaxT();
	// Returns a pointer to the current psis.
	std::complex<double> * getPsi();
	int getNElec();

	Potentials::Potential* getPotPointer() { return pot; }
};