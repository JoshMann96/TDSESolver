#pragma once
#include "CORECommonHeader.h"
#include "Measurers.h"
#include "Potentials.h"
#include "WfcRhoTools.h"
#include "KineticOperator.h"

// A cyclic integer class that wraps around a maximum value. Useful for managing the local wavefunction and potential history.
class cyclic_int
{
protected:
	int val, max;
public:
	cyclic_int() : val(0), max(0) {};
	cyclic_int(int max) : val(0), max(max) {};
	cyclic_int(int val, int max) : val(val), max(max) {};
	inline void increment() { val = (val + 1) % max; };
	inline cyclic_int& operator++() { increment(); return *this; }; //prefix
	inline cyclic_int operator++(int) { cyclic_int c = *this; increment(); return c; }; //postfix
	inline cyclic_int operator+(int n) { cyclic_int c(max); c.val = (val + n) % max; return c; };
	inline cyclic_int& operator+=(int n) { val = (val + n) % max; return *this; };
	inline cyclic_int& operator=(int n) { val = n % max; return *this; };
	inline operator int() const { return val; };
};

// Tool for tracking progress, calling a callback function when a certain percentage of the task is done.
class ProgressTracker
{
	private:
		std::function<void(int)> progCallback;
		int percDone = 0, nSteps = -1;
	public:
	ProgressTracker(std::function<void(int)> callback) : progCallback(callback) {};
	void update(int step) {
		if(nSteps < 0)
			throw std::runtime_error("ProgressTracker::update: Number of steps not set!");
		while (step*(long)100 / nSteps > percDone) {
			if (progCallback != nullptr)
				progCallback(percDone);
			percDone++;
		}
	};
	void reset(int nSteps) {
		this->nSteps = nSteps;
		percDone = 0;
	};
};

// Manages simulation by calling Potentials, Measuruers, and KineticOperators (with corresponding numerical methods for time integration) 
// 		using task parallelism as well as the parallelism used within each object.
// Stores a brief history of the potential and wavefunctions so that measurements can be done in parallel.
// For linear systems the potential calcualtion may be done in parallel as well.
class SimulationManager
{
private:
	Potentials::PotentialManager * pot;
	Measurers::MeasurementManager * meas;
	WfcToRho::Weight* wght = nullptr;
	WfcToRho::Density* dens = nullptr;

	KineticOperators::KineticOperator* kin;

	double *ts, dt, dx;
	double **vs, **rhos, *spatialDamp;
	int nPts, nElec, calcDensity = 0;
	cyclic_int index;
	int* step;
	std::complex<double> *scratch1, *scratch2;

	WfcToRho::NormalizationScheme normScheme = WfcToRho::UNNORMALIZED;

	int calculatePotential(double* rho, std::complex<double>* psi, double t, double* v);
	int updatePotential(int idx);
	int measure(int idx);

	std::complex<double> **psis;

	void freePsis(){
		for(int i = 0; i < 4; i++){
			if(psis[i]){
				sq_free(psis[i]);
				psis[i] = nullptr;
			}
		}
	}

	double* weights = nullptr;

	void calcWeights();

	const int HISTORY_LENGTH = 4; // TODO: See how much history is necessary for each scheme (psi/pot/mea asynchronous vs pot/mea asyncrhonous vs only mea asynchronous)

	ProgressTracker progTracker;

public:

	SimulationManager(int nPts, double dx, double dt, std::function<void(int)> callback = nullptr);
	~SimulationManager();

	// Adds a measurer to the simulation.
	void addMeasurer(Measurers::Measurer * nMeas);
	// Adds potential to simulation.
	void addPotential(Potentials::Potential * nPot);
	// Adds (rather, multiplies) spatial absorptive region.
	void addSpatialDamp(double* arr);

	// Setting functions relating to density calculation.
	void setWeight(WfcToRho::Weight* nwght) { wght = nwght; }
	void setDensity(WfcToRho::Density* ndens) { dens = ndens; }
	WfcToRho::Weight* getWeight () const { return wght; }
	WfcToRho::Density* getDensity () const { return dens; }

	void calcEnergies(int step, double* energies) const;

	double* getWeightValues() const { return weights; }
	double* getRho(int curStep) const;

	// Sets the kinetic operator to be used in the simulation.
	void setKineticOperator(KineticOperators::KineticOperator* nkin) { kin = nkin; }

	// Split-step iteration schemes
	void runOS_U2TU(int nSteps); // Assumes linear potential
	void runOS_UW2TUW(int nSteps); // Appropriate for nonlinear potentials

	// Finite difference methods
	void runFD_L(int nSteps); // Assumes linear potential
	//void runFD_NL(int nSteps); // Appropriate for nonlinear potential

	// Attemps to find ground state.
	// Same as above, using pseudospectral method (may take up a lot of memory for pseudospectral methods)
	void findEigenStates(double emin, double emax);
	// Find steady state from inhomogeneous BCs (FINITE DIFFERENCE METHODS ONLY)
	void findInhomogeneousSteadyStates_OBSOLETE(double threshold, int nElec, double* kl, double* kr, bool verbose = false);
	void findInhomogeneousEigenStates(int nElec, double* energies);
	// Sets the wave function of the simulation.
	void setPsi(std::complex<double>* npsi, WfcToRho::NormalizationScheme norm = WfcToRho::UNNORMALIZED);

	void iterateIndex();

	// Returns the number of points in the simulation.
	int getNumPoints() const {return nPts;};
	// Returns the dx or dt spacing.
	double getDX() const {return dx;};
	double getDT() const {return dt;};
	// Returns a pointer to the current psis.
	std::complex<double>* getPsi() const {return psis[index];};
	double* getRho() const {return rhos[index];};
	int getNElec() const {return nElec;};
	int* getNElecPtr() {return &nElec;};

	bool canAsyncCalcPot() const { return pot->getComplexity() != Potentials::PotentialComplexity::WAVEFUNCTION_DEPENDENT; }

	double** getWeightsPtr() { return &weights; }
	KineticOperators::KineticOperator** getKin() { return &kin; }

	Potentials::Potential* getPotPointer() const { return pot; }

	int findElectricalSurfaceCentroidRule(int minPos, int maxPos);
};