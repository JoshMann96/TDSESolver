#pragma once
#include "CORECommonHeader.h"
#include "Measurers.h"
#include "Potentials.h"
#include "WfcRhoTools.h"
#include "KineticOperator.h"

// Manages simulation by controling potentials, measurements, and TDSE iterator(s) for multiple electrons at the same time.
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

class SimulationManager
{
private:
	Potentials::PotentialManager * pot;
	Measurers::MeasurementManager * meas;
	WfcToRho::Weight* wght = nullptr;
	WfcToRho::Density* dens = nullptr;

	KineticOperators::KineticOperator* kin;
	KineticOperators::KineticOperator_PSM * kin_psm;
	KineticOperators::KineticOperator_FDM* kin_fdm;

	double *ts, dt, dx;
	double **vs, **rhos, *spatialDamp;
	int nPts, index, nElec, calcDensity = 0;
	int* step;
	std::complex<double> *scratch1, *scratch2;

	WfcToRho::NormalizationScheme normScheme = WfcToRho::UNNORMALIZED;

	int calculatePotential(double* rho, std::complex<double>* psi, double t, double* v);
	int updatePotential(int idx);
	int stepItPAR(int idx0, int idx1);
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
	WfcToRho::Weight* getWeight() { return wght; }
	WfcToRho::Density* getDensity() { return dens; }

	void calcEnergies(int step, double* energies);

	double* getWeightValues(){ return weights; }
	double* getRho(int curStep);

	// Sets the kinetic operator to be used in the simulation. Differentiates between pseudospectral and finite difference methods.
	void setKineticOperator_PSM(KineticOperators::KineticOperator_PSM* nkin) { kin = nkin; kin_psm = nkin; }
	void setKineticOperator_FDM(KineticOperators::KineticOperator_FDM* nkin) { kin = nkin; kin_fdm = nkin; }

	// Split-step iteration schemes
	void runOS_U2TU(int nSteps); // Assumes linear potential
	// Attemps to find ground state.
	// Same as above, using pseudospectral method (may take up a lot of memory for pseudospectral methods)
	void findEigenStates(double fermie, double w, double maxT, double rate);
	// Find steady state from inhomogeneous BCs (FINITE DIFFERENCE METHODS ONLY)
	void findInhomogeneousSteadyStates_OBSOLETE(double threshold, int nElec, double* kl, double* kr, bool verbose = false);
	void findInhomogeneousEigenStates(int nElec, double* energies);
	// Sets the wave function of the simulation.
	void setPsi(std::complex<double>* npsi, WfcToRho::NormalizationScheme norm = WfcToRho::UNNORMALIZED);

	void iterateIndex();
	int getIndex();
	int getNextIndex();
	int getPrevIndex();
	int getPrevPrevIndex();

	// Returns the number of points in the simulation.
	int getNumPoints();
	// Returns the dx or dt spacing.
	double getDX();
	double getDT();
	// Returns a pointer to the current psis.
	std::complex<double>* getPsi();
	double* getRho();
	int getNElec();
	int* getNElecPtr();
	double** getWeightsPtr() { return &weights; }
	KineticOperators::KineticOperator** getKin() { return &kin; }

	Potentials::Potential* getPotPointer() { return pot; }

	int findElectricalSurfaceCentroidRule(int minPos, int maxPos);
};