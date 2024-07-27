#pragma once
#include "CORECommonHeader.h"
#include "Measurers.h"
#include "Potentials.h"
#include "KineticOperator.h"
// Manages simulation by controling potentials, measurements, and TDSE iterator(s) for multiple electrons at the same time.

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

	double *ts, maxT, dt, dx;
	double **vs, **rhos, *spatialDamp;
	int nPts, index, nElec, numSteps, calcDensity = 0;
	int* step;
	std::function <void(int)> progCallback;
	std::complex<double> *scratch1, *scratch2;

	int updatePotential(std::complex<double>* psi, int idx, double* rho);
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

public:
	SimulationManager(int nPts, double dx, double dt, double maxT, std::function<void(int)> callback);
	~SimulationManager();
	// Adds a measurer to the simulation.
	void addMeasurer(Measurers::Measurer * nMeas);
	// Adds potential to simulation.
	void addPotential(Potentials::Potential * nPot);
	// Adds spatial absorptive region.
	void addSpatialDamp(double* arr);

	void setWeight(WfcToRho::Weight* nwght) { wght = nwght; }
	void setDensity(WfcToRho::Density* ndens) { dens = ndens; }
	WfcToRho::Weight* getWeight() { return wght; }
	WfcToRho::Density* getDensity() { return dens; }

	void calcEnergies(int step, double* energies);

	double* getWeightValues(){ return weights; }
	double* getRho(int curStep);

	void setKineticOperator_PSM(KineticOperators::KineticOperator_PSM* nkin) { kin = nkin; kin_psm = nkin; }
	void setKineticOperator_FDM(KineticOperators::KineticOperator_FDM* nkin) { kin = nkin; kin_fdm = nkin; }

	void runOS_U2TU();
	void runOS_UW2TUW();

	// Finishes initialization of manager (REQUIRED BEFORE RUNNING).
	void finishInitialization();
	// Attemps to find ground state.
	// Same as above, using pseudospectral method (may take up a lot of memory)
	void findEigenStates(double fermie, double w, double maxT, double rate);
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
	int getNumSteps();
	// Returns a pointer to the current psis.
	std::complex<double> * getPsi();
	double * getRho();
	int getNElec();
	int* getNElecPtr();
	double** getWeightsPtr() { return &weights; }
	KineticOperators::KineticOperator** getKin() { return &kin; }

	Potentials::Potential* getPotPointer() { return pot; }
};