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
	KineticOperators::KineticOperator* kin;
	KineticOperators::KineticOperator_PSM * kin_psm;
	KineticOperators::KineticOperator_FDM* kin_fdm;
	double *ts, maxT, dt, dx;
	double **vs, *spatialDamp;
	int nPts, index, nelec;
	std::function <int(int)> progCallback;
	std::complex<double> **psis, *scratch1, *scratch2, *tpsi;
	double getTotalEnergy(std::complex<double> * psi, double * v);

	int getVPAR(int idx, int idxPsi);
	int stepItPAR(int idx0, int idx1);
	int measPAR(int idx);

public:
	SimulationManager(int nPts, double dx, double dt, double max, std::function<int(int)> callback);
	~SimulationManager();
	// Adds a measurer to the simulation.
	void addMeasurer(Measurers::Measurer * nMeas);
	// Adds potential to simulation.
	void addPotential(Potentials::Potential * nPot);
	// Adds spatial absorptive region.
	void addSpatialDamp(double* arr);

	void setKineticOperator_PSM(KineticOperators::KineticOperator_PSM* nkin) { kin = nkin; kin_psm = nkin; }
	void setKineticOperator_FDM(KineticOperators::KineticOperator_FDM* nkin) { kin = nkin; kin_fdm = nkin; }

	void runOS_U2TU(int idx);
	void runOS_UW2TUW(int idx);

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
	// Returns a pointer to the current psis.
	std::complex<double> * getPsi();
	int getNElec();
	KineticOperators::KineticOperator* getKin() { return kin; }

	Potentials::Potential* getPotPointer() { return pot; }
};