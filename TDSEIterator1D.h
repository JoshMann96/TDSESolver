#pragma once
// Solves the TDSE by iteration.
class TDSEIterator1D
{
private:
	std::complex<double> cnAl, cnBe, eulAl, eulBe, mpAl, mpBe;
	std::complex<double> *cnb, *cncp, *cnd, *cndp;
	std::complex<double> *dl, *d, *du, *du2, *a, *x, *outData, *temp, *dl0, *d0, *du0;
	std::complex<double> *timePhase;
	std::complex<double>* spatialDamp;
	std::future<int> *fs = new std::future<int>[16];
	int *ia, *ja, *perm;
	int *iparm, *pt, *ipiv;
	int nPts;
	int mtype = 3, maxfct = 1, mnum = 1, phase = 13, nrhs = 1, msglvl = 0, error = 0;
	double dx, dt;
	double h2al, h2be, h2dt;
	ThreadPool* pool;
	std::vector<std::future<int>> res;
	int first_stepOS_U2TU = 1;
	DFTI_DESCRIPTOR_HANDLE dftiHandle = 0;
	std::complex<double>* osKineticPhase, *osPotentialPhase;
	double* tempD;
public:
	// Initializes the TDSE iterator.
	TDSEIterator1D(double dx, double dt, int simSize);
	// Initializes the iterator for the Crank-Nicolson method.
	void initializeCN();
	// OBSOLETE
	void initializeCNPA_OBS();
	// Adds time phase to parts of solver for absorption.
	void addTimePhase(std::complex<double> * arr);
	// Add spatial damping for absorption.
	void addSpatialDamp(double* arr);
	// Steps using Crank-Nicolson method (Intel MKL).
	void stepCN(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ);
	void stepCN(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ, double ndt);
	void multiStepCN(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ, int nelec);
	void stepOS_U2TU(std::complex<double>* __restrict psi0, double* __restrict v0, std::complex<double>* targ, int nelec);
	void stepOS_UW2T(std::complex<double>* __restrict psi0, double* __restrict v0, std::complex<double>* targ, int nelec);
	void stepOS_UW(std::complex<double>* __restrict psi0, double* __restrict v0, std::complex<double>* targ, int nelec);
	// Uses iterative refinement to obtain a less noisy solution -- very slow (~7-8 times slower), and simulation is limited by noise from initial state anyway.
	void stepCNREF(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ);
	void stepCNREF(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ, double ndt);
	void stepH2Euler(std::complex<double> *__restrict psi0, double *__restrict v0, std::complex<double> * targ);
	// Steps using Euler method.
	void stepEuler(std::complex<double> * psi0, double * v, std::complex<double> * targ);
	// Steps using CN method in imaginary time.
	void stepCNImT(std::complex<double> *__restrict psi0, double *__restrict v, std::complex<double> * targ, double rate);
	// Steps using Euler method in imaginary time.
	void stepEulerImT(std::complex<double> * psi0, double * v, std::complex<double> * targ);
	// OBSOLETE
	void stepCNPA_OBS(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ);
	// Steps using Crank-Nicolson by directly implementing the Thomas algorithm (sketchy).
	void stepCNSL(std::complex<double> * psi0, double * v0, double * vf, std::complex<double> * targ);
};