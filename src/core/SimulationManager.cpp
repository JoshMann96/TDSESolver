#include "SimulationManager.h"
#include "PhysCon.h"
#include "MathTools.h"

//callback sends progress int 0-100 (can be nullptr for no callback)
SimulationManager::SimulationManager(int nPts, double dx, double dt, double maxT, std::function<void(int)> callback)
	: maxT(maxT), dx(dx), nPts(nPts), numSteps(std::ceil(maxT / dt))
{
	SimulationManager::dt = maxT / (numSteps - 1);

	pot = new Potentials::PotentialManager(nPts);
	meas = new Measurers::MeasurementManager("");
	psis = (std::complex<double>**) sq_malloc(sizeof(std::complex<double>*)*4);
	for(int i = 0; i < 4; i++)
		psis[i] = nullptr;

	vs = (double**) sq_malloc(sizeof(double*)*4);
	rhos = (double**) sq_malloc(sizeof(double*)*4);
	ts = (double*) sq_malloc(sizeof(double)*4);
	for (int i = 0; i < 4; i++){
		vs[i] = (double*) sq_malloc(sizeof(double) * nPts);
		rhos[i] = (double*) sq_malloc(sizeof(double) * nPts);
	}
	std::fill_n(ts, 4, 0.0);

	step = (int*) sq_malloc(sizeof(int) * 4);
	std::fill_n(step, 4, 0);

	scratch1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
	scratch2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
	spatialDamp = (double*) sq_malloc(sizeof(double) * nPts);
	std::fill_n(spatialDamp, nPts, 1.0);
	/*
	SimulationManager::maxT = maxT; SimulationManager::dt = dt; SimulationManager::dx = dx; SimulationManager::nPts = nPts;
	SimulationManager::mpiRoot = mpiRoot; SimulationManager::mpiUpdateTag = mpiUpdateTag; SimulationManager::mpiJob = mpiJob;
	*/

	progCallback = callback;

	index = 0;
	nElec = 0;
}

SimulationManager::~SimulationManager()
{
	for(int i = 0; i < 4; i++){
		sq_free(vs[i]);
		sq_free(rhos[i]);
	}
	sq_free(vs);
	sq_free(rhos);
	
	freePsis();
	sq_free(psis);

	sq_free(ts);

	sq_free(scratch1);
	sq_free(scratch2);
	sq_free(spatialDamp);
	
	if(weights)
		sq_free(weights);

	delete meas;
	delete pot;	
}

void SimulationManager::addMeasurer(Measurers::Measurer* m) {
	meas->addMeasurer(m);
}

void SimulationManager::addPotential(Potentials::Potential* p) {
	pot->addPotential(p);
	if(p->getComplexity() == Potentials::PotentialComplexity::WAVEFUNCTION_DEPENDENT)
		calcDensity = 1;
}

void SimulationManager::addSpatialDamp(double* arr) {
	vtls::seqMulArrays(nPts, arr, spatialDamp);
}

void SimulationManager::finishInitialization() {
	pot->finishAddingPotentials();
}

void SimulationManager::calcEnergies(int curStep, double* energies) {
		for(int i = 0; i < 4; i++){
			if(curStep == step[i]){ //look for the present step's index
				double* rho = (double*) sq_malloc(sizeof(double)*nPts);
				for(int j = 0; j < nElec; j++){
					vtls::normSqr(nPts, &psis[i][j*nPts], rho);
					energies[j] = vtlsInt::rSumMul(nPts, rho, vs[i], dx) + kin->evaluateKineticEnergy(&psis[i][j*nPts]);
					//potential energy + kinetic energy
				}
				sq_free(rho);

				return;
			}
		}

		throw std::runtime_error("SimulationManager::calcEnergies: Step not found!");
}

void SimulationManager::calcWeights(){
	if (nElec < 1)
		throw std::runtime_error("SimulationManager::calcEnergies: Number of electrons is not finite! Failed to initialize.");

	if(weights)
		sq_free(weights); weights = nullptr;
	weights = (double*) sq_malloc(sizeof(double)*nElec);
	double* energies = (double*) sq_malloc(sizeof(double)*nElec);

	calcEnergies(step[index], energies);

	wght->calcWeights(nElec, energies, weights);

	sq_free(energies);
}

void SimulationManager::findEigenStates(double emin, double emax, double maxT, double rate) {
	pot->getVBare(0.0, vs[index]);

	std::complex<double>* states = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nPts);

	kin->findEigenStates(vs[index], emin, emax, states, &nElec);

	freePsis();
	psis[0] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

	vtls::copyArray(nPts * nElec, states, psis[0]);

	sq_free(states);

	for (int i = 0; i < nElec; i++)
		vtls::normalizeSqrNorm(nPts, &psis[0][i * nPts], dx);

	for (int i = 1; i < 4; i++) {
		psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
		vtls::copyArray(nPts * nElec, psis[0], psis[i]);
	}

	calcWeights();
	if(calcDensity)
		dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
}

void SimulationManager::setPsi(std::complex<double>* npsi) {
	if (!nElec) {
		nElec = 1;
		freePsis();
		for (int i = 0; i < 4; i++) {
			psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
		}
	}
	vtls::copyArray(nPts, npsi, psis[index]);
	vtls::normalizeSqrNorm(nPts, psis[index], dx);
}

int SimulationManager::updatePotential(std::complex<double>* psi, int idx, double* rho) {
	auto strt = std::chrono::high_resolution_clock::now();
	if(calcDensity)
		dens->calcRho(nPts, nElec, dx, weights, psi, rho);
	pot->getV(rho, psi, ts[idx], vs[idx]);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

int SimulationManager::measure(int idx) {
	auto strt = std::chrono::high_resolution_clock::now();
	meas->measure(step[idx], psis[idx], vs[idx], ts[idx]);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

//Run simulation using operator splitting Fourier method (applies potential as linear)
void SimulationManager::runOS_U2TU() {
	updatePotential(psis[prevIndex()], prevIndex(), rhos[prevIndex()]);
	kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nElec);
	iterateIndex();

	int percDone = 0;
	std::future<int> f1;
	auto rM = &SimulationManager::measure;
	while (step[prevPrevIndex()] < numSteps) {
		f1 = std::async(rM, this, prevPrevIndex());

		updatePotential(psis[prevIndex()], prevIndex(), rhos[prevIndex()]);
		kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nElec);

		f1.get();

		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > percDone) {
			if (progCallback != NULL)
				progCallback(percDone);
			percDone++;
		}
	}
	if (progCallback != NULL)
		progCallback(percDone);
}

//Run simulation using operator splitting Fourier method (applies potential as nonlinear, second potential phase is recalculated after propagation phase)
void SimulationManager::runOS_UW2TUW() {
	std::complex<double>* tpsi = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

	updatePotential(psis[prevIndex()], prevIndex(), rhos[prevIndex()]);
	kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nElec);
	ts[prevIndex()] += dt / 2.0;
	updatePotential(tpsi, prevIndex(), rhos[prevIndex()]);
	kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nElec);
	iterateIndex();

	int percDone = 0;
	std::future<int> f1;
	auto rM = &SimulationManager::measure;
	while (step[prevPrevIndex()] < numSteps) {
		f1 = std::async(rM, this, prevPrevIndex());

		updatePotential(psis[prevIndex()], prevIndex(), rhos[prevIndex()]);
		kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nElec);
		ts[prevIndex()] += dt / 2.0;
		updatePotential(tpsi, prevIndex(), rhos[prevIndex()]);
		kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nElec);
		f1.get();

		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > percDone) {
			if (progCallback != NULL)
				progCallback(percDone);
			percDone++;
		}
	}

	sq_free(tpsi);
}

void SimulationManager::iterateIndex() {
	step[nextIndex()] = step[index] + 1;
	ts[nextIndex()] = step[nextIndex()] * dt;

	index++;
	if (index > 3)
		index = 0;
}

int SimulationManager::prevIndex() {
	if (index == 0)
		return 3;
	else
		return index - 1;
}

int SimulationManager::prevPrevIndex() {
	if (index <= 1)
		return 2 + index;
	else
		return index - 2;
}

int SimulationManager::nextIndex() {
	if (index == 3)
		return 0;
	else
		return index + 1;
}

int SimulationManager::getNumPoints() {
	return nPts;
}

double SimulationManager::getDX() {
	return dx;
}

double SimulationManager::getDT() {
	return dt;
}

double SimulationManager::getMaxT() {
	return maxT;
}

int SimulationManager::getNumSteps(){
	return numSteps;
}

std::complex<double> * SimulationManager::getPsi() {
	return psis[index];
}

double * SimulationManager::getRho(){
	return rhos[index];
}

int SimulationManager::getNElec() {
	return nElec;
}

int* SimulationManager::getNElecPtr(){
	return &nElec;
}