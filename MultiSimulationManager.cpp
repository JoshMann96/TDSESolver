#include "stdafx.h"

MultiSimulationManager::MultiSimulationManager(int nPts, double dx, double dt, double maxT)
{
	pot = new Potentials::PotentialManager(nPts);
	meas = new Measurers::MeasurementManager("");
	//psis = new std::complex<double>*[4];
	vs = new double*[4];
	ts = new double[4];
	for (int i = 0; i < 4; i++) {
		//psis[i] = new std::complex<double>[nPts];
		vs[i] = new double[nPts];
		ts[i] = i * dt;
	}
	scratch1 = new std::complex<double>[nPts];
	scratch2 = new std::complex<double>[nPts];
	spatialDamp = new double[nPts];
	std::fill_n(spatialDamp, nPts, 1.0);
	MultiSimulationManager::maxT = maxT; MultiSimulationManager::dt = dt; MultiSimulationManager::dx = dx; MultiSimulationManager::nPts = nPts;

	index = 0;
	nelec = 0;
}

MultiSimulationManager::~MultiSimulationManager()
{
	meas->terminate();
}

void MultiSimulationManager::addMeasurer(Measurers::Measurer* m) {
	meas->addMeasurer(m);
}

void MultiSimulationManager::addPotential(Potentials::Potential* p) {
	pot->addPotential(p);
}

void MultiSimulationManager::addSpatialDamp(double* arr) {
	vtls::seqMulArrays(nPts, arr, spatialDamp);
}

void MultiSimulationManager::finishInitialization() {
	pot->finishAddingPotentials();
	//it->initializeCN();
	//it->initializeCNPA();
}

double MultiSimulationManager::getTotalEnergy(std::complex<double> * psi, double * v) {
	double totE = 0.0;
	for (int i = 0; i < nelec; i++) {
		vtls::secondDerivative(nPts, &psi[i*nPts], scratch1, dx);
		vtls::scaMulArray(nPts, -PhysCon::hbar*PhysCon::hbar / (2.0*PhysCon::me), scratch1);
		vtls::seqMulArrays(nPts, v, &psi[i*nPts], scratch2);
		vtls::addArrays(nPts, scratch2, scratch1);
		for (int j = 0; j < nPts; j++)
			scratch2[j] = std::conj(psi[i*nPts+j]);
		totE += std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));
	}
	return totE;
}


void MultiSimulationManager::findEigenStates(double emin, double emax, double maxT, double rate) {
	pot->getV(0.0, vs[index]);

	kin->findEigenStates(vs[index], emin, emax, &(psis[index]), &nelec);

	vtls::copyArray(nPts * nelec, psis[index], psis[nextIndex()]);
	vtls::copyArray(nPts * nelec, psis[index], psis[prevPrevIndex()]);
	vtls::copyArray(nPts * nelec, psis[index], psis[prevIndex()]);
}

void MultiSimulationManager::setPsi(std::complex<double>* npsi) {
	if (!nelec) {
		psis = new std::complex<double>*[4];
		nelec = 1;
		for (int i = 0; i < 4; i++) {
			psis[i] = new std::complex<double>[nPts];
		}

	}
	vtls::copyArray(nPts, npsi, psis[index]);
	vtls::normalizeSqrNorm(nPts, psis[index], dx);
}

int MultiSimulationManager::getVPAR(int idx, int idxPsi) {
	auto strt = std::chrono::high_resolution_clock::now();
	pot->getV(psis[idxPsi], ts[idx], vs[idx]);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

int MultiSimulationManager::measPAR(int idx) {
	auto strt = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < nelec; i++)
		meas->measure(&psis[idx][i*nPts], vs[idx], ts[idx]);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

//Run simulation using operator splitting Fourier method (applies potential as linear)
void MultiSimulationManager::runOS_U2TU(ProgressTracker* prg, int idx) {
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
	kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nelec);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	//auto t1 = std::chrono::high_resolution_clock::now();
	//auto t2 = std::chrono::high_resolution_clock::now();
	//auto t3 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
	std::future<int> f1;
	auto rM = &MultiSimulationManager::measPAR;
	auto rLog = &ProgressTracker::update;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());

		//t1 = std::chrono::high_resolution_clock::now();
		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
		//t2 = std::chrono::high_resolution_clock::now();
		kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nelec);
		//t3 = std::chrono::high_resolution_clock::now();
		//auto dur2 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		//auto dur3 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
		//std::cout << "Potential: " << dur2.count() << "\nStep: " << dur3.count() << "\n" << std::endl;
		f1.get();
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tf = std::chrono::system_clock::now();
			dur = tf - t0;
			std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count() * (maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
			numPrints++;
		}
	}
	std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count() * (maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
	meas->terminate();
}

//Run simulation using operator splitting Fourier method (applies potential as nonlinear, second potential phase is recalculated after propagation phase)
void MultiSimulationManager::runOS_UW2TUW(ProgressTracker* prg, int idx) {
	tpsi = new std::complex<double>[nPts * nelec];
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
	kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nelec);
	pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()]);
	kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nelec);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	//auto t1 = std::chrono::high_resolution_clock::now();
	//auto t2 = std::chrono::high_resolution_clock::now();
	//auto t3 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
	std::future<int> f1;
	auto rM = &MultiSimulationManager::measPAR;
	auto rLog = &ProgressTracker::update;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());

		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
		kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nelec);
		pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()]);
		kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nelec);

		//t1 = std::chrono::high_resolution_clock::now();
		//t2 = std::chrono::high_resolution_clock::now();
		//t3 = std::chrono::high_resolution_clock::now();
		//auto dur2 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		//auto dur3 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
		//std::cout << "Potential: " << dur2.count() << "\nStep: " << dur3.count() << "\n" << std::endl;
		f1.get();
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tf = std::chrono::system_clock::now();
			dur = tf - t0;
			std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count() * (maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
			numPrints++;
		}
	}
	std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count() * (maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
	meas->terminate();
}

void MultiSimulationManager::iterateIndex() {
	ts[nextIndex()] = ts[index] + dt;
	index++;
	if (index > 3)
		index = 0;
}

int MultiSimulationManager::prevIndex() {
	if (index == 0)
		return 3;
	else
		return index - 1;
}

int MultiSimulationManager::prevPrevIndex() {
	if (index <= 1)
		return 2 + index;
	else
		return index - 2;
}

int MultiSimulationManager::nextIndex() {
	if (index == 3)
		return 0;
	else
		return index + 1;
}

int MultiSimulationManager::getNumPoints() {
	return nPts;
}

double MultiSimulationManager::getDX() {
	return dx;
}

double MultiSimulationManager::getDT() {
	return dt;
}

double MultiSimulationManager::getMaxT() {
	return maxT;
}

std::complex<double> * MultiSimulationManager::getPsi() {
	return psis[index];
}

int MultiSimulationManager::getNElec() {
	return nelec;
}