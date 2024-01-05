#include "stdafx.h"

MultiSimulationManager::MultiSimulationManager(int nPts, double dx, double dt, double maxT, int mpiJob)
	: maxT(maxT), dt(dt), dx(dx), nPts(nPts), mpiJob(mpiJob)
{
	pot = new Potentials::PotentialManager(nPts);
	meas = new Measurers::MeasurementManager("");
	psis = new std::complex<double>*[4];
	vs = new double*[4];
	ts = new double[4];
	for (int i = 0; i < 4; i++) {
		vs[i] = (double*) fftw_malloc(sizeof(double) * nPts);
		ts[i] = i * dt;
	}
	scratch1 = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
	scratch2 = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
	spatialDamp = (double*) fftw_malloc(sizeof(double) * nPts);
	std::fill_n(spatialDamp, nPts, 1.0);
	/*
	MultiSimulationManager::maxT = maxT; MultiSimulationManager::dt = dt; MultiSimulationManager::dx = dx; MultiSimulationManager::nPts = nPts;
	MultiSimulationManager::mpiRoot = mpiRoot; MultiSimulationManager::mpiUpdateTag = mpiUpdateTag; MultiSimulationManager::mpiJob = mpiJob;
	*/

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
	pot->finishAddingPotentials(kin);
	//it->initializeCN();
	//it->initializeCNPA();
}

//SHOULD USE KIN METHOD OR WFCRHOTOOLS
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
	pot->getV(0.0, vs[index], kin);

	std::complex<double>* states = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nPts);

	kin->findEigenStates(vs[index], emin, emax, states, &nelec);

	psis[0] = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec);

	vtls::copyArray(nPts * nelec, states, psis[0]);

	if (states)
		fftw_free(states); states = NULL;

	for (int i = 0; i < nelec; i++)
		vtls::normalizeSqrNorm(nPts, &psis[0][i * nPts], dx);

	for (int i = 1; i < 4; i++) {
		psis[i] = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec);
		vtls::copyArray(nPts * nelec, psis[0], psis[i]);
	}

}

void MultiSimulationManager::setPsi(std::complex<double>* npsi) {
	if (!nelec) {
		psis = new std::complex<double>*[4];
		nelec = 1;
		for (int i = 0; i < 4; i++) {
			psis[i] = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
		}

	}
	vtls::copyArray(nPts, npsi, psis[index]);
	vtls::normalizeSqrNorm(nPts, psis[index], dx);
}

int MultiSimulationManager::getVPAR(int idx, int idxPsi) {
	auto strt = std::chrono::high_resolution_clock::now();
	pot->getV(psis[idxPsi], ts[idx], vs[idx], kin);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

int MultiSimulationManager::measPAR(int idx) {
	auto strt = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < nelec; i++)
		meas->measure(&psis[idx][i*nPts], vs[idx], ts[idx], kin);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

//Run simulation using operator splitting Fourier method (applies potential as linear)
void MultiSimulationManager::runOS_U2TU(int idx) {
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
	kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nelec);
	iterateIndex();
	//auto t0 = std::chrono::system_clock::now();
	//auto tf = std::chrono::system_clock::now();
	//auto t1 = std::chrono::high_resolution_clock::now();
	//auto t2 = std::chrono::high_resolution_clock::now();
	//auto t3 = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> dur;
	//int numPrints = 0;
	int percDone = 0;
	std::future<int> f1;
	auto rM = &MultiSimulationManager::measPAR;
	//auto rLog = &ProgressTracker::update;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());

		//t1 = std::chrono::high_resolution_clock::now();
		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
		//t2 = std::chrono::high_resolution_clock::now();
		kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nelec);
		//t3 = std::chrono::high_resolution_clock::now();
		//auto dur2 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		//auto dur3 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
		//std::cout << "Potential: " << dur2.count() << "\nStep: " << dur3.count() << "\n" << std::endl;
		f1.get();
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > percDone) {
			//tf = std::chrono::system_clock::now();
			//dur = tf - t0;
			//PROGRESS UPDATE
			MPI_Ssend(&percDone, 1, MPI_INT, MPI_Root_Proc, MPITag::UpdateSent, MPI_COMM_WORLD);
			percDone++;
			//std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count() * (maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
			//numPrints++;
		}
	}
	MPI_Ssend(&percDone, 1, MPI_INT, MPI_Root_Proc, MPITag::UpdateSent, MPI_COMM_WORLD);
	//std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count() * (maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
	meas->terminate();
}

//Run simulation using operator splitting Fourier method (applies potential as nonlinear, second potential phase is recalculated after propagation phase)
void MultiSimulationManager::runOS_UW2TUW(int idx) {
	tpsi = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec);
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
	kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nelec);
	pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()], kin);
	kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nelec);
	iterateIndex();
	//auto t0 = std::chrono::system_clock::now();
	//auto tf = std::chrono::system_clock::now();
	//auto t1 = std::chrono::high_resolution_clock::now();
	//auto t2 = std::chrono::high_resolution_clock::now();
	//auto t3 = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> dur;
	//int numPrints = 0;
	int percDone = 0;
	std::future<int> f1;
	auto rM = &MultiSimulationManager::measPAR;
	//auto rLog = &ProgressTracker::update;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());

		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
		kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nelec);
		pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()], kin);
		kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nelec);

		//t1 = std::chrono::high_resolution_clock::now();
		//t2 = std::chrono::high_resolution_clock::now();
		//t3 = std::chrono::high_resolution_clock::now();
		//auto dur2 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		//auto dur3 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
		//std::cout << "Potential: " << dur2.count() << "\nStep: " << dur3.count() << "\n" << std::endl;
		f1.get();
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > percDone) {
			//PROGRESS UPDATE
			MPI_Ssend(&percDone, 1, MPI_INT, MPI_Root_Proc, MPITag::UpdateSent, MPI_COMM_WORLD);
			percDone++;
			//tf = std::chrono::system_clock::now();
			//dur = tf - t0;
			//std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count() * (maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
			//numPrints++;
		}
	}
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