#include "SimulationManager.h"

//callback sends progress int 0-100 (can be nullptr for no callback)
SimulationManager::SimulationManager(int nPts, double dx, double dt, double maxT, std::function<int(int)> callback)
	: maxT(maxT), dt(dt), dx(dx), nPts(nPts)
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
	SimulationManager::maxT = maxT; SimulationManager::dt = dt; SimulationManager::dx = dx; SimulationManager::nPts = nPts;
	SimulationManager::mpiRoot = mpiRoot; SimulationManager::mpiUpdateTag = mpiUpdateTag; SimulationManager::mpiJob = mpiJob;
	*/

	progCallback = callback;

	index = 0;
	nelec = 0;
}

SimulationManager::SimulationManager(int nPts, double dx, double dt, double maxT, std::function<int(int)> callback, void* args)
	: args(args)
{
	SimulationManager(nPts, dx, dt, maxT, callback);

}

SimulationManager::~SimulationManager()
{
	meas->terminate();
}

void SimulationManager::addMeasurer(Measurers::Measurer* m) {
	meas->addMeasurer(m);
}

void SimulationManager::addPotential(Potentials::Potential* p) {
	std::cout << "00" << std::endl;
	pot->addPotential(p);
	std::cout << "01" << std::endl;
}

void SimulationManager::addSpatialDamp(double* arr) {
	vtls::seqMulArrays(nPts, arr, spatialDamp);
}

void SimulationManager::finishInitialization() {
	pot->finishAddingPotentials(kin);
	//it->initializeCN();
	//it->initializeCNPA();
}

//SHOULD USE KIN METHOD OR WFCRHOTOOLS
double SimulationManager::getTotalEnergy(std::complex<double> * psi, double * v) {
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


void SimulationManager::findEigenStates(double emin, double emax, double maxT, double rate) {
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

void SimulationManager::setPsi(std::complex<double>* npsi) {
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

int SimulationManager::getVPAR(int idx, int idxPsi) {
	auto strt = std::chrono::high_resolution_clock::now();
	pot->getV(psis[idxPsi], ts[idx], vs[idx], kin);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

int SimulationManager::measPAR(int idx) {
	auto strt = std::chrono::high_resolution_clock::now();
	meas->measureMany(psis[idx], vs[idx], ts[idx], kin, nelec, nPts);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

//Run simulation using operator splitting Fourier method (applies potential as linear)
void SimulationManager::runOS_U2TU(int idx) {
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
	kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nelec);
	iterateIndex();

	int percDone = 0;
	std::future<int> f1;
	auto rM = &SimulationManager::measPAR;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());

		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
		kin_psm->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], spatialDamp, psis[index], nelec);

		f1.get();

		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > percDone) {
			//PROGRESS UPDATE
			//MPI_Ssend(&percDone, 1, MPI_INT, MPI_Root_Proc, MPITag::UpdateSent, MPI_COMM_WORLD);
			if (progCallback != NULL)
				progCallback(percDone);
			percDone++;
		}
	}
	//MPI_Ssend(&percDone, 1, MPI_INT, MPI_Root_Proc, MPITag::UpdateSent, MPI_COMM_WORLD);
	if (progCallback != NULL)
		progCallback(percDone);
	meas->terminate();
}

//Run simulation using operator splitting Fourier method (applies potential as nonlinear, second potential phase is recalculated after propagation phase)
void SimulationManager::runOS_UW2TUW(int idx) {
	tpsi = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec);
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
	kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nelec);
	pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()], kin);
	kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nelec);
	iterateIndex();

	int percDone = 0;
	std::future<int> f1;
	auto rM = &SimulationManager::measPAR;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());

		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()], kin);
		kin_psm->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], spatialDamp, tpsi, nelec);
		pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()], kin);
		kin_psm->stepOS_UW(tpsi, vs[prevIndex()], spatialDamp, psis[index], nelec);
		f1.get();

		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > percDone) {
			//PROGRESS UPDATE
			//MPI_Ssend(&percDone, 1, MPI_INT, MPI_Root_Proc, MPITag::UpdateSent, MPI_COMM_WORLD);
			if (progCallback != NULL)
				progCallback(percDone);
			percDone++;
		}
	}
	meas->terminate();
}

void SimulationManager::iterateIndex() {
	ts[nextIndex()] = ts[index] + dt;
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

std::complex<double> * SimulationManager::getPsi() {
	return psis[index];
}

int SimulationManager::getNElec() {
	return nelec;
}

void* SimulationManager::getargs(){
	return args;
}