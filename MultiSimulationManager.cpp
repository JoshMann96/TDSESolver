#include "stdafx.h"

MultiSimulationManager::MultiSimulationManager(int nPts, double dx, double dt, double maxT)
{
	pot = new Potentials::PotentialManager(nPts);
	meas = new Measurers::MeasurementManager("");
	reg = new AbsorptiveRegions::AbsorptiveRegionManager();
	it = new TDSEIterator1D(dx, dt, nPts);
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

void MultiSimulationManager::addAbsorptiveRegion(AbsorptiveRegions::AbsorptiveRegion* r) {
	reg->addAbsorptiveRegion(r);
}

void MultiSimulationManager::addTimeRotation(std::complex<double> * arr) {
	it->addTimePhase(arr);
}

void MultiSimulationManager::addSpatialDamp(double* arr) {
	it->addSpatialDamp(arr);
}

void MultiSimulationManager::finishInitialization() {
	pot->finishAddingPotentials();
	it->initializeCN();
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


void MultiSimulationManager::findGroundState(int nSteps, double prec) {
	psis = new std::complex<double>*[4];
	nelec = 1;
	for (int i = 0; i < 4; i++) {
		psis[i] = new std::complex<double>[nPts];
	}

	pot->getV(0.0, vs[index]);
	double curE = getTotalEnergy(psis[index], vs[index]);
	double prevE;
	double rate = 1.0;
	double prog;
	do {
		prevE = curE;
		for (int i = 0; i < nSteps / 2; i++) {
			it->stepCNImT(psis[index], vs[index], psis[nextIndex()], rate);
			vtls::normalizeSqrNorm(nPts, psis[nextIndex()], dx);
			it->stepCNImT(psis[nextIndex()], vs[index], psis[index], rate);
			vtls::normalizeSqrNorm(nPts, psis[index], dx);
		}
		curE = getTotalEnergy(psis[index], vs[index]);
		prog = std::log10(prec / std::abs(curE - prevE));
		rate = -prog / std::log10(nSteps);
	} while (std::abs(curE - prevE) > prec);
	vtls::copyArray(nPts, psis[index], psis[nextIndex()]);
	vtls::copyArray(nPts, psis[index], psis[prevPrevIndex()]);
	vtls::copyArray(nPts, psis[index], psis[prevIndex()]);
}

void MultiSimulationManager::findGroundStatePrintProgress(int nSteps, double prec) {
	psis = new std::complex<double>*[4];
	nelec = 1;
	for (int i = 0; i < 4; i++) {
		psis[i] = new std::complex<double>[nPts];
	}

	pot->getV(0.0, vs[index]);
	double curE = getTotalEnergy(psis[index], vs[index]);
	double prevE;
	double rate = 1.0;
	double prog;
	do {
		prevE = curE;
		for (int i = 0; i < nSteps / 2; i++) {
			it->stepCNImT(psis[index], vs[index], psis[nextIndex()], rate);
			vtls::normalizeSqrNorm(nPts, psis[nextIndex()], dx);
			it->stepCNImT(psis[nextIndex()], vs[index], psis[index], rate);
			vtls::normalizeSqrNorm(nPts, psis[index], dx);
		}
		curE = getTotalEnergy(psis[index], vs[index]);
		prog = std::log10(prec / std::abs(curE - prevE));
		rate = -prog / std::log10(nSteps);
		std::cout << "\rImaginary Time Propagation: 10^" << std::setw(4) << (int)prog << std::flush;
	} while (std::abs(curE - prevE) > prec);
	vtls::copyArray(nPts, psis[index], psis[nextIndex()]);
	vtls::copyArray(nPts, psis[index], psis[prevPrevIndex()]);
	vtls::copyArray(nPts, psis[index], psis[prevIndex()]);
	std::cout << "\rImaginary Time Propagation: DONE          " << std::endl;
}

void MultiSimulationManager::findGroundStateWithWaveFuncDepPot(int nSteps, double prec) {
	psis = new std::complex<double>*[4];
	nelec = 1;
	for (int i = 0; i < 4; i++) {
		psis[i] = new std::complex<double>[nPts];
	}

	pot->getV(psis[index], 0.0, vs[index]);
	double curE = getTotalEnergy(psis[index], vs[index]);
	double prevE;
	double rate = 1.0;
	double prog;
	do {
		prevE = curE;
		for (int i = 0; i < nSteps / 2; i++) {
			it->stepCNImT(psis[index], vs[index], psis[nextIndex()], rate);
			vtls::normalizeSqrNorm(nPts, psis[nextIndex()], dx);
			pot->getV(psis[nextIndex()], 0.0, vs[index]);
			it->stepCNImT(psis[nextIndex()], vs[index], psis[index], rate);
			vtls::normalizeSqrNorm(nPts, psis[index], dx);
			pot->getV(psis[index], 0.0, vs[index]);
		}
		curE = getTotalEnergy(psis[index], vs[index]);
		prog = std::log10(prec / std::abs(curE - prevE));
		rate = -prog / std::log10(nSteps);
		std::cout << "\rImaginary Time Propagation: 10^" << std::setw(4) << (int)prog << std::flush;
	} while (std::abs(curE - prevE) > prec);
	vtls::copyArray(nPts, psis[index], psis[nextIndex()]);
	vtls::copyArray(nPts, psis[index], psis[prevPrevIndex()]);
	vtls::copyArray(nPts, psis[index], psis[prevIndex()]);
	std::cout << "\rImaginary Time Propagation: DONE          " << std::endl;
}

void MultiSimulationManager::findEigenState(double energ, double tolerance) {
	psis = new std::complex<double>*[4];
	nelec = 1;
	for (int i = 0; i < 4; i++) {
		psis[i] = new std::complex<double>[nPts];
	}

	pot->getV(0.0, vs[index]);
	double * d = new double[nPts];
	double * e = new double[nPts];
	int numEigens = 0, nsplit = 0;
	double * eigs = new double[nPts];
	int * iblock = new int[nPts];
	int * isplit = new int[nPts];
	for (int i = 0; i < nPts; i++)
		d[i] = -PhysCon::hbar*PhysCon::hbar / (2.0*PhysCon::me)*(-2.0) / (dx*dx) + vs[index][i];
	for (int i = 0; i < nPts; i++)
		e[i] = -PhysCon::hbar*PhysCon::hbar / (2.0*PhysCon::me) / (dx*dx);

	do {
		LAPACKE_dstebz('V', 'B', nPts, energ - tolerance, energ + tolerance, 0, 0, 0.0, d, e, &numEigens, &nsplit, eigs, iblock, isplit);
		//std::complex<double> * states = new std::complex<double>[numEigens*nPts];
		tolerance *= 2.0;
		if (numEigens == 0)
			std::cout << "Did not find any eigenvalues in the provided range. Doubling tolerance." << std::endl;
	} while (numEigens == 0);
	//std::cout << "Found " << numEigens << " eigenvalue(s)." << std::endl;
	tolerance /= 2.0;
	std::complex<double> * states = new std::complex<double>[numEigens*nPts];
	int * isupz = new int[numEigens * 2];
	int *ifailv = new int[numEigens];
	LAPACKE_zstein(LAPACK_COL_MAJOR, nPts, d, e, numEigens, eigs, iblock, isplit, states, nPts, ifailv);
	std::complex<double> * psi = new std::complex<double>[nPts];// , *psi2 = new std::complex<double>[nPts], *psi3 = new std::complex<double>[nPts];
	int idx = 0;
	double dif = std::abs(eigs[0] - energ);
	for (int i = 1; i < numEigens; i++) {
		if (std::abs(eigs[i] - energ) < dif) {
			dif = std::abs(eigs[i] - energ);
			idx = i;
		}
	}
	vtls::copyArray(nPts, &states[nPts*idx], psi);

	/*
	std::fill_n(psi, nPts, 0.0);
	for (int i = 0; i < numEigens; i++)
		vtls::addArrays(nPts, &states[nPts*i], psi);
		*/

		/*
		vtls::scaAddArray(nPts, -energ, vs[index]);
		for (int i = 0; i < 10000; i++) {
			vtls::normalizeSqrNorm(nPts, psi, dx);
			it->stepCN(psi, vs[index], vs[index], psi2);
			for (int k = 0; k < nPts; k++)
				psi[k] = std::real(psi2[k]);
		}
		vtls::scaAddArray(nPts, energ, vs[index]);
		*/

		/*
		vtls::scaAddArray(nPts, -energ, vs[index]);
		for (int j = 0; j < 100; j++) {
			for (int i = 0; i < 100; i++) {
				vtls::normalizeSqrNorm(nPts, psi, dx);
				it->stepH2Euler(psi, vs[index], psi2);
				it->stepH2Euler(psi2, vs[index], psi);
			}
			std::cout << j << std::endl;
		}
		vtls::scaAddArray(nPts, energ, vs[index]);
		*/

	vtls::normalizeSqrNorm(nPts, psi, dx);
	vtls::copyArray(nPts, psi, psis[nextIndex()]);
	vtls::copyArray(nPts, psi, psis[prevPrevIndex()]);
	vtls::copyArray(nPts, psi, psis[prevIndex()]);
	vtls::copyArray(nPts, psi, psis[index]);

	if (d)
		delete[] d; d = NULL;
	if (e)
		delete[] e; e = NULL;
	if (eigs)
		delete[] eigs; eigs = NULL;
	if (isupz)
		delete[] isupz; isupz = NULL;
	if (states)
		delete[] states; states = NULL;
	if (psi)
		delete[] psi; psi = NULL;
	if (iblock)
		delete[] iblock; iblock = NULL;
	if (isplit)
		delete[] isplit; isplit = NULL;
	if (ifailv)
		delete[] ifailv; ifailv = NULL;

}

void MultiSimulationManager::findMetallicInitialState_PSM(double fermie, double w, double maxT, double rate) {
	pot->getV(0.0, vs[index]);
	if (nPts > 46340) {
		std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (findMetallicInitialState_PSW)" << std::endl;
		throw NULL;
	}
	double* ap = new double[nPts * (nPts + 1) / 2];
	double* states = new double[nPts * nPts];
	double* eigs = new double[nPts];
	int* ifail = new int[nPts];
	int numEigens = 0;
	double tmul = PhysCon::hbar * PhysCon::hbar * 2 * PhysCon::pi * PhysCon::pi / (PhysCon::me * dx * dx);
	for (int j = 0; j < nPts; j++) {
		for (int i = 0; i <= j; i++) {
			if (i == j)
				ap[i + (long)j * (j + 1) / 2] = tmul * (1.0 / 12.0 + 1.0 / (6.0 * nPts * nPts)) + vs[index][j];
			else
				ap[i + (long)j * (j + 1) / 2] = tmul * std::pow(-1, i - j) / (2.0 * std::pow(nPts * std::sin(PhysCon::pi * (i - j) / nPts), 2));
		}
	}
	//ap is stored upper triangular matrix, eigs needs to be initialized to have size npts, vecs npts^2... 
	LAPACKE_dspevx(LAPACK_COL_MAJOR, 'V', 'V', 'U', nPts, ap, -fermie - w, -w, 0, 0, 2 * LAPACKE_dlamch('S'), &numEigens, eigs, states, nPts, ifail);

	psis = new std::complex<double>*[4];
	nelec = numEigens;
	for (int i = 0; i < 4; i++) {
		psis[i] = new std::complex<double>[nPts * nelec];
	}
	vtls::copyArray(nPts * nelec, states, psis[index]);
	for (int i = 0; i < numEigens; i++)
		vtls::normalizeSqrNorm(nPts, &psis[index][i * nPts], dx);

	vtls::copyArray(nPts * nelec, psis[index], psis[nextIndex()]);
	vtls::copyArray(nPts * nelec, psis[index], psis[prevPrevIndex()]);
	vtls::copyArray(nPts * nelec, psis[index], psis[prevIndex()]);

	if (eigs)
		delete[] eigs; eigs = NULL;
	if (ap)
		delete[] ap; ap = NULL;
	if (states)
		delete[] states; states = NULL;
	if (ifail)
		delete[] ifail; ifail = NULL;
	
	std::cout << "Created " << nelec << " electrons." << std::endl;
}

void MultiSimulationManager::findMetallicInitialState(double fermie, double w, double maxT, double rate) {
	//https://software.intel.com/en-us/mkl-developer-reference-c-symmetric-eigenvalue-problems-lapack-computational-routines#E20540C0-34FD-4CD2-B803-B6096460B4F2
	pot->getV(0.0, vs[index]);
	double * d = new double[nPts];
	double * e = new double[nPts];
	int numEigens = 0, nsplit = 0;
	double * eigs = new double[nPts];
	int * iblock = new int[nPts];
	int * isplit = new int[nPts];
	for (int i = 0; i < nPts; i++)
		d[i] = -PhysCon::hbar*PhysCon::hbar / (2.0*PhysCon::me)*(-2.0) / (dx*dx) + vs[index][i];
	for (int i = 0; i < nPts; i++)
		e[i] = -PhysCon::hbar*PhysCon::hbar / (2.0*PhysCon::me) / (dx*dx);

	LAPACKE_dstebz('V', 'B', nPts, -fermie - w, -w, 0, 0, 0.0, d, e, &numEigens, &nsplit, eigs, iblock, isplit);
	std::complex<double> * states = new std::complex<double>[numEigens*nPts];
	int * isupz = new int[numEigens * 2];
	int *ifailv = new int[numEigens];
	LAPACKE_zstein(LAPACK_COL_MAJOR, nPts, d, e, numEigens, eigs, iblock, isplit, states, nPts, ifailv);

	psis = new std::complex<double>*[4];
	nelec = numEigens;
	for (int i = 0; i < 4; i++) {
		psis[i] = new std::complex<double>[nPts*nelec];
	}
	vtls::copyArray(nPts*nelec, states, psis[index]);
	for (int i = 0; i < numEigens; i++)
		vtls::normalizeSqrNorm(nPts, &psis[index][i*nPts], dx);

	vtls::copyArray(nPts*nelec, psis[index], psis[nextIndex()]);
	vtls::copyArray(nPts*nelec, psis[index], psis[prevPrevIndex()]);
	vtls::copyArray(nPts*nelec, psis[index], psis[prevIndex()]);

	if (d)
		delete[] d; d = NULL;
	if (e)
		delete[] e; e = NULL;
	if (eigs)
		delete[] eigs; eigs = NULL;
	if (isupz)
		delete[] isupz; isupz = NULL;
	if (states)
		delete[] states; states = NULL;
	if (iblock)
		delete[] iblock; iblock = NULL;
	if (isplit)
		delete[] isplit; isplit = NULL;
	if (ifailv)
		delete[] ifailv; ifailv = NULL;

	std::cout << "Created " << nelec << " electrons." << std::endl;
}

void MultiSimulationManager::findMetallicInitialState_HOD(double fermie, double w, double maxT, double rate, int order) {
	if (nPts > 46340) {
		std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (findMetallicInitialState_HOD)" << std::endl;
		throw NULL;
	}
	pot->getV(0.0, vs[index]);
	//calculate derivative coefficients
	double* coeffMat = new double[(2 * order + 1) * (2 * order + 1)];
	for (int j = -order; j <= order; j++)
		coeffMat[(j + order) * (2 * order + 1)] = 1;
	for (int i = 1; i <= 2 * order; i++) {
		for (int j = -order; j <= order; j++) {
			coeffMat[i + (j + order) * (2 * order + 1)] = std::pow(j, i);
		}
	}
	double* b = new double[2 * order + 1];
	std::fill_n(b, 2 * order + 1, 0);
	b[2] = 2;
	int* ipiv = new int[2 * order + 1];
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, 2 * order + 1, 2 * order + 1, coeffMat, 2 * order + 1, ipiv);
	LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', 2 * order + 1, 1, coeffMat, 2 * order + 1, ipiv, b, 2 * order + 1);


	//build hamiltonian (upper only)
	double* ham = new double[nPts * (order + 1)];

	//i=order means diagonal, i=0 farthest from diag
	double dconst = -PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me) / (dx * dx);
	for (int i = 0; i <= order; i++) {
		for (int j = 0; j < nPts; j++) {
			ham[i + j * (order + 1)] = dconst * b[i];
			if (i == order)
				ham[i + j * (order + 1)] += vs[index][i];
		}
	}

	//get eigenvalues & vectors
	int numeigs;
	double* eigs = new double[nPts];
	int* ifail = new int[nPts];
	double* psires = new double[nPts*nPts];
	double* q = new double[nPts*nPts];
	LAPACKE_dsbevx(LAPACK_COL_MAJOR, 'V', 'V', 'U', nPts, order, ham, order + 1, q, nPts, -fermie - w, -w, 0, 0, 2 * LAPACKE_dlamch('S'), &numeigs, eigs, psires, nPts, ifail);

	psis = new std::complex<double>*[4];
	nelec = numeigs;
	for (int i = 0; i < 4; i++) {
		psis[i] = new std::complex<double>[nPts * nelec];
	}
	vtls::copyArray(nPts * nelec, psires, psis[index]);
	for (int i = 0; i < numeigs; i++)
		vtls::normalizeSqrNorm(nPts, &psis[index][i * nPts], dx);

	vtls::copyArray(nPts * nelec, psis[index], psis[nextIndex()]);
	vtls::copyArray(nPts * nelec, psis[index], psis[prevPrevIndex()]);
	vtls::copyArray(nPts * nelec, psis[index], psis[prevIndex()]);

	if (coeffMat)
		delete[] coeffMat; coeffMat = NULL;
	if (b)
		delete[] b; b = NULL;
	if (eigs)
		delete[] eigs; eigs = NULL;
	if (ipiv)
		delete[] ipiv; ipiv = NULL;
	if (ham)
		delete[] ham; ham = NULL;
	if (ifail)
		delete[] ifail; ifail = NULL;
	if (psires)
		delete[] psires; psires = NULL;
	if (q)
		delete[] q; q = NULL;

	std::cout << "Created " << nelec << " electrons." << std::endl;
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

void MultiSimulationManager::run() {
	pot->getV(psis[prevIndex()], ts[index], vs[index]);
	for (int i = 0; i < nelec; i++)
		meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
	iterateIndex();
	while (ts[prevPrevIndex()] < maxT) {
		pot->getV(psis[prevIndex()], ts[index], vs[index]);
		it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		for(int i = 0; i < nelec; i++)
			meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
		iterateIndex();
	}
	meas->terminate();
}

void MultiSimulationManager::runLogProgress(ProgressTracker * prg, int idx) {
	int numPrints = 0;
	pot->getV(psis[prevIndex()], ts[index], vs[index]);
	for (int i = 0; i < nelec; i++)
		meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;
	auto rLog = &ProgressTracker::update;
	while (ts[prevPrevIndex()] < maxT) {
		pot->getV(psis[prevIndex()], ts[index], vs[index]);
		it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		for (int i = 0; i < nelec; i++)
			meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tf = std::chrono::system_clock::now();
			dur = tf - t0;
			std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
			numPrints++;
		}
	}
	meas->terminate();
}

void MultiSimulationManager::runPrintProgress() {
	pot->getV(psis[prevIndex()], ts[index], vs[index]);
	for (int i = 0; i < nelec; i++)
		meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
	while (ts[prevPrevIndex()] <= maxT) {
		pot->getV(psis[prevIndex()], ts[index], vs[index]);
		it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		for (int i = 0; i < nelec; i++)
			meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tf = std::chrono::system_clock::now();
			dur = tf - t0;
			std::cout << "\r" << std::setw(3) << (int)(ts[prevPrevIndex()] / maxT * 100.0) << "%"
				<< ", Estimated Time Remaining: " << std::setw(7) << (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]) << " s" << std::flush;
			numPrints++;
		}
	}
	std::cout << std::endl;
	tf = std::chrono::system_clock::now();
	dur = tf - t0;
	std::cout << "Simulation Time: " << dur.count() << " s" << std::endl;
	meas->terminate();
}


int MultiSimulationManager::getVPAR(int idx, int idxPsi) {
	auto strt = std::chrono::high_resolution_clock::now();
	pot->getV(psis[idxPsi], ts[idx], vs[idx]);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

int MultiSimulationManager::stepItPAR(int idx0, int idx1) {
	auto strt = std::chrono::high_resolution_clock::now();
	it->multiStepCN(psis[idx0], vs[idx0], vs[idx1], psis[idx1], nelec);
	//it->stepOS_U2TU(psis[idx0], vs[idx0], psis[idx1], nelec);
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


void MultiSimulationManager::runParallel() {
	for (int k = 0; k < 2; k++) {
		pot->getV(psis[index], ts[index], vs[index]);
		for (int i = 0; i < nelec; i++)
			meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
		iterateIndex();
	}
	pot->getV(psis[index], ts[index], vs[index]);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
	std::future<int> f1, f2, f3;
	auto rV = &MultiSimulationManager::getVPAR;
	auto rS = &MultiSimulationManager::stepItPAR;
	auto rM = &MultiSimulationManager::measPAR;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rV, this, index, prevPrevIndex());
		f2 = std::async(rS, this, prevPrevIndex(), prevIndex());
		f3 = std::async(rM, this, prevPrevIndex());
		f1.get(); f2.get(); f3.get();
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tf = std::chrono::system_clock::now();
			dur = tf - t0;
			std::cout << "\r" << std::setw(3) << (int)(ts[prevPrevIndex()] / maxT * 100.0) << "%"
				<< ", Estimated Time Remaining: " << std::setw(7) << (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]) << " s" << std::flush;
			numPrints++;
		}
	}
	std::cout << std::endl;
	tf = std::chrono::system_clock::now();
	dur = tf - t0;
	std::cout << "Simulation Time: " << dur.count() << " s" << std::endl;
	meas->terminate();
}

void MultiSimulationManager::runParallelLogProgress(ProgressTracker * prg, int idx) {
	for (int k = 0; k < 2; k++) {
		pot->getV(psis[index], ts[index], vs[index]);
		for (int i = 0; i < nelec; i++)
			meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
		iterateIndex();
	}
	pot->getV(psis[index], ts[index], vs[index]);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
	std::future<int> f1, f2, f3, f4;
	auto rV = &MultiSimulationManager::getVPAR;
	auto rS = &MultiSimulationManager::stepItPAR;
	auto rM = &MultiSimulationManager::measPAR;
	auto rLog = &ProgressTracker::update;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rV, this, index, prevPrevIndex());
		f2 = std::async(rS, this, prevPrevIndex(), prevIndex());
		f3 = std::async(rM, this, prevPrevIndex());
		f1.get(); f2.get(); f3.get();
		//std::cout << f1.get() << "/" << f2.get() << "/" << f3.get() << std::endl;
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tf = std::chrono::system_clock::now();
			dur = tf - t0;
			std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
			numPrints++;
		}
	}
	std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
	meas->terminate();
}

void MultiSimulationManager::runSCFLogProgress(ProgressTracker * prg, int idx, double thresh) {
	pot->getV(psis[index], ts[index], vs[index]);
	iterateIndex();
	pot->getV(psis[index], ts[index], vs[index]);
	it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
	std::future<int> f1, f2, f3, f4;
	auto rM = &MultiSimulationManager::measPAR;
	auto rLog = &ProgressTracker::update;
	double dev = 0, odev = 0;
	std::complex<double> * comp = new std::complex<double>[nPts*nelec];
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());
		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
		//vtls::copyArray(nPts, vs[index], vs[nextIndex()]);
		pot->getV(psis[prevIndex()], ts[index], vs[index]);
		it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		//it->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		vtls::copyArray(nPts*nelec, psis[index], comp);
		dev = __DBL_MAX__;
		do {
			pot->getV(psis[index], ts[index], vs[index]);
			it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
			//it->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
			odev = dev;
			dev = 0;
			for (int i = 0; i < nPts*nelec; i++)
				dev += std::pow(std::abs(psis[index][i] - comp[i]), 2);
			dev *= dx;
			vtls::copyArray(nPts*nelec, psis[index], comp);
			std::cout << dev / thresh << std::endl;
			if (dev>odev) {
				std::cout << "Simulation struggling to converge." << std::endl;
				vtlsPrnt::printGraph(nPts, vs[index]);
				throw("Not Converging");
			}
		} while (dev > thresh);
		f1.get();
		iterateIndex();
		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tf = std::chrono::system_clock::now();
			dur = tf - t0;
			std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
			numPrints++;
		}
	}
	std::async(std::launch::async, rLog, prg, idx, ts[prevPrevIndex()] / maxT, (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]));
	meas->terminate();
	delete[] comp;
}

void MultiSimulationManager::runPertLogProgress(ProgressTracker* prg, int idx, int nit) {
	pot->getV(psis[index], ts[index], vs[index]);
	iterateIndex();
	pot->getV(psis[index], ts[index], vs[index]);
	it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
	std::future<int> f1;
	auto rM = &MultiSimulationManager::measPAR;
	auto rLog = &ProgressTracker::update;
	while (ts[prevPrevIndex()] <= maxT) {
		f1 = std::async(rM, this, prevPrevIndex());
		pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
		//vtls::copyArray(nPts, vs[index], vs[nextIndex()]);
		pot->getV(psis[prevIndex()], ts[index], vs[index]);
		it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		//it->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		for(int j = 0; j < nit - 1; j++){
			pot->getV(psis[index], ts[index], vs[index]);
			it->multiStepCN(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
			//it->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], vs[index], psis[index], nelec);
		}
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

//Run simulation using operator splitting Fourier method (applies potential as linear)
void MultiSimulationManager::runOS_U2TULogProgress(ProgressTracker* prg, int idx) {
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
	it->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], psis[index], nelec);
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
		it->stepOS_U2TU(psis[prevIndex()], vs[prevIndex()], psis[index], nelec);
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
void MultiSimulationManager::runOS_UW2TUWLogProgress(ProgressTracker* prg, int idx) {
	tpsi = new std::complex<double>[nPts * nelec];
	iterateIndex();
	pot->getV(psis[prevIndex()], ts[prevIndex()], vs[prevIndex()]);
	it->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], tpsi, nelec);
	pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()]);
	it->stepOS_UW(tpsi, vs[prevIndex()], psis[index], nelec);
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
		it->stepOS_UW2T(psis[prevIndex()], vs[prevIndex()], tpsi, nelec);
		pot->getV(tpsi, ts[prevIndex()], vs[prevIndex()]);
		it->stepOS_UW(tpsi, vs[prevIndex()], psis[index], nelec);

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

void MultiSimulationManager::runParallelOMPOld() {
	//omp_set_num_threads(3);
	pot->getV(ts[index], vs[index]);
	meas->measure(psis[index], vs[index], ts[index]);
	iterateIndex();
	pot->getV(ts[index], vs[index]);
	iterateIndex();
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;
	int numPrints = 0;
#pragma omp parallel
	{
		while (ts[prevPrevIndex()] <= maxT) {
#pragma omp sections
			{
#pragma omp section
				{
					pot->getV(ts[index], vs[index]);
				}
#pragma omp section
				{
					it->stepCN(psis[prevPrevIndex()], vs[prevPrevIndex()], vs[prevIndex()], psis[prevIndex()]);
				}
#pragma omp section
				{
					meas->measure(psis[prevPrevIndex()], vs[prevPrevIndex()], ts[prevPrevIndex()]);
				}
			}
#pragma omp single
			{
				iterateIndex();
				if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
					tf = std::chrono::system_clock::now();
					dur = tf - t0;
					std::cout << "\r" << std::setw(3) << (int)(ts[prevPrevIndex()] / maxT * 100.0) << "%"
						<< ", Estimated Time Remaining: " << std::setw(7) << (int)(dur.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]) << " s" << std::flush;
					numPrints++;
				}
			}
		}
	}
	std::cout << std::endl;
	tf = std::chrono::system_clock::now();
	dur = tf - t0;
	std::cout << "Simulation Time: " << dur.count() << " s" << std::endl;
	meas->terminate();
}

void MultiSimulationManager::runTimed() {
	double tpot = 0.0;
	double tmeas = 0.0;
	double tstep = 0.0;
	auto t0 = std::chrono::system_clock::now();
	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur;

	auto t0o = std::chrono::system_clock::now();
	auto tfo = std::chrono::system_clock::now();
	std::chrono::duration<double> duro;
	int numPrints = 0;

	pot->getV(ts[index], vs[index]);
	for (int i = 0; i < nelec; i++)
		meas->measure(&psis[index][i*nPts], vs[index], ts[index]);
	iterateIndex();
	pot->getV(ts[index], vs[index]);
	iterateIndex();

	while (ts[prevPrevIndex()] < maxT) {
		t0 = std::chrono::system_clock::now();
		pot->getV(ts[index], vs[index]);
		tf = std::chrono::system_clock::now();
		dur = tf - t0;
		tpot += dur.count();

		t0 = std::chrono::system_clock::now();
		it->multiStepCN(psis[prevPrevIndex()], vs[prevPrevIndex()], vs[prevIndex()], psis[prevIndex()], nelec);
		tf = std::chrono::system_clock::now();
		dur = tf - t0;
		tstep += dur.count();

		t0 = std::chrono::system_clock::now();
		for (int i = 0; i < nelec; i++)
			meas->measure(&psis[prevPrevIndex()][i*nPts], vs[prevPrevIndex()], ts[prevPrevIndex()]);
		tf = std::chrono::system_clock::now();
		dur = tf - t0;
		tmeas += dur.count();

		iterateIndex();

		if (ts[prevPrevIndex()] / maxT * 100.0 > numPrints) {
			tfo = std::chrono::system_clock::now();
			duro = tfo - t0o;
			std::cout << "\r" << std::setw(3) << (int)(ts[prevPrevIndex()] / maxT * 100.0) << "%"
				<< ", Estimated Time Remaining: " << std::setw(7) << (int)(duro.count()*(maxT - ts[prevPrevIndex()]) / ts[prevPrevIndex()]) << " s" << std::flush;
			numPrints++;
		}
	}
	meas->terminate();
	std::cout << std::endl;
	std::cout << "Potential Function Total Time: " << tpot << " s" << std::endl;
	std::cout << "Measurement Total Time: " << tmeas << " s" << std::endl;
	std::cout << "Crank-Nicholson Total Time: " << tstep << " s" << std::endl;
	std::cout << "Overall Total Time: " << (tpot + tmeas + tstep) << " s" << std::endl;
	std::cin.ignore();
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