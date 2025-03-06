#include "SimulationManager.h"
#include "PhysCon.h"
#include "MathTools.h"

//callback sends progress int 0-100 (can be nullptr for no callback)
SimulationManager::SimulationManager(int nPts, double dx, double dt, std::function<void(int)> callback)
	: dx(dx), nPts(nPts), dt(dt), progTracker(callback), nElec(0)
{
	index = cyclic_int(0, HISTORY_LENGTH);

	pot = new Potentials::PotentialManager(nPts);
	meas = new Measurers::MeasurementManager("");
	psis = (std::complex<double>**) sq_malloc(sizeof(std::complex<double>*)*HISTORY_LENGTH);
	for(int i = 0; i < HISTORY_LENGTH; i++)
		psis[i] = nullptr;

	vs = (double**) sq_malloc(sizeof(double*)*HISTORY_LENGTH);
	rhos = (double**) sq_malloc(sizeof(double*)*HISTORY_LENGTH);
	ts = (double*) sq_malloc(sizeof(double)*HISTORY_LENGTH);
	for (int i = 0; i < HISTORY_LENGTH; i++){
		vs[i] = (double*) sq_malloc(sizeof(double) * nPts);
		rhos[i] = (double*) sq_malloc(sizeof(double) * nPts);
	}
	std::fill_n(ts, HISTORY_LENGTH, 0.0);

	step = (int*) sq_malloc(sizeof(int) * HISTORY_LENGTH);
	std::fill_n(step, HISTORY_LENGTH, 0);

	scratch1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
	scratch2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
	spatialDamp = (double*) sq_malloc(sizeof(double) * nPts);
	std::fill_n(spatialDamp, nPts, 1.0);
	/*
	SimulationManager::maxT = maxT; SimulationManager::dt = dt; SimulationManager::dx = dx; SimulationManager::nPts = nPts;
	SimulationManager::mpiRoot = mpiRoot; SimulationManager::mpiUpdateTag = mpiUpdateTag; SimulationManager::mpiJob = mpiJob;
	*/
}

SimulationManager::~SimulationManager()
{
	for(int i = 0; i < HISTORY_LENGTH; i++){
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
	if (m->needsDensity())
		calcDensity = 1;
}

void SimulationManager::addPotential(Potentials::Potential* p) {
	pot->addPotential(p);
	if(p->getComplexity() == Potentials::PotentialComplexity::WAVEFUNCTION_DEPENDENT)
		calcDensity = 1;
}

void SimulationManager::addSpatialDamp(double* arr) {
	vtls::seqMulArrays(nPts, arr, spatialDamp);
}

void SimulationManager::calcEnergies(int curStep, double* energies) {
		for(int i = 0; i < HISTORY_LENGTH; i++){
			if(curStep == step[i]){ //look for the present step's index
				double* rho = (double*) sq_malloc(sizeof(double)*nPts);
				for(int j = 0; j < nElec; j++){
					vtls::normSqr(nPts, &psis[i][j*nPts], rho);
					energies[j] = vtlsInt::rSumMul(nPts, rho, vs[i], dx)/vtlsInt::rSum(nPts, rho,dx) + kin->evaluateKineticEnergy(&psis[i][j*nPts]);
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

	if (wght == nullptr)
		throw std::runtime_error("SimulationManager::calcEnergies: No weight function set!");

	if(weights)
		sq_free(weights); weights = nullptr;
	weights = (double*) sq_malloc(sizeof(double)*nElec);
	double* energies = (double*) sq_malloc(sizeof(double)*nElec);

	calcEnergies(step[index], energies);

	wght->calcWeights(nElec, energies, weights, normScheme);

	sq_free(energies);
}

void SimulationManager::findEigenStates(double emin, double emax) {
	normScheme = WfcToRho::NormalizationScheme::NORMALIZED;
	
	pot->getVBare(0.0, vs[index]);

	std::complex<double>* states;

	kin->findEigenStates(vs[index], emin, emax, &states, &nElec);

	freePsis();
	psis[0] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

	vtls::copyArray(nPts * nElec, states, psis[0]);

	sq_free(states);

	for (int i = 0; i < nElec; i++)
		vtls::normalizeSqrNorm(nPts, &psis[0][i * nPts], dx);

	for (int i = 1; i < HISTORY_LENGTH; i++) {
		psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
		vtls::copyArray(nPts * nElec, psis[0], psis[i]);
	}

	calcWeights();
	if(calcDensity)
		dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
}

// DEPRECATED
void SimulationManager::findInhomogeneousSteadyStates_OBSOLETE(double threshold, int nElec, double* kl, double* kr, bool verbose){
	this->nElec = nElec;

	freePsis();
	for(int i = 0; i < HISTORY_LENGTH; i++){
		psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
		std::fill_n(psis[i], nPts*nElec, 0.0);
		pot->getVBare(0.0, vs[i]);
	}

	double err;
	double* temp1 = (double*) sq_malloc(sizeof(double) * nPts * nElec);
	double* temp2 = (double*) sq_malloc(sizeof(double) * nPts * nElec);
	bool converged = false;
	int i = 0;
	while(!converged){
		//kin_fdm->projectHistory(psis[prevIndex()], kl, kr, vs[index], nElec); // THIS MUST BE FIXED TO USE THIS FUNCTION AGAIN
		kin_fdm->step(psis[index - 1], vs[index], spatialDamp, psis[index], nElec);

		vtls::normSqr(nPts*nElec, psis[index], temp1);
		vtls::normSqr(nPts*nElec, psis[index - 1], temp2);
		vtls::scaMulAddArrays(nPts*nElec, -1.0, temp1, temp2); // temp2 = old - new
		vtls::abs(nPts*nElec, temp2, temp2); // temp2 = |old - new|
		// ? (sum of |old - new|) / (sum of |new|) < threshold
		err = vtlsInt::rSum(nPts*nElec, temp2, 1.0) / vtlsInt::rSum(nPts*nElec, temp1, 1.0);
		if(verbose && i % 10 == 0)
			std::cout << "Error: " << err << std::endl;
		converged = vtlsInt::rSum(nPts*nElec, temp2, 1.0) / vtlsInt::rSum(nPts*nElec, temp1, 1.0) < threshold;

		vtls::copyArray(nPts*nElec, psis[index], psis[index - 1]);
		i++;
	}

	for(int i = 0; i < HISTORY_LENGTH; i++)
		vtls::copyArray(nPts*nElec, psis[index], psis[i]);

	sq_free(temp1);
	sq_free(temp2);
	calcWeights();
	if(calcDensity)
		dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
	for(int i = 0; i < HISTORY_LENGTH; i++)
		vtls::copyArray(nPts, rhos[index], rhos[i]);
}

void SimulationManager::findInhomogeneousEigenStates(int nElec, double* energies){
	this->nElec = nElec;

	freePsis();
	for(int i = 0; i < HISTORY_LENGTH; i++){
		psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
		std::fill_n(psis[i], nPts*nElec, 0.0);
		pot->getVBare(0.0, vs[i]);
	}

	kin_fdm->findInhomogeneousEigenStates(vs[index], energies, psis[index], nElec);
	for (int i = 1; i < HISTORY_LENGTH; i++) 
		vtls::copyArray(nPts * nElec, psis[index], psis[i]);

	calcWeights();
	if(calcDensity)
		dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);

	normScheme = WfcToRho::NormalizationScheme::UNNORMALIZED;
}

void SimulationManager::setPsi(std::complex<double>* npsi, WfcToRho::NormalizationScheme norm) {
	normScheme = norm;

	if (!nElec) {
		nElec = 1;
		freePsis();
		for (int i = 0; i < HISTORY_LENGTH; i++) {
			psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
		}
	}
	vtls::copyArray(nPts, npsi, psis[index]);

	if(normScheme == WfcToRho::NormalizationScheme::NORMALIZED)
		vtls::normalizeSqrNorm(nPts, psis[index], dx);
}

int SimulationManager::calculatePotential(double* rho, std::complex<double>* psi, double t, double* v){
	auto strt = std::chrono::high_resolution_clock::now();
	if(calcDensity)
		dens->calcRho(nPts, nElec, dx, weights, psi, rho);
	pot->getV(rho, psi, t, v);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

int SimulationManager::updatePotential(int idx) {return calculatePotential(rhos[idx], psis[idx], ts[idx], vs[idx]);}

int SimulationManager::measure(int idx) {
	auto strt = std::chrono::high_resolution_clock::now();
	meas->measure(step[idx], psis[idx], vs[idx], ts[idx]);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

//Run simulation using operator splitting Fourier method (applies potential as linear)
void SimulationManager::runOS_U2TU(int nSteps) {
	auto rMeasure = &SimulationManager::measure;
	auto rUpdatePotential = &SimulationManager::updatePotential;
	std::future<int> fM, fUP;

	// initialize progress tracker
	progTracker.reset(nSteps);

	bool asyncCalc = canAsyncCalcPot();
	if(!asyncCalc)
		std::cout << "Warning: Potential is not wavefunction independent! It is recommended to use runOS_UW2TUW to more accurately account for the nonlinearity." << std::endl;

	for(int i = 0; i < nSteps; i++){
		if(asyncCalc){
			// evaluate potential n+1
			if (i == 0)
				updatePotential(index);
			else
				fUP.get();
			fUP = std::async(rUpdatePotential, this, index + 1);
		}
		else{
			// evaluate potential n
			updatePotential(index);
		}

		// evalute n->n+1
		kin_psm->stepOS_U2TU(psis[index], vs[index], spatialDamp, psis[index + 1], nElec);
		
		// measure step n while n+1->n+2 begins
		if(i != 0)
			fM.get();
		fM = std::async(rMeasure, this, index);
		
		progTracker.update(i);

		iterateIndex();
	}
	
	// collect remaining futures
	if(asyncCalc)
		fUP.get();
	fM.get();

	progTracker.update(nSteps);
}

//Run simulation using operator splitting Fourier method (applies potential as nonlinear, second potential phase is recalculated after propagation phase)
void SimulationManager::runOS_UW2TUW(int nSteps) {
	// variables for the midpoint of step
	std::complex<double>* tpsi = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
	double* trho = (double*) sq_malloc(sizeof(double) * nPts);
	double* tv = (double*) sq_malloc(sizeof(double) * nPts);

	auto rMeasure = &SimulationManager::measure;
	auto rUpdatePotential = &SimulationManager::updatePotential;
	std::future<int> fM, fUP;

	// initialize progress tracker
	progTracker.reset(nSteps);

	for(int i = 0; i < nSteps; i++){
		// step n->n+1/2
		updatePotential(index);
		kin_psm->stepOS_UW2T(psis[index], vs[index], spatialDamp, tpsi, nElec);

		// measure step n while n+1/2->n+1 begins
		if(i != 0)
			fM.get();
		fM = std::async(rMeasure, this, index);

		// step n+1/2->n+1
		calculatePotential(trho, tpsi, ts[index] + dt / 2.0, tv);
		kin_psm->stepOS_UW(tpsi, tv, spatialDamp, psis[index + 1], nElec);
		
		progTracker.update(i);

		iterateIndex();
	}
	
	fM.get();
	progTracker.update(nSteps);

	sq_free(tpsi);
	sq_free(trho);
	sq_free(tv);
}

void SimulationManager::runFD_L(int nSteps){
	double* tpot = (double*) sq_malloc(sizeof(double) * nPts);
	auto rMeasure = &SimulationManager::measure;
	auto rUpdatePotential = &SimulationManager::updatePotential;
	std::future<int> fM, fUP;

	// initialize progress tracker
	progTracker.reset(nSteps);

	bool asyncCalc = canAsyncCalcPot();
	if(!asyncCalc)
		std::cout << "Warning: Potential is not wavefunction independent! It is recommended to use runFD_NL to more accurately account for the nonlinearity." << std::endl;

	for(int i = 0; i < nSteps; i++){
		if(asyncCalc){
			// evaluate potential n+1
			if (i == 0){ // evaluate potential n, n+1 if needed
				updatePotential(index);
				updatePotential(index + 1);
			}
			else
				fUP.get();
			fUP = std::async(rUpdatePotential, this, index + 2);
		}
		else{
			// evaluate potential n
			if(i == 0)
				updatePotential(index);
			updatePotential(index + 1);
		}
		// potentials n, n+1 are now calculated
		// calculate averaged potential
		vtls::addArrays(nPts, vs[index], vs[index + 1], tpot);
		vtls::scaMulArray(nPts, 0.5, tpot);

		// evalute n->n+1
		kin_fdm->step(psis[index], tpot, spatialDamp, psis[index+1], nElec);
		
		// measure step n while n+1->n+2 begins
		if(i != 0)
			fM.get();
		fM = std::async(rMeasure, this, index);
		
		progTracker.update(i);

		iterateIndex();
	}
}

void SimulationManager::iterateIndex() {
	step[index+1] = step[index] + 1;
	ts[index+1] = step[index+1] * dt;

	index++;
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

std::complex<double>* SimulationManager::getPsi() {
	return psis[index];
}

double* SimulationManager::getRho(){
	return rhos[index];
}

int SimulationManager::getNElec() {
	return nElec;
}

int* SimulationManager::getNElecPtr(){
	return &nElec;
}

int SimulationManager::findElectricalSurfaceCentroidRule(int minPos, int maxPos){
	/*
	* Calculate the electrical centroid of the electron density using first-order perturbation theory.
	*/
	std::complex<double>* mat = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nElec*(nElec-1)/2);
	std::complex<double>* xpsi = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	auto matIndex = [this](int i, int j){return i*nElec+j-((i+1)*(i+2))/2;}; //helper function for packing/unpacking matrix

	std::complex<double>* psi = psis[index];

	double* energies = (double*) sq_malloc(sizeof(double)*nElec);
	calcEnergies(step[index], energies);

	double* idxs = (double*) sq_malloc(sizeof(double)*nPts);

	// calculate matrix elements
	for(int i = 0; i < nElec-1; i++){
		for(int k = 0; k < nPts; k++)
			xpsi[k] = ((double)k) * psi[i*nPts+k];
		for(int j = i+1; j < nElec; j++)
			mat[matIndex(i,j)] = vtlsInt::rSumMul(nPts, xpsi, &psi[j*nPts], dx) / (energies[j]-energies[i]);
	}

	// calculate density change
	double* drho = (double*) sq_malloc(sizeof(double)*nPts);
	std::complex<double>* ppsi_nc = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* ppsi_cc = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* temp = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	std::fill_n(drho, nPts, 0.0);
	for(int i = 0; i < nElec-1; i++){
		std::fill_n(ppsi_nc, nPts, 0.0);
		std::fill_n(ppsi_cc, nPts, 0.0);
		for(int j = i+1; j < nElec; j++){
			vtls::scaMulArray(nPts, mat[matIndex(i,j)]*weights[i], &psi[j*nPts], temp);
			vtls::addArrays(nPts, temp, ppsi_cc);

			vtls::scaMulArray(nPts, mat[matIndex(i,j)]*weights[j], &psi[j*nPts], temp);
			vtls::addArrays(nPts, temp, ppsi_nc);
		}
		for(int k = 0; k < nPts; k++)
			drho[k] += std::real( std::conj(psi[i*nPts+k])*ppsi_cc[k] - psi[i*nPts+k]*std::conj(ppsi_nc[k]) );
	}

	double xsum = 0.0, sum = 0.0;
	for(int i = minPos; i < maxPos; i++){
		xsum += i*drho[i];
		sum += drho[i];
	}

	std::cout << "Surface position found to be at index ~" << xsum/sum << " of " << nPts << std::endl;

	sq_free(mat);
	sq_free(xpsi);
	sq_free(drho);
	sq_free(ppsi_nc);
	sq_free(ppsi_cc);
	sq_free(temp);
	sq_free(energies);

	return (int)(xsum/sum);
}