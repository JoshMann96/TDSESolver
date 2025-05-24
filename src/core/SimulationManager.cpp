#include "SimulationManager.h"
#include "PhysCon.h"
#include "MathTools.h"

SimulationManager::SimulationManager(size_t nPts, double xMin, double dx, double dt, std::function<void(double)> callback, size_t numCallbackCalls)
	: dx(dx), nPts(nPts), dt(dt), progTracker(callback, numCallbackCalls), nElec(0)
{
	index = cyclic_int<size_t>(0, HISTORY_LENGTH);

	pot = new Potentials::PotentialManager(nPts);
	meas = new Measurers::MeasurementManager();
	psis = (std::complex<double>**) sq_malloc(sizeof(std::complex<double>*)*HISTORY_LENGTH);
	for(size_t i = 0; i < HISTORY_LENGTH; i++)
		psis[i] = nullptr;
	wavefunctionInitialized = false;

	vs = (double**) sq_malloc(sizeof(double*)*HISTORY_LENGTH);
	rhos = (double**) sq_malloc(sizeof(double*)*HISTORY_LENGTH);
	ts = (double*) sq_malloc(sizeof(double)*HISTORY_LENGTH);
	for (size_t i = 0; i < HISTORY_LENGTH; i++){
		vs[i] = (double*) sq_malloc(sizeof(double) * nPts);
		rhos[i] = (double*) sq_malloc(sizeof(double) * nPts);
	}
	std::fill_n(ts, HISTORY_LENGTH, 0.0);

	step = (size_t*) sq_malloc(sizeof(size_t) * HISTORY_LENGTH);
	std::fill_n(step, HISTORY_LENGTH, 0);

	scratch1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
	scratch2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
	spatialDamp = (double*) sq_malloc(sizeof(double) * nPts);
	std::fill_n(spatialDamp, nPts, 1.0);
	/*
	SimulationManager::maxT = maxT; SimulationManager::dt = dt; SimulationManager::dx = dx; SimulationManager::nPts = nPts;
	SimulationManager::mpiRoot = mpiRoot; SimulationManager::mpiUpdateTag = mpiUpdateTag; SimulationManager::mpiJob = mpiJob;
	*/

	x = (double*) sq_malloc(sizeof(double) * nPts);
	for(size_t i = 0; i < nPts; i++)
		x[i] = xMin + i*dx;
}

SimulationManager::~SimulationManager()
{
	for(size_t i = 0; i < HISTORY_LENGTH; i++){
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
	sq_free(step);

	sq_free(x);
	
	if(weights)
		sq_free(weights);

	delete meas;
	delete pot;	
}

void SimulationManager::addMeasurer(Measurers::Measurer* m) {
	meas->addMeasurer(m);
	if (m->needsDensity())
		calcDensity = true;
}

void SimulationManager::addPotential(Potentials::Potential* p) {
	pot->addPotential(p);
	if(p->getComplexity() == Potentials::PotentialComplexity::WAVEFUNCTION_DEPENDENT)
		calcDensity = true;
}

void SimulationManager::addSpatialDamp(const double* arr) {
	vtls::seqMulArrays(nPts, arr, spatialDamp);
}

void SimulationManager::calcEnergies(size_t curStep, double* energies) const {
	assert(wavefunctionInitialized);
	for(size_t i = 0; i < HISTORY_LENGTH; i++){
		if(curStep == step[i]){ //look for the present step's index
			double* rho = (double*) sq_malloc(sizeof(double)*nPts);
			for(size_t j = 0; j < nElec; j++){
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
		throw std::runtime_error("SimulationManager::calcWeights: Number of electrons is not finite! Failed to initialize.");

	if(weights)
		sq_free(weights); weights = nullptr;
	weights = (double*) sq_malloc(sizeof(double)*nElec);

	if (wght == nullptr){
		std::cout << "No weight function set! Using default of 1.0 for all states." << std::endl;
		for (size_t i = 0; i < nElec; i++)
			weights[i] = 1.0;
	}
	else{
		assert(wavefunctionInitialized);

		double* energies = (double*) sq_malloc(sizeof(double)*nElec);
		calcEnergies(step[index], energies);
		wght->calcWeights(nElec, energies, weights, normScheme);
		sq_free(energies);
	}
}

void SimulationManager::findEigenStates(double emin, double emax) {
	normScheme = WfcToRho::NormalizationScheme::NORMALIZED;
	
	pot->getVBare(0.0, vs[index]);

	std::complex<double>* states;

	kin->findEigenStates(vs[index], emin, emax, &states, &sq_malloc, &nElec);

	freePsis();
	psis[0] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

	vtls::copyArray(nPts * nElec, states, psis[0]);

	sq_free(states);

	for (size_t i = 0; i < nElec; i++)
		vtls::normalizeSqrNorm(nPts, &psis[0][i * nPts], dx);

	for (size_t i = 1; i < HISTORY_LENGTH; i++) {
		psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
		vtls::copyArray(nPts * nElec, psis[0], psis[i]);
	}

	wavefunctionInitialized = true;

	calcWeights();
	if(calcDensity){
		if(dens == nullptr)
			throw std::runtime_error("SimulationManager::findEigenStates: Density not set!");
		dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
	}
}

void SimulationManager::findInhomogeneousEigenStates(size_t nElec, const double* energies){
	KineticOperators::KineticOperator_FDM* kin_fdm = dynamic_cast<KineticOperators::KineticOperator_FDM*>(kin);
	if(kin_fdm == nullptr)
		throw std::runtime_error("SimulationManager::findInhomogeneousEigenStates: Kinetic operator is not a finite difference method!");
	
	this->nElec = nElec;

	freePsis();
	for(size_t i = 0; i < HISTORY_LENGTH; i++){
		psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
		std::fill_n(psis[i], nPts*nElec, 0.0);
		pot->getVBare(0.0, vs[i]);
	}

	kin_fdm->findInhomogeneousEigenStates(vs[index], energies, psis[index], nElec);
	for (size_t i = 1; i < HISTORY_LENGTH; i++) 
		vtls::copyArray(nPts * nElec, psis[index], psis[index + i]);

	wavefunctionInitialized = true;

	try{
	calcWeights();
	}
	catch(std::runtime_error& e){
		std::complex<double>* smallWfcs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * 1000 * nElec);
		double* realWfcs = (double*) sq_malloc(sizeof(double) * 1000 * nElec);
		for(int i = 0; i < nElec; i++)
			vtls::linearInterpolateNoEdge(nPts, &psis[index][i*nPts], 1000, smallWfcs + i*1000);
		vtls::copyArrayRe(1000*nElec, smallWfcs, realWfcs);
		plotting::GNUPlotter gp(1000, nElec, x, realWfcs);
		sq_free(smallWfcs);
		sq_free(realWfcs);
		throw e;
	}
	if(calcDensity){
		if(dens == nullptr)
			throw std::runtime_error("SimulationManager::findInhomogeneousEigenStates: Density not set!");
		dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
	}

	normScheme = WfcToRho::NormalizationScheme::UNNORMALIZED;
}

void SimulationManager::setPsi(std::complex<double>* npsi, WfcToRho::NormalizationScheme norm) {
	normScheme = norm;

	if (!nElec) {
		nElec = 1;
		freePsis();
		for (size_t i = 0; i < HISTORY_LENGTH; i++) {
			psis[i] = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
		}
	}
	vtls::copyArray(nPts, npsi, psis[index]);

	if(normScheme == WfcToRho::NormalizationScheme::NORMALIZED)
		vtls::normalizeSqrNorm(nPts, psis[index], dx);

	if(calcDensity){
		if(dens == nullptr)
			throw std::runtime_error("SimulationManager::setPsi: Density not set!");
		dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
	}

	wavefunctionInitialized = true;
}

size_t SimulationManager::calculatePotential(double* rho, const std::complex<double>* psi, double t, double* v){
	auto strt = std::chrono::high_resolution_clock::now();
	if(calcDensity){
		if(dens == nullptr)
			throw std::runtime_error("SimulationManager::calculatePotential: Density not set!");
		dens->calcRho(nPts, nElec, dx, weights, psi, rho);
	}
	pot->getV(rho, psi, t, v);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

size_t SimulationManager::calculatePotentialFromRawRho(double* rho, const std::complex<double>* psi, double t, double* v){
	auto strt = std::chrono::high_resolution_clock::now();
	dens->calcRho(nPts, nElec, dx, rho);
	pot->getV(rho, nullptr, t, v);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

size_t SimulationManager::updatePotential(int idx) {return calculatePotential(rhos[idx], psis[idx], ts[idx], vs[idx]);}

size_t SimulationManager::measure(int idx) {
	auto strt = std::chrono::high_resolution_clock::now();
	meas->measure(step[idx], psis[idx], vs[idx], ts[idx]);
	auto end = std::chrono::high_resolution_clock::now();
	auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - strt);
	return dur.count();
}

//Run simulation using operator splitting Fourier method (applies potential as linear)
void SimulationManager::runEPS_U2TU(size_t nSteps) {
	assert(wavefunctionInitialized);

	KineticOperators::KineticOperator_PSM* kin_psm = dynamic_cast<KineticOperators::KineticOperator_PSM*>(kin);
	if(kin_psm == nullptr)
		throw std::runtime_error("SimulationManager::runEPS_U2TU: Kinetic operator is not a pseudospectral method!");
	assert(nElec > 0);

	auto rMeasure = &SimulationManager::measure;
	auto rUpdatePotential = &SimulationManager::updatePotential;
	std::future<size_t> fM, fUP;

	// initialize progress tracker
	progTracker.reset(nSteps);

	bool asyncCalc = canAsyncCalcPot();
	if(!asyncCalc)
		std::cout << "Warning: Potential is not wavefunction independent! It is recommended to use runEPS_UW2TUW to more accurately account for the nonlinearity." << std::endl;

	for(size_t i = 0; i < nSteps; i++){
		// evaluate potential n
		if(asyncCalc){
			if (i == 0)
				updatePotential(index); // directly calculate
			else
				fUP.get(); // gather result
			fUP = std::async(rUpdatePotential, this, index + 1); // start n+1
		}
		else{
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

	// perform last measurement
	if (!asyncCalc)
		updatePotential(index);
	measure(index);

	progTracker.update(nSteps);
}

//Run simulation using operator splitting Fourier method (applies potential as nonlinear, second potential phase is recalculated after propagation phase)
void SimulationManager::runEPS_UW2TUW(size_t nSteps) {
	assert(wavefunctionInitialized);

	KineticOperators::KineticOperator_PSM* kin_psm = dynamic_cast<KineticOperators::KineticOperator_PSM*>(kin);
	if(kin_psm == nullptr)
		throw std::runtime_error("SimulationManager::runEPS_UW2TUW: Kinetic operator is not a pseudospectral method!");
	assert(nElec > 0);

	// variables for the midpoint of step
	std::complex<double>* tpsi = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
	double* trho = (double*) sq_malloc(sizeof(double) * nPts);
	double* tv = (double*) sq_malloc(sizeof(double) * nPts);

	auto rMeasure = &SimulationManager::measure;
	std::future<size_t> fM;

	// initialize progress tracker
	progTracker.reset(nSteps);

	for(size_t i = 0; i < nSteps; i++){
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
	
	// collect remaining future
	fM.get();

	// perform last measurement
	updatePotential(index);
	measure(index);

	progTracker.update(nSteps);

	sq_free(tpsi);
	sq_free(trho);
	sq_free(tv);
}

void SimulationManager::runCN_L(size_t nSteps){
	assert(wavefunctionInitialized);
	
	KineticOperators::CrankNicolson* kin_cn = dynamic_cast<KineticOperators::CrankNicolson*>(kin);
	if(kin_cn == nullptr)
		throw std::runtime_error("SimulationManager::runCN_L: Kinetic operator is not CrankNicolson!");
	assert(nElec > 0);

	double* tv = (double*) sq_malloc(sizeof(double) * nPts);
	auto rMeasure = &SimulationManager::measure;
	auto rUpdatePotential = &SimulationManager::updatePotential;
	std::future<size_t> fM, fUP;

	// initialize progress tracker
	progTracker.reset(nSteps);

	bool asyncCalc = canAsyncCalcPot();
	if(!asyncCalc)
		std::cout << "Warning: Potential is not wavefunction independent! It is recommended to use runCN_NL to more accurately account for the nonlinearity." << std::endl;

	for(size_t i = 0; i < nSteps; i++){
		// evaluate potential n+1 (n is already calculated)
		if(asyncCalc){
			if (i == 0){ // evaluate potential n, n+1 at beginning
				updatePotential(index);
				updatePotential(index + 1);
			}
			else
				fUP.get();
			fUP = std::async(rUpdatePotential, this, index + 2); // get n+2 going
		}
		else{
			if(i == 0)
				updatePotential(index);
			updatePotential(index + 1);
		}
		// potentials n, n+1 are now calculated
		// calculate averaged potential
		vtls::addArrays(nPts, vs[index], vs[index + 1], tv);
		vtls::scaMulArray(nPts, 0.5, tv);

		// evalute n->n+1
		kin_cn->step(psis[index], tv, spatialDamp, psis[index+1], nElec);
		
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

	// perform last measurement
	if (!asyncCalc)
		updatePotential(index);
	measure(index);

	progTracker.update(nSteps);

	sq_free(tv);
}

void SimulationManager::runCN_NL(size_t nSteps){
	assert(wavefunctionInitialized);
	
	KineticOperators::CrankNicolson* kin_cn = dynamic_cast<KineticOperators::CrankNicolson*>(kin);
	if(kin_cn == nullptr)
		throw std::runtime_error("SimulationManager::runCN_NL: Kinetic operator is not CrankNicolson!");
	assert(nElec > 0);

	double* tv = (double*) sq_malloc(sizeof(double) * nPts);
	auto rMeasure = &SimulationManager::measure;
	std::future<size_t> fM;

	// initialize progress tracker
	progTracker.reset(nSteps);

	for(size_t i = 0; i < nSteps; i++){
		// evaluate potential n
		updatePotential(index);

		// estimate the density at the next step using the present potential (virtual step)
		kin_cn->stepVirtual(psis[index], vs[index], spatialDamp, psis[index+1], nElec);

		// evaluate estimated potential n+1
		// 	try to calculate the raw density from the device then post-process, otherwise calculate rho normally
		if (kin_cn->calcRawRhoByDevice(weights, rhos[index + 1], true)) // if the device can calculate raw density
			calculatePotentialFromRawRho(rhos[index + 1], psis[index + 1], ts[index] + dt, vs[index+1]);
		else // otherwise calculate density on CPU
			updatePotential(index + 1);

		// averaged potential
		vtls::addArrays(nPts, vs[index], vs[index + 1], tv);
		vtls::scaMulArray(nPts, 0.5, tv);

		// true step n -> n+1
		kin_cn->step(psis[index], tv, spatialDamp, psis[index + 1], nElec);

		// measure step n while n+1->n+2 begins
		if(i != 0)
			fM.get();
		fM = std::async(rMeasure, this, index);
		
		progTracker.update(i);

		iterateIndex();
	}

	// collect remaining futures
	fM.get();

	// perform last measurement
	updatePotential(index);
	measure(index);

	progTracker.update(nSteps);

	sq_free(tv);
}

void SimulationManager::iterateIndex() {
	step[index+1] = step[index] + 1;
	ts[index+1] = step[index+1] * dt;

	index++;
}

size_t SimulationManager::findElectricalSurfaceCentroidRule(size_t minPos, size_t maxPos){
	/*
	* Calculate the electrical centroid of the electron density using first-order perturbation theory.
	*/
	assert(wavefunctionInitialized);

	std::complex<double>* mat = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nElec*(nElec-1)/2);
	std::complex<double>* xpsi = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	auto matIndex = [this](size_t i, size_t j){return i*nElec+j-((i+1)*(i+2))/2;}; //helper function for packing/unpacking matrix

	std::complex<double>* psi = psis[index];

	double* energies = (double*) sq_malloc(sizeof(double)*nElec);
	calcEnergies(step[index], energies);

	double* idxs = (double*) sq_malloc(sizeof(double)*nPts);

	// calculate matrix elements
	for(size_t i = 0; i < nElec-1; i++){
		for(size_t k = 0; k < nPts; k++)
			xpsi[k] = ((double)k) * psi[i*nPts+k];
		for(size_t j = i+1; j < nElec; j++)
			mat[matIndex(i,j)] = vtlsInt::rSumMul(nPts, xpsi, &psi[j*nPts], dx) / (energies[j]-energies[i]);
	}

	// calculate density change
	double* drho = (double*) sq_malloc(sizeof(double)*nPts);
	std::complex<double>* ppsi_nc = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* ppsi_cc = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* temp = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
	std::fill_n(drho, nPts, 0.0);
	for(size_t i = 0; i < nElec-1; i++){
		std::fill_n(ppsi_nc, nPts, 0.0);
		std::fill_n(ppsi_cc, nPts, 0.0);
		for(size_t j = i+1; j < nElec; j++){
			vtls::scaMulArray(nPts, mat[matIndex(i,j)]*weights[i], &psi[j*nPts], temp);
			vtls::addArrays(nPts, temp, ppsi_cc);

			vtls::scaMulArray(nPts, mat[matIndex(i,j)]*weights[j], &psi[j*nPts], temp);
			vtls::addArrays(nPts, temp, ppsi_nc);
		}
		for(size_t k = 0; k < nPts; k++)
			drho[k] += std::real( std::conj(psi[i*nPts+k])*ppsi_cc[k] - psi[i*nPts+k]*std::conj(ppsi_nc[k]) );
	}

	double xsum = 0.0, sum = 0.0;
	for(size_t i = minPos; i < maxPos; i++){
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

	return (size_t)(xsum/sum);
}