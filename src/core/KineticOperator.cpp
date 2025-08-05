#include "KineticOperator.h"
#include "PhysCon.h"
#include "fftw3.h"
#include <omp.h>
#include <fftw3.h>
#include "CORECommonHeader.h"
#include "blas.h"
#include "MathTools.h"

namespace KineticOperators {

	GenDisp_PSM::~GenDisp_PSM(){
		if (osKineticPhase)
			sq_free(osKineticPhase);
		if (osPotentialPhase)
			sq_free(osPotentialPhase);
		if (opMat)
			sq_free(opMat);
		if (osKineticEnergy)
			sq_free(osKineticEnergy);
		if (temp1)
			sq_free(temp1);
		if (temp2)
			sq_free(temp2);
		if (temp3)
			sq_free(temp3);
		if (groupVel)
			sq_free(groupVel);
		if (psik)
			sq_free(psik);
		

		mtx.lock();
		if(fftwOneForward)
			fftw_destroy_plan(fftwOneForward);
		if(fftwOneBackward)
			fftw_destroy_plan(fftwOneBackward);
		if(fftwAllForward)
			fftw_destroy_plan(fftwAllForward);
		if(fftwAllBackward)
			fftw_destroy_plan(fftwAllBackward);
		mtx.unlock();
	}

	void GenDisp_PSM::stepOS_U2TU(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (size_t i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//Perform FFT, apply full momentum-space phase contribution, invert FFT
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);

		//Apply half of phase contribution from potential, apply DFT normalization
#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osPotentialPhase[j];
			}
		}

	}

	void GenDisp_PSM::stepOS_UW2T(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (size_t i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//Perform FFT, apply full momentum-space phase contribution, invert FFT
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);
	}

	void GenDisp_PSM::stepOS_UW(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) {
		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (size_t i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}
	}

	void GenDisp_PSM::initializeAllFFT(size_t nElec) {
		if (firstStepAll || this->plan_nElec != nElec) {
			plan_nElec = nElec;
			/*DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nElec);
			DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
			DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
			DftiCommitDescriptor(dftiHandle);*/

			//initialize FFTW for performance, find best algo
			std::complex<double>* test = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

			mtx.lock();
			if(fftwAllForward)
				fftw_destroy_plan(fftwAllForward); fftwAllForward=NULL;
			if(fftwAllBackward)
				fftw_destroy_plan(fftwAllBackward);fftwAllBackward=NULL;
			
			fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			assert(nPts <= INT_MAX);
			int nPts_int = static_cast<int>(nPts);

			fftwAllForward = fftw_plan_many_dft(1, &nPts_int, nElec, reinterpret_cast<fftw_complex*>(test), &nPts_int, 1, nPts_int, reinterpret_cast<fftw_complex*>(test), &nPts_int, 1, nPts, FFTW_FORWARD, fftwPlanPolicy);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts_int, nElec, reinterpret_cast<fftw_complex*>(test), &nPts_int, 1, nPts_int, reinterpret_cast<fftw_complex*>(test), &nPts_int, 1, nPts, FFTW_BACKWARD, fftwPlanPolicy);

			mtx.unlock();

			if(!fftwAllForward || !fftwAllBackward)
				throw std::runtime_error("FFTW \"all\" plan creation failed");
			
			sq_free(test);

			if (firstStepAll) {
				if (osPotentialPhase)
					sq_free(osPotentialPhase); osPotentialPhase = nullptr;
				if (osKineticPhase)
					sq_free(osKineticPhase); osKineticPhase = nullptr;

				//initialize phase multipliers
				osPotentialPhase = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
				osKineticPhase = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
				for (size_t i = 0; i < nPts; i++)
					osKineticPhase[i] = std::exp(-PhysCon::im * dt / PhysCon::hbar * osKineticEnergy[i]);
			}

			firstStepAll = false;
		}
	}

	void GenDisp_PSM::initializeOneFFT() {
		if (firstStepOne) {
			fftw_complex *temp = (fftw_complex*) sq_malloc(sizeof(fftw_complex) * nPts);

			mtx.lock();

			if(fftwOneForward)
				fftw_destroy_plan(fftwOneForward); fftwOneForward=NULL;
			if(fftwOneBackward)
				fftw_destroy_plan(fftwOneBackward); fftwOneBackward=NULL;

			fftw_plan_with_nthreads(1);
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			assert(nPts <= INT_MAX);
			int nPts_int = static_cast<int>(nPts);

			fftwOneForward = fftw_plan_dft(1, &nPts_int, temp, temp, FFTW_FORWARD, FFTW_ESTIMATE);
			fftwOneBackward = fftw_plan_dft(1, &nPts_int, temp, temp, FFTW_BACKWARD, FFTW_ESTIMATE);

			mtx.unlock();

			if(!fftwOneForward || !fftwOneBackward)
				throw std::runtime_error("FFTW \"one\" plan creation failed");

			firstStepOne = false;

			sq_free(temp);
		}
	}

	void GenDisp_PSM::executeAllFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void GenDisp_PSM::executeAllFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
#pragma omp parallel for
		for(size_t i = 0; i < nElec; i++)
			vtls::scaMulArray(nPts, 1.0/nPts, &targ[i*nPts]);
	}

	void GenDisp_PSM::executeOneFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwOneForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void GenDisp_PSM::executeOneFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwOneBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
		vtls::scaMulArray(nPts, 1.0/nPts, targ);
	}

	void GenDisp_PSM::calcOpMat() {
		if(needMat){
			needMat = false;
			if (nPts > 46340) 
				throw std::runtime_error("Long datatype is required for grids of size nPts>46340. Rewrite this code (GenDisp_PSM::calcOpMat)");

			initializeOneFFT();

			if(opMat)
				sq_free(opMat); opMat = nullptr;
			opMat = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts*(nPts+1))/2);

			std::complex<double>* kinDiags = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

			vtls::copyArray(nPts, osKineticEnergy, kinDiags);

			executeOneFFTBackward(kinDiags);

			for (size_t d = 0; d < nPts; d++) {
				std::complex<double> cv = kinDiags[d];
				for (size_t i = 0; i < nPts - d; i++)
					opMat[(i * i + (2 * d + 3) * i + d * (d + 1)) / 2] = cv;
			}

			sq_free(kinDiags);

		}
	}

	void GenDisp_PSM::findEigenStates(const double* v, double emin, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs) {
		assert(nPts <= std::sqrt(LAPACK_INT_MAX));
			
		calcOpMat();
		for (size_t i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];


		*states = (std::complex<double>*)allocator(sizeof(std::complex<double>) * nPts * nPts);

		dcomplex * work = (dcomplex *)sq_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)sq_malloc(sizeof(double)*7*nPts);
		lapack_int * iwork3 = (lapack_int *)sq_malloc(sizeof(lapack_int)*5*nPts);
		double * eigs = (double *)sq_malloc(sizeof(double)*nPts);
		lapack_int * ifail = (lapack_int *)sq_malloc(sizeof(lapack_int)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		lapack_int info;
		lapack_int nPts_int = static_cast<lapack_int>(nPts);
		lapack_int nEigs_int;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts_int, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0, 0, &prec, &nEigs_int, eigs, reinterpret_cast<dcomplex *>(*states), &nPts_int, work, work2, iwork3, ifail, &info);

		if (info != 0) {
			std::cerr << "Error in LAPACK_zhpevx: " << info << std::endl;
			throw std::runtime_error("LAPACK_zhpevx failed");
		}

		freeOpMat();

		*nEigs = static_cast<size_t>(nEigs_int);
		nElec = *nEigs;

		sq_free(work);
		sq_free(work2);
		sq_free(iwork3);
		sq_free(eigs);
		sq_free(ifail);
	}

	void GenDisp_PSM::findGroundState(const double* v, size_t maxStates, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs) {
		initializeAllFFT(maxStates);

		double stop_thresh = 1e-8;
		double conv_thresh = 1e-10;
		size_t max_its = 1000;
		size_t its = 0;

		std::complex<double> *vecs = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts * maxStates);
		vtls::Orthonormalizer<lapack_complex_double> ortho(nPts, maxStates);

		double *energies = (double*)sq_malloc(sizeof(double) * maxStates);
		double *energysqr = (double*)sq_malloc(sizeof(double) * maxStates);
		double *vProp = (double*)sq_malloc(sizeof(double) * nPts);
		double *kProp = (double*)sq_malloc(sizeof(double) * nPts);
		double my_dt, maxEnergy;
		double relconv = 10.0*stop_thresh;

		// fill with noise in reciprocal space
		for (size_t i = 0; i < nPts * maxStates; i++)
			vecs[i] = (std::rand()/(double)RAND_MAX)-0.5;
		executeAllFFTBackward(vecs);
		ortho.orthonormalize(reinterpret_cast<lapack_complex_double*>(vecs));
		for (size_t i = 0; i < maxStates; i++)
			energies[i] = evaluateEnergy(&vecs[i*nPts], v);

		// get max energy, calculate dt
		maxEnergy = std::max(vtls::max(maxStates, energies), std::abs(vtls::min(maxStates, energies)));
		my_dt = dt;

		double* dens = (double*)sq_malloc(sizeof(double) * nPts);
		double* temp = (double*)sq_malloc(sizeof(double) * nPts);

		while(relconv > stop_thresh && its < max_its) {
			// calculate imaginary propagators
			for (size_t i = 0; i < nPts; i++) {
				vProp[i] = std::exp(- 0.5 * my_dt / PhysCon::hbar * v[i]);
				kProp[i] = std::exp(- my_dt / PhysCon::hbar * osKineticEnergy[i].real());
			}

			// propagate 10 times
			for (size_t j = 0; j < 10; j++) {
				// 1/2 potential
		#pragma omp parallel for
				for (size_t i = 0; i < maxStates; i++)
					vtls::seqMulArrays(nPts, vProp, &vecs[i * nPts]);
				
				// kinetic
				executeAllFFTForward(vecs);
		#pragma omp parallel for
				for (size_t i = 0; i < maxStates; i++)
					vtls::seqMulArrays(nPts, kProp, &vecs[i * nPts]);
				executeAllFFTBackward(vecs);

				// 1/2 potential
		#pragma omp parallel for
				for (size_t i = 0; i < maxStates; i++)
					vtls::seqMulArrays(nPts, vProp, &vecs[i * nPts]);
			}

			// orthonormalize
			ortho.orthonormalize(reinterpret_cast<lapack_complex_double*>(vecs));

			// calculate new energies
			for (size_t i = 0; i < maxStates; i++)
				energies[i] = evaluateEnergy(&vecs[i*nPts], v);
			for (size_t i = 0; i < maxStates; i++)
				energysqr[i] = evaluateEnergySquared(&vecs[i*nPts], v);
			// RMS
			double rms = 0.0;
			for (size_t i = 0; i < maxStates; i++)
				rms += std::abs(energysqr[i] - energies[i] * energies[i]);
			rms = std::sqrt(rms / maxStates);

			// get max energy, calculate dt
			maxEnergy = std::max(vtls::max(maxStates, energies), std::abs(vtls::min(maxStates, energies)));
			// calculate relative convergence
			relconv = rms / maxEnergy;
			// new time step (TODO: THIS NEEDS MORE WORK FOR BETTER/FASTER CONVERGENCE)
			if (its < max_its / 2)
				my_dt = std::min((100.0*std::log(relconv/conv_thresh) + 1.0) * dt, PhysCon::hbar/maxEnergy);
			else
				my_dt = dt;

			its++;

			if (its % 50 == 0)
				std::cout << "\tIteration " << its << "/" << max_its << "\n\t\tRMS <H> error: " << relconv << std::endl;
		}

		// order states by energy
		size_t *idxs = (size_t*)sq_malloc(sizeof(size_t) * maxStates);
		vtls::sort_idxs(maxStates, energies, idxs);
		// find num eigenstates according to emax
		*nEigs = maxStates;
		for (size_t i = 0; i < maxStates; i++){
			if (energies[i] > emax) {
				*nEigs = i;
				break;
			}
		}

		this->nElec = *nEigs;
		if (*nEigs == maxStates)
			std::cout << "Found " << *nEigs << " eigenstates with energy below " << emax << ", no states above." << std::endl;
		else
			std::cout << "Found " << *nEigs << " eigenstates with energy below " << emax << std::endl;

		// copy results
		*states = (std::complex<double>*)allocator(sizeof(std::complex<double>) * nPts * (*nEigs));

		for (size_t i = 0; i < *nEigs; i++)
			vtls::copyArray(nPts, &(vecs[idxs[i] * nPts]), &((*states)[i * nPts]));

		sq_free(dens);
		sq_free(temp);

		sq_free(idxs);
		sq_free(vecs);
		sq_free(energies);
		sq_free(energysqr);
		sq_free(vProp);
		sq_free(kProp);
	}

	double GenDisp_PSM::evaluateEnergy(const std::complex<double>* psi, const double* v) {
		initializeOneFFT();

		if(!temp1)
			temp1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);
		if(!temp2)
			temp2 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);

		// kinetic energy
		vtls::copyArray(nPts, psi, temp1);
		//DftiComputeForward(dftiHandleKin, temp1);
		executeOneFFTForward(temp1);
		vtls::seqMulArrays(nPts, osKineticEnergy, temp1, temp2);
		for (size_t i = 0; i < nPts; i++)
			temp1[i] = std::conj(temp1[i]);

		double res = std::real(vtlsInt::innerProduct(nPts, temp1, temp2, 1.0) / vtls::getNorm(nPts, temp1, 1.0));

		// potential energy
		vtls::seqMulArrays(nPts, v, psi, temp1);
		res += std::real(vtlsInt::conjugateInnerProduct(nPts, psi, temp1, 1.0) / vtls::getNorm(nPts, psi, 1.0));

		return res;
	}

	double GenDisp_PSM::evaluateEnergySquared(const std::complex<double>* psi, const double* v) {
		initializeOneFFT();

		if(!temp1)
			temp1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);
		if(!temp2)
			temp2 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);
		if(!temp3)
			temp3 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);

		// calculate \psi* (T+V)(T+V) \psi
		// T \psi
		vtls::copyArray(nPts, psi, temp1);
		executeOneFFTForward(temp1);
		vtls::seqMulArrays(nPts, osKineticEnergy, temp1);
		executeOneFFTBackward(temp1);
		// V \psi
		vtls::seqMulArrays(nPts, v, psi, temp2);
		vtls::addArrays(nPts, temp2, temp1);

		// temp1 = (T+V)\psi
		// T (T+V) \psi
		vtls::copyArray(nPts, temp1, temp2);
		executeOneFFTForward(temp2);
		vtls::seqMulArrays(nPts, osKineticEnergy, temp2);
		executeOneFFTBackward(temp2);
		// V (T+V) \psi
		vtls::seqMulArrays(nPts, v, temp1, temp3);
		vtls::addArrays(nPts, temp3, temp2);

		// temp2 = (T+V)(T+V) \psi

		double res = std::real(vtlsInt::conjugateInnerProduct(nPts, psi, temp2, 1.0) / vtls::getNorm(nPts, psi, 1.0));

		return res;
	}

	void GenDisp_PSM::calcRawCurrent(const std::complex<double>* psi, const double* weights, double* current, size_t nElec) {
		assert(groupVel != nullptr);
		initializeAllFFT(nElec);
		if (!psik)
			psik = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

		// Fourier transform
		vtls::copyArrayConj(nPts * nElec, -1.0, psi, psik);
		executeAllFFTForward(psik);
		// apply group velocity
		for(size_t i = 0; i < nElec; i++)
			vtls::seqMulArrays(nPts, groupVel, &psik[i * nPts]);
		// inverse Fourier transform
		executeAllFFTBackward(psik);
		// individual currents
		vtls::seqMulArrays(nPts * nElec, psi, psik);
		// sum over all states, apsply weight
		std::fill_n(current, nPts, 0.0);
		for(size_t i = 0; i < nElec; i++)
			vtls::scaMulAddArraysRe(nPts, weights[i], &psik[i * nPts], current);
	}


	GenDisp_PSM_FreeElec::GenDisp_PSM_FreeElec(size_t nPts, double dx, double dt, double m_eff, uint fftwPlanPolicy) : GenDisp_PSM(nPts, dx, dt, fftwPlanPolicy) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);//new std::complex<double>[nPts];

		double dphs = PhysCon::hbar*PhysCon::hbar / (2.0 * PhysCon::me*m_eff) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[0] = 0;
		for (size_t i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dphs * (double)(i * i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			sq_free(osKineticEnergy); osKineticEnergy = nullptr;
	}

	GenDisp_PSM_Series::GenDisp_PSM_Series(size_t nPts, double dx, double dt, size_t nPoly, const double* polyCoeffs, uint fftwPlanPolicy) : GenDisp_PSM(nPts, dx, dt, fftwPlanPolicy) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		osKineticEnergy[0] = 0;
		for (size_t i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dk * (double)(i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}
		vtls::polyEval(nPts, osKineticEnergy, nPoly, polyCoeffs, osKineticEnergy);

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			sq_free(osKineticEnergy); osKineticEnergy = nullptr;
	}

	GenDisp_PSM_MathExpr::GenDisp_PSM_MathExpr(size_t nPts, double dx, double dt, std::string expr, uint fftwPlanPolicy) : GenDisp_PSM(nPts, dx, dt, fftwPlanPolicy) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
		double* ks = (double*) sq_malloc(nPts * sizeof(double));

		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		ks[0] = 0.0;
		for (size_t i = 1; i < nPts / 2 + 1; i++) {
			ks[i] = dk * (double)(i);
			ks[nPts - i] = -dk * (double)(i);
		}
		vtls::evalMathExpr(nPts, "k", ks, expr, osKineticEnergy);

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			sq_free(osKineticEnergy); osKineticEnergy = nullptr;
		if (ks)
			sq_free(ks); ks = nullptr;
	}

	NonUnifGenDisp_PSM::~NonUnifGenDisp_PSM(){
		if (osKineticEnergy)
			sq_free(osKineticEnergy);
		if (osPotentialPhase)
			sq_free(osPotentialPhase);
		if (opMat)
			sq_free(opMat);
		if (tempPsi)
			sq_free(tempPsi);
		if (tempPsiCum)
			sq_free(tempPsiCum);
		if (osKineticMask)
			sq_free(osKineticMask);
		if (norms)
			sq_free(norms);
		if (groupVel)
			sq_free(groupVel);
		if (psik)
			sq_free(psik);
		if (tempCur)
			sq_free(tempCur);

		mtx.lock();
		if(fftwOneForward)
			fftw_destroy_plan(fftwOneForward);
		if(fftwOneBackward)
			fftw_destroy_plan(fftwOneBackward);
		if(fftwAllForward)
			fftw_destroy_plan(fftwAllForward);
		if(fftwAllBackward)
			fftw_destroy_plan(fftwAllBackward);
		mtx.unlock();
	}

	void NonUnifGenDisp_PSM::stepOS_U2TU(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (size_t i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//get original norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (size_t i = 0; i < nElec; i++)
				norms[i] = vtls::getNorm(nPts, &targ[i*nPts], dx);

		//Apply each order of kinetic exponential
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

		vtls::copyArray(nPts * nElec, targ, tempPsiCum);
		for (size_t o = 1; o <= expOrder; o++) {
			for (size_t d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
#pragma omp parallel for
				for (size_t i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nElec]);

				//back to real space
				//DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTBackward(&tempPsi[d * nPts * nElec]);
				//apply mask
#pragma omp parallel for
				for (size_t i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);

				//back to recip space
				//DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTForward(&tempPsi[d * nPts * nElec]);
				//apply rest of kinetic energy
#pragma omp parallel for
				for (size_t i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
			}
			//reset cumulative psi to first contribution
			vtls::copyArray(nPts * nElec, tempPsi, tempPsiCum);

			//combine all the other new psi components
			for (size_t d = 1; d < nDisp; d++)
				vtls::addArrays(nPts * nElec, &tempPsi[d * nPts * nElec], tempPsiCum);

			//apply factor (becomes factorial with multiple applications)
			vtls::scaMulArray(nPts * nElec, (-PhysCon::im * dt / PhysCon::hbar) / (double)o, tempPsiCum);

			//and add contribution to result
			vtls::addArrays(nPts * nElec, tempPsiCum, targ);
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);
		//restore norm if requested, else apply DFT normalization
		if (forceNorm)
#pragma omp parallel for
			for (size_t i = 0; i < nElec; i++)
				vtls::setNorm(nPts, &targ[i * nPts], dx, norms[i]);

		//Apply half of phase contribution from potential
#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osPotentialPhase[j];
			}
		}
	}

	void NonUnifGenDisp_PSM::stepOS_UW2T(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (size_t i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//get original norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (size_t i = 0; i < nElec; i++)
				norms[i] = vtls::getNorm(nPts, &targ[i * nPts], dx);

		//Apply each order of kinetic exponential
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);
		vtls::copyArray(nPts * nElec, targ, tempPsiCum);
		for (size_t o = 1; o <= expOrder; o++) {
			for (size_t d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
				for (size_t i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
				//back to real space
				//DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTBackward(&tempPsi[d * nPts * nElec]);
				//apply mask
				for (size_t i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
				//back to recip space
				//DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTForward(&tempPsi[d * nPts * nElec]);
				//apply rest of kinetic energy
				for (size_t i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
			}
			//reset cumulative psi to first contribution
			vtls::copyArray(nPts * nElec, tempPsi, tempPsiCum);
			//combine all the other new psi components
			for (size_t d = 1; d < nDisp; d++)
				vtls::addArrays(nPts * nElec, &tempPsi[d * nPts * nElec], tempPsiCum);
			//apply factor (becomes factorial with multiple applications)
			vtls::scaMulArray(nPts * nElec, (-PhysCon::im * dt / PhysCon::hbar) / (double)o, tempPsiCum);
			//and add contribution to result
			vtls::addArrays(nPts * nElec, tempPsiCum, targ);
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);

		//restore norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (size_t i = 0; i < nElec; i++)
				vtls::setNorm(nPts, &targ[i * nPts], dx, norms[i]);
	}

	void NonUnifGenDisp_PSM::stepOS_UW(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) {
		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (size_t i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for
		for (size_t i = 0; i < nElec; i++) {
			for (size_t j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}
	}

	void NonUnifGenDisp_PSM::initializeAllFFT(size_t nElec) {
		if (firstStepAll || NonUnifGenDisp_PSM::nElec != nElec) {
			NonUnifGenDisp_PSM::nElec = nElec;

			std::complex<double>* temp = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

			mtx.lock();

			if(fftwAllForward)
				fftw_destroy_plan(fftwAllForward); fftwAllForward=NULL;
			if(fftwAllBackward)
				fftw_destroy_plan(fftwAllBackward); fftwAllBackward=NULL;

			fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;
			
			assert(nPts <= INT_MAX);
			int nPts_int = static_cast<int>(nPts);

			fftwAllForward = fftw_plan_many_dft(1, &nPts_int, nElec, reinterpret_cast<fftw_complex*>(temp), &nPts_int, 1, nPts, reinterpret_cast<fftw_complex*>(temp), &nPts_int, 1, nPts, FFTW_FORWARD, fftwPlanPolicy);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts_int, nElec, reinterpret_cast<fftw_complex*>(temp), &nPts_int, 1, nPts, reinterpret_cast<fftw_complex*>(temp), &nPts_int, 1, nPts, FFTW_BACKWARD, fftwPlanPolicy);
			
			mtx.unlock();
			
			if(!fftwAllForward || !fftwAllBackward)
				throw std::runtime_error("FFTW \"all\" plan creation failed");

			sq_free(temp);

			if (osPotentialPhase)
				sq_free(osPotentialPhase); osPotentialPhase = nullptr;
			if (tempPsi)
				sq_free(tempPsi); tempPsi = nullptr;
			if (tempPsiCum)
				sq_free(tempPsiCum); tempPsiCum = nullptr;
			if (norms)
				sq_free(norms); norms = nullptr;

			osPotentialPhase = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
			tempPsi = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec * nDisp);
			tempPsiCum = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
			norms = (double*) sq_malloc(sizeof(double)*nElec);

			firstStepAll = false;
		}
	}

	void NonUnifGenDisp_PSM::initializeOneFFT() {
		if (firstStepOne) {
			std::complex<double>* temp = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);

			mtx.lock();

			//initialize FFTW for performance, find best algo
			if(fftwOneForward)
				fftw_destroy_plan(fftwOneForward); fftwOneForward=NULL;
			if(fftwOneBackward)
				fftw_destroy_plan(fftwOneBackward); fftwOneBackward=NULL;

			fftw_plan_with_nthreads(1);
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			assert(nPts <= INT_MAX);
			int nPts_int = static_cast<int>(nPts);

			fftwOneForward = fftw_plan_dft(1, &nPts_int, reinterpret_cast<fftw_complex*>(temp), reinterpret_cast<fftw_complex*>(temp), FFTW_FORWARD, FFTW_ESTIMATE);
			fftwOneBackward = fftw_plan_dft(1, &nPts_int, reinterpret_cast<fftw_complex*>(temp), reinterpret_cast<fftw_complex*>(temp), FFTW_BACKWARD, FFTW_ESTIMATE);

			mtx.unlock();

			if(!fftwOneForward || !fftwOneBackward)
				throw std::runtime_error("FFTW \"one\" plan creation failed");

			firstStepOne = false;

			sq_free(temp);
		}
	}

	void NonUnifGenDisp_PSM::executeAllFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void NonUnifGenDisp_PSM::executeAllFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
#pragma omp parallel for
		for(size_t i = 0; i < nElec; i++)
			vtls::scaMulArray(nPts, 1.0/nPts, &targ[i*nPts]);
	}

	void NonUnifGenDisp_PSM::executeOneFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwOneForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void NonUnifGenDisp_PSM::executeOneFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwOneBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
		vtls::scaMulArray(nPts, 1.0/nPts, targ);
	}

	void NonUnifGenDisp_PSM::calcOpMat() {
		if (needMat) {
			needMat = 0;
			if (nPts > 46340) {
				std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (NonUnifGenDisp_PSM::calcOpMat)" << std::endl;
				return;
			}

			initializeOneFFT();
			if(opMat)
				sq_free(opMat); opMat = nullptr;
			opMat = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1)) / 2);
			std::fill_n(opMat, (nPts * (nPts + 1)) / 2, 0.0);

			std::complex<double>* kinDiags = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
			std::complex<double>* kinMat = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1)) / 2);
			std::complex<double>* temp = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1)) / 2);

			for (size_t d = 0; d < nDisp; d++) {
				vtls::copyArray(nPts, &osKineticEnergy[d*nPts], kinDiags);
				//DftiComputeBackward(dftiHandleMat, kinDiags);
				executeOneFFTBackward(kinDiags);

				for (size_t dk = 0; dk < nPts; dk++) {
					std::complex<double> cv = kinDiags[dk];
					for (size_t i = 0; i < nPts - dk; i++)
						kinMat[(i * i + (2 * dk + 3) * i + dk * (dk + 1)) / 2] = cv;
				}

				vtls::mulHermitDiagHermit(nPts, kinMat, &osKineticMask[d * nPts], temp);
				vtls::addArrays((nPts * (nPts + 1)) / 2, temp, opMat);
			}

			if (kinDiags)
				sq_free(kinDiags); kinDiags = nullptr;
			if (kinMat)
				sq_free(kinMat); kinMat = nullptr;
			if (temp)
				sq_free(temp); temp = nullptr;

		}
	}

	void NonUnifGenDisp_PSM::findEigenStates(const double* v, double emin, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs) {
		calcOpMat();
		for (size_t i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];

		*states = (std::complex<double>*)allocator(sizeof(std::complex<double>) * nPts * nPts);

		dcomplex * work = (dcomplex *)sq_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)sq_malloc(sizeof(double)*7*nPts);
		lapack_int * iwork3 = (lapack_int *)sq_malloc(sizeof(lapack_int)*5*nPts);
		double * eigs = (double *)sq_malloc(sizeof(double)*nPts);
		lapack_int * ifail = (lapack_int *)sq_malloc(sizeof(lapack_int)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		lapack_int info;
		lapack_int nPts_int = static_cast<lapack_int>(nPts);
		lapack_int nEigs_int;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts_int, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0, 0, &prec, &nEigs_int, eigs, reinterpret_cast<dcomplex *>(*states), &nPts_int, work, work2, iwork3, ifail, &info);

		if (info != 0) {
			std::cerr << "Error in LAPACK_zhpevx: " << info << std::endl;
			throw std::runtime_error("LAPACK_zhpevx failed");
		}

		freeOpMat();

		nElec = nEigs[0];

		if (work)
			sq_free(work); work = nullptr;
		if (work2)
			sq_free(work2); work2 = nullptr;
		if (iwork3)
			sq_free(iwork3); iwork3 = nullptr;
		if (eigs)
			sq_free(eigs); eigs = nullptr;
		if (ifail)
			sq_free(ifail); ifail = nullptr;
	}

	double NonUnifGenDisp_PSM::evaluateEnergy(const std::complex<double>* psi, const double* v) {
		initializeOneFFT();

		std::complex<double>* temp1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);
		std::complex<double>* temp2 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);
		std::complex<double>* temp3 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);

		vtls::copyArray(nPts, psi, temp1);
		//DftiComputeForward(dftiHandleKin, temp1);
		executeOneFFTForward(temp1);

		std::fill_n(temp2, nPts, 0.0);
		for (size_t d = 0; d < nDisp; d++) {
			vtls::seqMulArrays(nPts, &osKineticEnergy[d*nPts], temp1, temp3);
			//DftiComputeBackward(dftiHandleKin, temp3);
			executeOneFFTBackward(temp3);
			vtls::seqMulArrays(nPts, &osKineticMask[d*nPts], temp3);
			//DftiComputeForward(dftiHandleKin, temp3);
			executeOneFFTForward(temp3);
			vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], temp3);
			vtls::addArrays(nPts, temp3, temp2);
		}
		for (size_t i = 0; i < nPts; i++)
			temp1[i] = std::conj(temp1[i]);

		double res = std::real(vtlsInt::innerProduct(nPts, temp1, temp2, 1.0) / vtls::getNorm(nPts, temp1, 1.0));

		// potential energy
		vtls::copyArray(nPts, psi, temp1);
		vtls::seqMulArrays(nPts, v, temp1);
		res += std::real(vtlsInt::conjugateInnerProduct(nPts, psi, temp1, 1.0) / vtls::getNorm(nPts, temp1, 1.0));

		sq_free(temp1);
		sq_free(temp2);
		sq_free(temp3);

		return res;
	}

	void NonUnifGenDisp_PSM::calcRawCurrent(const std::complex<double>* psi, const double* weights, double* current, size_t nElec) {
		assert(groupVel != nullptr);
		initializeAllFFT(nElec);
		if (!psik)
			psik = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts * nElec);
		if (!tempCur)
			tempCur = (double*)sq_malloc(sizeof(double) * nPts);

		std::fill_n(current, nPts, 0.0);
		for(size_t d = 0; d < nDisp; d++){
			// Fourier transform
			vtls::copyArrayConj(nPts * nElec, psi, psik);
			executeAllFFTForward(psik);
			// apply group velocity
			for (size_t i = 0; i < nElec; i++)
				vtls::seqMulArrays(nPts, &groupVel[d*nPts], &psik[i * nPts]);
			// inverse Fourier transform
			executeAllFFTBackward(psik);
			// individual currents
			vtls::seqMulArrays(nPts * nElec, psi, psik);
			// sum over all states, apply weight
			std::fill_n(tempCur, nPts, 0.0);
			for (size_t i = 0; i < nElec; i++)
				vtls::scaMulAddArraysRe(nPts, weights[i], &psik[i * nPts], tempCur);
			vtls::seqMulAddArrays(nPts, &osKineticMask[d * nPts], tempCur, current);
		}
	}


	NonUnifGenDisp_PSM_EffMassBoundary::NonUnifGenDisp_PSM_EffMassBoundary(size_t nPts, double dx, double dt, size_t expOrder, bool forceNormalization, double meff_l, double meff_r, double transRate, size_t transPos, double edgeRate, uint fftwPlanPolicy) : NonUnifGenDisp_PSM(nPts, dx, dt, 2, expOrder, forceNormalization, fftwPlanPolicy) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*2);
		double* mask = (double*)sq_malloc(sizeof(double)*nPts*2);

		double dphs = PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me * meff_r) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[0] = 0;
		for (size_t i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dphs * (double)(i * i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}

		dphs = PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me * meff_l) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[nPts] = 0;
		for (size_t i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[nPts + i] = dphs * (double)(i * i);
			osKineticEnergy[nPts + nPts - i] = osKineticEnergy[i];
		}

		if (edgeRate != 0.0) {
			for (size_t i = 0; i < nPts; i++) {
				mask[i] = (1.0 / (1.0 + std::exp(-dx * transRate * (i - transPos))) + 1.0 / (1.0 + std::exp(dx * edgeRate * i))) / (1.0 + std::exp(dx * edgeRate * (i - nPts)));
				mask[i + nPts] = 1.0 - mask[i];
			}
		}
		else {
			for (size_t i = 0; i < nPts; i++) {
				mask[i] = 1.0 / (1.0 + std::exp(-dx * transRate * (i - transPos)));
				mask[i + nPts] = 1.0 - mask[i];
			}
		}

		NonUnifGenDisp_PSM::set_osKineticEnergy(osKineticEnergy, mask);
		if (osKineticEnergy)
			sq_free(osKineticEnergy); osKineticEnergy = nullptr;
		if (mask)
			sq_free(mask); mask = nullptr;
	}

	NonUnifGenDisp_PSM_MathExprBoundary::NonUnifGenDisp_PSM_MathExprBoundary(size_t nPts, double dx, double dt, size_t expOrder, bool forceNormalization, size_t nDisp, std::vector<std::string> exprs, double* transRates, size_t* transPoss, uint fftwPlanPolicy) : NonUnifGenDisp_PSM(nPts, dx, dt, nDisp, expOrder, forceNormalization, fftwPlanPolicy) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nDisp);
		double* mask = (double*) sq_malloc(nPts * nDisp * sizeof(double));
		std::fill_n(mask, nPts * nDisp, 1.0);
		double* ks = (double*) sq_malloc(nPts * sizeof(double));

		// generate ks
		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		ks[0] = 0.0;
		for (size_t i = 1; i < nPts / 2 + 1; i++) {
			ks[i] = dk * (double)(i);
			ks[nPts - i] = -dk * (double)(i);
		}

		// generate osKineticEnergy
		for(size_t i = 0; i < nDisp; i++)
			vtls::evalMathExpr(nPts, "k", ks, exprs[i], &osKineticEnergy[i*nPts]);

		// generate sigmoid masks
		for (size_t i = 0; i < nDisp; i++) {
			for (size_t j = 0; j < nPts; j++) {
				for(size_t k = 0; k < i; k++)
					mask[i * nPts + j] *= 1.0 / (1.0 + std::exp(-dx * transRates[k] * (j - transPoss[k])));
				for (size_t k = i; k < nDisp - 1; k++) 
					mask[i * nPts + j] *= 1.0 / (1.0 + std::exp(dx * transRates[k] * (j - transPoss[k])));
			}
		}

		//set for base class
		NonUnifGenDisp_PSM::set_osKineticEnergy(osKineticEnergy, mask);
		if (osKineticEnergy)
			sq_free(osKineticEnergy); osKineticEnergy = nullptr;
		if (mask)
			sq_free(mask); mask = nullptr;
	}


	void KineticOperator_FDM::projectHistory(const std::complex<double>* psi, const std::complex<double>* phsL, const std::complex<double>* phsR, const double* v, size_t nElec) {
		std::complex<double>* bcwfs = ( std::complex<double>* )sq_malloc(sizeof(std::complex<double>) * nElec);

		cblas_zcopy(nElec, &psi[0], nPts, bcwfs, 1);
		lbc->fillHistory(bcwfs, phsL, v[0]);

		cblas_zcopy(nElec, &psi[nPts-1], nPts, bcwfs, 1);
		rbc->fillHistory(bcwfs, phsR, v[nPts-1]);

		sq_free(bcwfs);
	}

	CrankNicolson::CrankNicolson(size_t nPts, double dx, double dt, double m_eff, FDBCs::BoundaryCondition* leftBC, FDBCs::BoundaryCondition* rightBC, bool useCuda) :
		KineticOperator_FDM(nPts, leftBC, rightBC), dx(dx), dt(dt), m_eff(m_eff), useCuda(useCuda) {
		d = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);
		ud= (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * (nPts-1));
		ld= (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * (nPts-1));


		lhsDiag0 = 1.0 + 0.5*PhysCon::im*PhysCon::hbar/PhysCon::me/m_eff*dt/(dx*dx);
		lhsOffDiag0 = -0.25*PhysCon::im*PhysCon::hbar/PhysCon::me/m_eff*dt/(dx*dx);

		rhsDiag0 = 1.0 - 0.5*PhysCon::im*PhysCon::hbar/PhysCon::me/m_eff*dt/(dx*dx);
		rhsOffDiag = 0.25*PhysCon::im*PhysCon::hbar/PhysCon::me/m_eff*dt/(dx*dx);

		std::fill_n(d, nPts, lhsDiag0);
		std::fill_n(ud, nPts-1, lhsOffDiag0);
		std::fill_n(ld, nPts-1, lhsOffDiag0);

		potCoef = 0.5*PhysCon::im*dt/PhysCon::hbar;
	}

	void CrankNicolson::_step(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec, bool isVirtual) {
		if(bct1 == nullptr)
			bct1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
		if(bct2 == nullptr)
			bct2 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
		if(lbct == nullptr)
			lbct = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
		if(rbct == nullptr)
			rbct = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
		#ifndef USE_CUDA
		if(useCuda)
			throw std::runtime_error("CUDA support not compiled in this build");
		#else // USE_CUDA
		if(useCuda && !cuSolver){
			cuSolver = new cudaTridiagonalSolverSystem(nPts, nElec);
			if(!r_d)
				r_d = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);

			std::fill_n(ud, nPts-1, lhsOffDiag0);
			std::fill_n(ld, nPts-1, lhsOffDiag0);
			ud[0] = lbc->getLHSAdjEle();
			ld[nPts-2] = rbc->getLHSAdjEle();
			cuSolver->setOffDiag(ld, ud, cudaTridiagonalSolverSystem::LHS);

			std::fill_n(ud, nPts-1, rhsOffDiag);
			std::fill_n(ld, nPts-1, rhsOffDiag);
			cuSolver->setOffDiag(ld, ud, cudaTridiagonalSolverSystem::RHS);
			cuSolver->setX(psi0);
		}
		#endif // USE_CUDA
			
		//prepare left BC
		cblas_zcopy(nElec, psi0, nPts, bct1, 1); // map first element of all wavefunctions to bct1
		cblas_zcopy(nElec, &psi0[1], nPts, bct2, 1); // map second element of all wavefunctions to bct2
	
		lbc->prepareStep(bct1, bct2, v[0]);
		lbc->getRHS(bct1, bct2, v[0], lbct, nElec);
		if(!isVirtual)
			lbc->finishStep(bct1, bct2, std::real(v[0]));

		//prepare right BC
		cblas_zcopy(nElec, &psi0[nPts-1], nPts, bct1, 1); // map last element of all wavefunctions to bct1
		cblas_zcopy(nElec, &psi0[nPts-2], nPts, bct2, 1); // map second to last element of all wavefunctions to bct2
		
		rbc->prepareStep(bct1, bct2, v[nPts-1]);
		rbc->getRHS(bct1, bct2, v[nPts-1], rbct, nElec);
		if(!isVirtual)
			rbc->finishStep(bct1, bct2, std::real(v[nPts-1]));

		//prepare LHS matrix
		std::fill_n(d, nPts, lhsDiag0);
		vtls::scaMulAddArrays(nPts-2, potCoef, &v[1], &d[1]); // d += potmul*v, leave BCs alone

		d[0] = lbc->getLHSEle();
		d[nPts-1] = rbc->getLHSEle();
		ud[0] = lbc->getLHSAdjEle();
		ld[nPts-2] = rbc->getLHSAdjEle();
		if(useCuda){
			#ifdef USE_CUDA
			cuSolver->setBdyCond(ud[0], lbct, cudaTridiagonalSolverSystem::LHS);
			cuSolver->setBdyCond(ld[nPts-2], rbct, cudaTridiagonalSolverSystem::RHS);

			// evaluate RHS
			std::fill_n(r_d, nPts, rhsDiag0);
			vtls::scaMulAddArrays(nPts, -potCoef, v, r_d); // r_d += potmul*v
			cuSolver->rhsProduct(r_d, isVirtual, false);
			#endif // USE_CUDA
		}
		else{
			std::fill_n(ud, nPts-1, lhsOffDiag0);
			std::fill_n(ld, nPts-1, lhsOffDiag0);
			ud[0] = lbc->getLHSAdjEle();
			ld[nPts-2] = rbc->getLHSAdjEle();

			// evaluate RHS
			#pragma omp parallel for collapse(2)
			for(size_t j = 0; j < nElec; j++)
				for(size_t k = 1; k < nPts-1; k++)
					targ[j*nPts+k] = (-potCoef*v[k]+rhsDiag0)*psi0[j*nPts+k] +
						(rhsOffDiag*psi0[j*nPts+k-1] + rhsOffDiag*psi0[j*nPts+k+1]);

			// apply RHS BC
			cblas_zcopy(nElec, lbct, 1, targ, nPts);
			cblas_zcopy(nElec, rbct, 1, &targ[nPts-1], nPts);
		}

		//SOLVE
		if(useCuda){
			#ifdef USE_CUDA
			cuSolver->solve(d, isVirtual, isVirtual);
			cuSolver->vectorHadamardProduct(spatialDamp, isVirtual); // apply spatial damping
			if(!isVirtual)
				cuSolver->gatherX(targ, false);
			#endif // USE_CUDA
		}
		else{
			lapack_int info;
			assert(nElec <= LAPACK_INT_MAX);
			lapack_int nElec_int = static_cast<lapack_int>(nElec);
			assert(nPts <= LAPACK_INT_MAX);
			lapack_int nPts_int = static_cast<lapack_int>(nPts);
			LAPACK_zgtsv(&nPts_int, &nElec_int, reinterpret_cast<dcomplex*>(ld), reinterpret_cast<dcomplex*>(d), reinterpret_cast<dcomplex*>(ud), reinterpret_cast<dcomplex*>(targ), &nPts_int, &info);
			
			if(info != 0) {
				std::cerr << "Error in LAPACK_zgtsv: " << info << std::endl;
				throw std::runtime_error("LAPACK_zgtsv failed");
			}

			for(size_t i = 0; i < nElec; i++) // apply spatial damping
				vtls::seqMulArrays(nPts, spatialDamp, &targ[i*nPts]);
		}

		// Using "expert" LAPACK driver (WORKSPACE MUST BE ALLOCATED, NOT A SIMPLE UNCOMMENT)
		// much slower
		//auto t1 = std::chrono::high_resolution_clock::now();
		//LAPACK_zgtsvx("N", "N", &nPts, &nElec,  reinterpret_cast<dcomplex*>(ld), reinterpret_cast<dcomplex*>(d), reinterpret_cast<dcomplex*>(ud), templ, tempd, tempu, tempu2, ipiv, reinterpret_cast<dcomplex*>(rhs), &nPts, reinterpret_cast<dcomplex*>(targ), &nPts, &rcond, ferr, berr, work, rwork, &info);
		//auto t2 = std::chrono::high_resolution_clock::now();
		//std::cout << "Time taken for LAPACK_zgtsvx: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;
	}

	void CrankNicolson::findEigenStates(const double* v, double emin, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs){
		double* hd = (double*)sq_malloc(sizeof(double)*nPts);
		double* hod= (double*)sq_malloc(sizeof(double)*(nPts-1));
		lapack_int  nSplit;
		lapack_int* iblock = (lapack_int*)sq_malloc(sizeof(lapack_int)*nPts);
		lapack_int* isplit = (lapack_int*)sq_malloc(sizeof(lapack_int)*nPts);
		lapack_int* iwork = (lapack_int*)sq_malloc(sizeof(lapack_int)*3*nPts);
		double* work = (double*)sq_malloc(sizeof(double)*5*nPts);
		double* eigs = (double*)sq_malloc(sizeof(double)*nPts);
		
		double energyScaler = PhysCon::me*m_eff*dx*dx/(PhysCon::hbar*PhysCon::hbar); // to make matrix nicely scaled
		emin = emin*energyScaler;
		emax = emax*energyScaler;

		std::fill_n(hd, nPts, 1.0);
		vtls::scaMulAddArrays(nPts, energyScaler, v, hd);
		std::fill_n(hod, nPts-1, -0.5);

		// get eigenvalues
		const char* cS = "S";
		double prec = 2.0*LAPACK_dlamch(cS);
		lapack_int info;
		assert(nPts <= LAPACK_INT_MAX);
		lapack_int nPts_int = static_cast<lapack_int>(nPts);
		lapack_int nEigs_int;
		LAPACK_dstebz("V", "B", &nPts_int, &emin, &emax, 0, 0, &prec, hd, hod, &nEigs_int, &nSplit, eigs, iblock, isplit, work, iwork, &info);

		if (info != 0) {
			std::cerr << "Error in LAPACK_dstebz: " << info << std::endl;
			throw std::runtime_error("LAPACK_dstebz failed");
		}

		*nEigs = static_cast<size_t>(nEigs_int);
		double* statesTemp = (double*)sq_malloc(sizeof(double)*nPts*(*nEigs));
		lapack_int* ifail = (lapack_int*)sq_malloc(sizeof(lapack_int)*(*nEigs));

		// get eigenvectors
		LAPACK_dstein(&nPts_int, hd, hod, &nEigs_int, eigs, iblock, isplit, statesTemp, &nPts_int, work, iwork, ifail, &info);

		if (info != 0) {
			std::cerr << "Error in LAPACK_dstein: " << info << std::endl;
			throw std::runtime_error("LAPACK_dstein failed");
		}

		// copy eigenvectors to states
		*states = (std::complex<double>*)allocator(sizeof(std::complex<double>)*nPts*(*nEigs));
		vtls::copyArray(nPts*(*nEigs), statesTemp, *states);
		
		sq_free(hd);
		sq_free(hod);
		sq_free(iblock);
		sq_free(isplit);
		sq_free(iwork);
		sq_free(work);
		sq_free(eigs);
		sq_free(statesTemp);
		sq_free(ifail);
	}

	void CrankNicolson::findInhomogeneousEigenStates(const double* v, const double* es, std::complex<double>* states, size_t nElec){
		std::complex<double>* lhs_d = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
		std::complex<double>* lhs_ld= (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
		std::complex<double>* lhs_ud= (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
		std::complex<double>* rhs   = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

		for(size_t i = 0; i < nElec; i++){
            fillInhomEigenMatrix(lhs_ld, lhs_ud, rhs, lhs_d, es[i], v);

            // check if the sytem is inhomogeneous
			if(std::abs(rhs[0]) < 1e-10 && std::abs(rhs[nPts-1]) < 1e-10)
				throw std::runtime_error("System must be inhomogeneous to use findInhomogeneousEigenStates");

			//SOLVE
			lapack_int info, one=1;
			assert(nPts <= LAPACK_INT_MAX);
			lapack_int nPts_int = static_cast<lapack_int>(nPts);

			LAPACK_zgtsv(&nPts_int, &one, reinterpret_cast<dcomplex*>(lhs_ld), reinterpret_cast<dcomplex*>(lhs_d), reinterpret_cast<dcomplex*>(lhs_ud), reinterpret_cast<dcomplex*>(rhs), &nPts_int, &info);
		
			if(info != 0) {
				std::cerr << "Error in LAPACK_zgtsv: " << info << std::endl;
				throw std::runtime_error("LAPACK_zgtsv failed");
			}

			vtls::copyArray(nPts, rhs, &states[i*nPts]);
		}

		// project history onto BCs
		std::complex<double>* phaseAdvancement = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nElec);
		for(size_t i = 0; i < nElec; i++)
			phaseAdvancement[i] = phaseAdvanceFromEnergy(es[i], dt);
		projectHistory(states, phaseAdvancement, phaseAdvancement, v, nElec);

		sq_free(lhs_d);
		sq_free(lhs_ld);
		sq_free(lhs_ud);
		sq_free(rhs);
		sq_free(phaseAdvancement);
    }

    void CrankNicolson::fillInhomEigenMatrix(std::complex<double> *lhs_ld, std::complex<double> *lhs_ud, std::complex<double> *rhs, std::complex<double> *lhs_d, double e, const double *v)
    {
        // (re)fill static matrix elements
        std::fill_n(lhs_ld, nPts - 1, -0.5);
        std::fill_n(lhs_ud, nPts - 1, -0.5);
        std::fill_n(rhs, nPts, 0.0);
        // fill main diagonal
        for (size_t j = 1; j < nPts - 1; j++)
            lhs_d[j] = 1.0 - 0.5 * (e - v[j]) / (PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me * m_eff * (dx * dx)));

        // apply BCs
        lhs_d[0] = lbc->getSteadyLHSEle(e - v[0]);
        lhs_d[nPts - 1] = rbc->getSteadyLHSEle(e - v[nPts - 1]);

        lhs_ud[0] = lbc->getSteadyLHSAdjEle(e - v[0]);
        lhs_ld[nPts - 2] = rbc->getSteadyLHSAdjEle(e - v[nPts - 1]);

        rhs[0] = lbc->getSteadyRHS(e - v[0]);
        rhs[nPts - 1] = rbc->getSteadyRHS(e - v[nPts - 1]);
    }

	double CrankNicolson::evaluateEnergy(const std::complex<double>* psi, const double* v){
		if(!tempPsi1)
			tempPsi1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

		double kineticCoef = PhysCon::hbar*PhysCon::hbar/(2.0*PhysCon::me*m_eff*dx*dx);

		// ignore left and right bdys, only use bdy values for derivatives
		vtls::seqMulArrays	 (nPts-2, v+1, psi+1, tempPsi1+1);
		vtls::scaMulAddArrays(nPts-2, 2.0*kineticCoef, psi+1, tempPsi1+1);
		vtls::scaMulAddArrays(nPts-2, -kineticCoef, psi+2, tempPsi1+1);
		vtls::scaMulAddArrays(nPts-2, -kineticCoef, psi, tempPsi1+1);

		return std::real(vtlsInt::conjugateInnerProduct(nPts-2, psi+1, tempPsi1+1, 1.0) / vtls::getNorm(nPts-2, psi+1, 1.0));
	}

	bool CrankNicolson::calcRawRhoByDevice(const double* weights, double* rho, bool virt){
		if(!useCuda)
			return false;
		#ifdef USE_CUDA
		cuSolver->calcRawRho(weights, rho, virt);
		#else // USE_CUDA
		throw std::runtime_error("CUDA support not compiled in this build");
		#endif // USE_CUDA
		
		return true;
	}

	bool CrankNicolson::calcRawCurByDevice(const double* weights, double* cur, bool virt){
		if(!useCuda)
			return false;
		#ifdef USE_CUDA
		cuSolver->calcRawCur(weights, cur, PhysCon::hbar/(PhysCon::me*m_eff*dx), virt);
		#else // USE_CUDA
		throw std::runtime_error("CUDA support not compiled in this build");
		#endif // USE_CUDA

		return true;
	}

	void CrankNicolson::calcRawCurrent(const std::complex<double>* psi, const double* weights, double* current, size_t nElec) {
		if(!tempPsi1)
			tempPsi1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);
		if(!tempPsi2)
			tempPsi2 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nPts);

		std::complex<double> dvdk = PhysCon::hbar/PhysCon::me/m_eff*PhysCon::im;

		std::fill_n(current, nPts, 0.0);
		for(size_t n = 0; n < nElec; n++) {
			vtls::copyArrayConj(nPts, &psi[n*nPts], tempPsi1);
			vtls::firstDerivative(nPts, tempPsi1, tempPsi2, dx);
			vtls::seqMulArrays(nPts, psi, tempPsi2);
			vtls::scaMulAddArraysRe(nPts, weights[n]*dvdk, tempPsi2, current);
		}
	}
}