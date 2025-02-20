#include "KineticOperator.h"
#include "PhysCon.h"
#include "fftw3.h"
#include <omp.h>
#include <fftw3.h>
#include "CORECommonHeader.h"
#include "blas.h"
#include "MathTools.h"

#define MULTIELEC_FFTW_POLICY FFTW_PATIENT

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

	void GenDisp_PSM::stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//Perform FFT, apply full momentum-space phase contribution, invert FFT
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);

		//Apply half of phase contribution from potential, apply DFT normalization
#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osPotentialPhase[j];
			}
		}

	}

	void GenDisp_PSM::stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//Perform FFT, apply full momentum-space phase contribution, invert FFT
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);
	}

	void GenDisp_PSM::stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) {
		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}
	}

	void GenDisp_PSM::initializeAllFFT(int nElec) {
		if (firstStepAll || GenDisp_PSM::nElec != nElec) {
			GenDisp_PSM::nElec = nElec;
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

			fftwAllForward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_FORWARD, MULTIELEC_FFTW_POLICY);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_BACKWARD, MULTIELEC_FFTW_POLICY);
			
			mtx.unlock();
			
			sq_free(test);

			if (firstStepAll) {
				if (osPotentialPhase)
					sq_free(osPotentialPhase); osPotentialPhase = nullptr;
				if (osKineticPhase)
					sq_free(osKineticPhase); osKineticPhase = nullptr;

				//initialize phase multipliers
				osPotentialPhase = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
				osKineticPhase = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
				for (int i = 0; i < nPts; i++)
					osKineticPhase[i] = std::exp(-PhysCon::im * dt / PhysCon::hbar * osKineticEnergy[i]);
			}

			firstStepAll = 0;
		}
	}

	void GenDisp_PSM::initializeOneFFT() {
		if (firstStepOne) {
			/*DftiCreateDescriptor(&dftiHandleKin, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleKin, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleKin);*/

			if (temp1)
				sq_free(temp1); temp1 = nullptr;
			if (temp2)
				sq_free(temp2); temp2 = nullptr;

			temp1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
			temp2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);

			mtx.lock();

			if(fftwOneForward)
				fftw_destroy_plan(fftwOneForward); fftwOneForward=NULL;
			if(fftwOneBackward)
				fftw_destroy_plan(fftwOneBackward); fftwOneBackward=NULL;

			fftw_plan_with_nthreads(1);
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			fftwOneForward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1), FFTW_FORWARD, FFTW_ESTIMATE);
			fftwOneBackward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2), FFTW_BACKWARD, FFTW_ESTIMATE);

			mtx.unlock();

			firstStepOne = 0;
		}
	}

	void GenDisp_PSM::executeAllFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void GenDisp_PSM::executeAllFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
#pragma omp parallel for
		for(int i = 0; i < nElec; i++)
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
			needMat = 0;
			if (nPts > 46340) 
				throw std::runtime_error("Long datatype is required for grids of size nPts>46340. Rewrite this code (GenDisp_PSM::calcOpMat)");

			initializeOneFFT();

			if(opMat)
				sq_free(opMat); opMat = nullptr;
			opMat = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts*(nPts+1))/2);

			std::complex<double>* kinDiags = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

			vtls::copyArray(nPts, osKineticEnergy, kinDiags);

			executeOneFFTBackward(kinDiags);

			for (int d = 0; d < nPts; d++) {
				std::complex<double> cv = kinDiags[d];
				for (int i = 0; i < nPts - d; i++)
					opMat[(i * i + (2 * d + 3) * i + d * (d + 1)) / 2] = cv;
			}

			sq_free(kinDiags); kinDiags = nullptr;

		}
	}

	void GenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>** states, int* nEigs) {
		if (nPts > 46340) {
			std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (GenDisp_PSM::findEigenStates)" << std::endl;
			throw -1;
		}
		*states = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nPts);
		calcOpMat();
		for (int i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];


		dcomplex * work = (dcomplex *)sq_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)sq_malloc(sizeof(double)*7*nPts);
		int * iwork3 = (int *)sq_malloc(sizeof(int)*5*nPts);
		double * eigs = (double *)sq_malloc(sizeof(double)*nPts);
		int * ifail = (int *)sq_malloc(sizeof(int)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		int info;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0, 0, &prec, nEigs, eigs, reinterpret_cast<dcomplex *>(*states), &nPts, work, work2, iwork3, ifail, &info);

		clearOpMat();

		nElec = *nEigs;

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

	double GenDisp_PSM::evaluateKineticEnergy(std::complex<double>* psi) {
		initializeOneFFT();

		vtls::copyArray(nPts, psi, temp1);
		//DftiComputeForward(dftiHandleKin, temp1);
		executeOneFFTForward(temp1);
		vtls::seqMulArrays(nPts, osKineticEnergy, temp1, temp2);
		for (int i = 0; i < nPts; i++)
			temp1[i] = std::conj(temp1[i]);

		return std::real(vtlsInt::rSumMul(nPts, temp1, temp2, 1.0) / vtls::getNorm(nPts, temp1, 1.0));
	}


	GenDisp_PSM_FreeElec::GenDisp_PSM_FreeElec(int nPts, double dx, double dt, double m_eff) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);//new std::complex<double>[nPts];

		double dphs = PhysCon::hbar*PhysCon::hbar / (2.0 * PhysCon::me*m_eff) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[0] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dphs * (double)(i * i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			sq_free(osKineticEnergy); osKineticEnergy = nullptr;
	}

	GenDisp_PSM_Series::GenDisp_PSM_Series(int nPts, double dx, double dt, int nPoly, double* polyCoeffs) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		osKineticEnergy[0] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dk * (double)(i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}
		vtls::polyEval(nPts, osKineticEnergy, nPoly, polyCoeffs, osKineticEnergy);

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			sq_free(osKineticEnergy); osKineticEnergy = nullptr;
	}

	GenDisp_PSM_MathExpr::GenDisp_PSM_MathExpr(int nPts, double dx, double dt, std::string expr) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
		double* ks = (double*) sq_malloc(nPts * sizeof(double));

		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		ks[0] = 0.0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
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
		if (temp1)
			sq_free(temp1);
		if (temp2)
			sq_free(temp2);
		if (temp3)
			sq_free(temp3);
		if (osKineticMask)
			sq_free(osKineticMask);
		if (norms)
			sq_free(norms);

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

	void NonUnifGenDisp_PSM::stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for collapse(2)
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//get original norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (int i = 0; i < nElec; i++)
				norms[i] = vtls::getNorm(nPts, &targ[i*nPts], dx);

		//Apply each order of kinetic exponential
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

		vtls::copyArray(nPts * nElec, targ, tempPsiCum);
		for (int o = 1; o <= expOrder; o++) {
			for (int d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
#pragma omp parallel for
				for (int i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nElec]);

				//back to real space
				//DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTBackward(&tempPsi[d * nPts * nElec]);
				//apply mask
#pragma omp parallel for
				for (int i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);

				//back to recip space
				//DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTForward(&tempPsi[d * nPts * nElec]);
				//apply rest of kinetic energy
#pragma omp parallel for
				for (int i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
			}
			//reset cumulative psi to first contribution
			vtls::copyArray(nPts * nElec, tempPsi, tempPsiCum);

			//combine all the other new psi components
			for (int d = 1; d < nDisp; d++)
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
			for (int i = 0; i < nElec; i++)
				vtls::setNorm(nPts, &targ[i * nPts], dx, norms[i]);

		//Apply half of phase contribution from potential
#pragma omp parallel for collapse(2)
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osPotentialPhase[j];
			}
		}
	}

	void NonUnifGenDisp_PSM::stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) {
		initializeAllFFT(nElec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//get original norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (int i = 0; i < nElec; i++)
				norms[i] = vtls::getNorm(nPts, &targ[i * nPts], dx);

		//Apply each order of kinetic exponential
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);
		vtls::copyArray(nPts * nElec, targ, tempPsiCum);
		for (int o = 1; o <= expOrder; o++) {
			for (int d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
				for (int i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
				//back to real space
				//DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTBackward(&tempPsi[d * nPts * nElec]);
				//apply mask
				for (int i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
				//back to recip space
				//DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nElec]);
				executeAllFFTForward(&tempPsi[d * nPts * nElec]);
				//apply rest of kinetic energy
				for (int i = 0; i < nElec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsi[i * nPts + d * nPts * nElec]);
			}
			//reset cumulative psi to first contribution
			vtls::copyArray(nPts * nElec, tempPsi, tempPsiCum);
			//combine all the other new psi components
			for (int d = 1; d < nDisp; d++)
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
			for (int i = 0; i < nElec; i++)
				vtls::setNorm(nPts, &targ[i * nPts], dx, norms[i]);
	}

	void NonUnifGenDisp_PSM::stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) {
		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for
		for (int i = 0; i < nElec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}
	}

	void NonUnifGenDisp_PSM::initializeAllFFT(int nElec) {
		if (firstStepAll || NonUnifGenDisp_PSM::nElec != nElec) {
			NonUnifGenDisp_PSM::nElec = nElec;
			/*DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nElec);
			DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
			DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
			DftiCommitDescriptor(dftiHandle);*/

			std::complex<double>* test = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nElec);

			mtx.lock();

			if(fftwAllForward)
				fftw_destroy_plan(fftwAllForward); fftwAllForward=NULL;
			if(fftwAllBackward)
				fftw_destroy_plan(fftwAllBackward); fftwAllBackward=NULL;

			fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;
			
			fftwAllForward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_FORWARD, MULTIELEC_FFTW_POLICY);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_BACKWARD, MULTIELEC_FFTW_POLICY);
			
			mtx.unlock();
			
			sq_free(test);

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

			firstStepAll = 0;
		}
	}

	void NonUnifGenDisp_PSM::initializeOneFFT() {
		if (firstStepOne) {
			/*DftiCreateDescriptor(&dftiHandleKin, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleKin, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleKin);*/

			if (temp1)
				sq_free(temp1); temp1 = nullptr;
			if (temp2)
				sq_free(temp2); temp2 = nullptr;
			if (temp3)
				sq_free(temp3); temp3 = nullptr;

			temp1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
			temp2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);
			temp3 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts);

			mtx.lock();

			//initialize FFTW for performance, find best algo
			if(fftwOneForward)
				fftw_destroy_plan(fftwOneForward); fftwOneForward=NULL;
			if(fftwOneBackward)
				fftw_destroy_plan(fftwOneBackward); fftwOneBackward=NULL;

			fftw_plan_with_nthreads(1);
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			fftwOneForward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1), FFTW_FORWARD, FFTW_ESTIMATE);
			fftwOneBackward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp2), FFTW_BACKWARD, FFTW_ESTIMATE);

			mtx.unlock();

			firstStepOne = 0;

			if (temp1)
				sq_free(temp1); temp1 = nullptr;
			if (temp2)
				sq_free(temp2); temp2 = nullptr;
			if (temp3)
				sq_free(temp3); temp3 = nullptr;
		}
	}

	void NonUnifGenDisp_PSM::executeAllFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void NonUnifGenDisp_PSM::executeAllFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
#pragma omp parallel for
		for(int i = 0; i < nElec; i++)
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

			for (int d = 0; d < nDisp; d++) {
				vtls::copyArray(nPts, &osKineticEnergy[d*nPts], kinDiags);
				//DftiComputeBackward(dftiHandleMat, kinDiags);
				executeOneFFTBackward(kinDiags);

				for (int dk = 0; dk < nPts; dk++) {
					std::complex<double> cv = kinDiags[dk];
					for (int i = 0; i < nPts - dk; i++)
						kinMat[(i * i + (2 * dk + 3) * i + dk * (dk + 1)) / 2] = cv;
				}

				vtls::mulTriagDiagTriag(nPts, kinMat, &osKineticMask[d * nPts], temp);
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

	void NonUnifGenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>** states, int* nEigs) {
		*states = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * nPts * nPts);
		calcOpMat();
		for (int i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];

		dcomplex * work = (dcomplex *)sq_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)sq_malloc(sizeof(double)*7*nPts);
		int * iwork3 = (int *)sq_malloc(sizeof(int)*5*nPts);
		double * eigs = (double *)sq_malloc(sizeof(double)*nPts);
		int * ifail = (int *)sq_malloc(sizeof(int)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		int info;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0, 0, &prec, nEigs, eigs, reinterpret_cast<dcomplex *>(*states), &nPts, work, work2, iwork3, ifail, &info);

		clearOpMat();

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

	double NonUnifGenDisp_PSM::evaluateKineticEnergy(std::complex<double>* psi) {
		initializeOneFFT();

		vtls::copyArray(nPts, psi, temp1);
		//DftiComputeForward(dftiHandleKin, temp1);
		executeOneFFTForward(temp1);

		std::fill_n(temp2, nPts, 0.0);
		for (int d = 0; d < nDisp; d++) {
			vtls::seqMulArrays(nPts, &osKineticEnergy[d*nPts], temp1, temp3);
			//DftiComputeBackward(dftiHandleKin, temp3);
			executeOneFFTBackward(temp3);
			vtls::seqMulArrays(nPts, &osKineticMask[d*nPts], temp3);
			//DftiComputeForward(dftiHandleKin, temp3);
			executeOneFFTForward(temp3);
			vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], temp3);
			vtls::addArrays(nPts, temp3, temp2);
		}
		for (int i = 0; i < nPts; i++)
			temp1[i] = std::conj(temp1[i]);

		return std::real(vtlsInt::rSumMul(nPts, temp1, temp2, 1.0) / vtls::getNorm(nPts, temp1, 1.0));
	}


	NonUnifGenDisp_PSM_EffMassBoundary::NonUnifGenDisp_PSM_EffMassBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, double meff_l, double meff_r, double transRate, int transPos, double edgeRate) : NonUnifGenDisp_PSM(nPts, dx, dt, 2, expOrder, forceNormalization) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*2);
		double* mask = (double*)sq_malloc(sizeof(double)*nPts*2);

		double dphs = PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me * meff_r) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[0] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dphs * (double)(i * i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}

		dphs = PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me * meff_l) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[nPts] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[nPts + i] = dphs * (double)(i * i);
			osKineticEnergy[nPts + nPts - i] = osKineticEnergy[i];
		}

		if (edgeRate != 0.0) {
			for (int i = 0; i < nPts; i++) {
				mask[i] = (1.0 / (1.0 + std::exp(-dx * transRate * (i - transPos))) + 1.0 / (1.0 + std::exp(dx * edgeRate * i))) / (1.0 + std::exp(dx * edgeRate * (i - nPts)));
				mask[i + nPts] = 1.0 - mask[i];
			}
		}
		else {
			for (int i = 0; i < nPts; i++) {
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

	NonUnifGenDisp_PSM_MathExprBoundary::NonUnifGenDisp_PSM_MathExprBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, int nDisp, std::vector<std::string> exprs, double* transRates, int* transPoss) : NonUnifGenDisp_PSM(nPts, dx, dt, nDisp, expOrder, forceNormalization) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nDisp);
		double* mask = (double*) sq_malloc(nPts * nDisp * sizeof(double));
		std::fill_n(mask, nPts * nDisp, 1.0);
		double* ks = (double*) sq_malloc(nPts * sizeof(double));

		// generate ks
		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		ks[0] = 0.0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			ks[i] = dk * (double)(i);
			ks[nPts - i] = -dk * (double)(i);
		}

		// generate osKineticEnergy
		for(int i = 0; i < nDisp; i++)
			vtls::evalMathExpr(nPts, "k", ks, exprs[i], &osKineticEnergy[i*nPts]);

		// generate sigmoid masks
		for (int i = 0; i < nDisp; i++) {
			for (int j = 0; j < nPts; j++) {
				for(int k = 0; k < i; k++)
					mask[i * nPts + j] *= 1.0 / (1.0 + std::exp(-dx * transRates[k] * (j - transPoss[k])));
				for (int k = i; k < nDisp - 1; k++) 
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


	void KineticOperator_FDM::projectHistory(std::complex<double>* psi, std::complex<double>* phsL, std::complex<double>* phsR, double* v, int nElec) {
		std::complex<double>* bcwfs = ( std::complex<double>* )sq_malloc(sizeof(std::complex<double>) * nElec);

		cblas_zcopy(nElec, &psi[0], nPts, bcwfs, 1);
		lbc->fillHistory(bcwfs, phsL, v[0]);

		cblas_zcopy(nElec, &psi[nPts-1], nPts, bcwfs, 1);
		rbc->fillHistory(bcwfs, phsR, v[nPts-1]);

		sq_free(bcwfs);
	}

	CrankNicolson::CrankNicolson(int nPts, double dx, double dt, double m_eff, FDBCs::BoundaryCondition* leftBC, FDBCs::BoundaryCondition* rightBC) :
		KineticOperator_FDM(nPts, leftBC, rightBC), dx(dx), dt(dt), m_eff(m_eff) {
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

			potmul = 0.5*PhysCon::im*dt/PhysCon::hbar;
	}

	void CrankNicolson::_step(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec, int isVirtual) {
		if(bct1 == nullptr)
			bct1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
		if(bct2 == nullptr)
			bct2 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
		if(lbct == nullptr)
			lbct = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
		if(rbct == nullptr)
			rbct = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);

		//prepare left BC
		cblas_zcopy(nElec, psi0, nPts, bct1, 1); // map first element of all wavefunctions to bct1
		cblas_zcopy(nElec, &psi0[1], nPts, bct2, 1); // map second element of all wavefunctions to bct2
	
		lbc->prepareStep(bct1, bct2, std::real(v[0]));
		lbc->getRHS(bct1, bct2, std::real(v[0]), lbct, nElec);
		if(!isVirtual)
			lbc->finishStep(bct1, bct2, std::real(v[0]));

		//prepare right BC
		cblas_zcopy(nElec, &psi0[nPts-1], nPts, bct1, 1); // map last element of all wavefunctions to bct1
		cblas_zcopy(nElec, &psi0[nPts-2], nPts, bct2, 1); // map second to last element of all wavefunctions to bct2
		
		rbc->prepareStep(bct1, bct2, std::real(v[nPts-1]));
		rbc->getRHS(bct1, bct2, std::real(v[nPts-1]), rbct, nElec);
		if(!isVirtual)
			rbc->finishStep(bct1, bct2, std::real(v[nPts-1]));

		//prepare LHS matrix
		std::fill_n(d, nPts, lhsDiag0);
		vtls::scaMulAddArrays(nPts-2, potmul, &v[1], &d[1]); // d += potmul*v, leave BCs alone
		std::fill_n(ud, nPts-1, lhsOffDiag0);
		std::fill_n(ld, nPts-1, lhsOffDiag0);

		d[0] = lbc->getLHSEle();
		d[nPts-1] = rbc->getLHSEle();
		ud[0] = lbc->getLHSAdjEle();
		ld[nPts-2] = rbc->getLHSAdjEle();

		// evaluate RHS
		#pragma omp parallel for collapse(2)
		for(int j = 0; j < nElec; j++)
			for(int k = 0; k < nPts-1; k++)
				targ[j*nPts+k] = (-potmul*v[k]+rhsDiag0)*psi0[j*nPts+k] +
					(rhsOffDiag*psi0[j*nPts+k-1] + rhsOffDiag*psi0[j*nPts+k+1]);

		// apply RHS BC
		cblas_zcopy(nElec, lbct, 1, targ, nPts);
		cblas_zcopy(nElec, rbct, 1, &targ[nPts-1], nPts);

		//SOLVE
		int info;
		LAPACK_zgtsv(&nPts, &nElec, reinterpret_cast<dcomplex*>(ld), reinterpret_cast<dcomplex*>(d), reinterpret_cast<dcomplex*>(ud), reinterpret_cast<dcomplex*>(targ), &nPts, &info);
	
		for(int i = 0; i < nElec; i++)
			vtls::seqMulArrays(nPts, spatialDamp, &targ[i*nPts]);
	}

	void CrankNicolson::findEigenStates(double* v, double emin, double emax, std::complex<double>** states, int* nEigs){
		double* hd = (double*)sq_malloc(sizeof(double)*nPts);
		double* hod= (double*)sq_malloc(sizeof(double)*(nPts-1));
		int  nSplit;
		int* iblock = (int*)sq_malloc(sizeof(int)*nPts);
		int* isplit = (int*)sq_malloc(sizeof(int)*nPts);
		int* iwork = (int*)sq_malloc(sizeof(int)*3*nPts);
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
		int info;
		LAPACK_dstebz("V", "B", &nPts, &emin, &emax, 0, 0, &prec, hd, hod, nEigs, &nSplit, eigs, iblock, isplit, work, iwork, &info);

		double* statesTemp = (double*)sq_malloc(sizeof(double)*nPts*(*nEigs));
		int* ifail = (int*)sq_malloc(sizeof(int)*(*nEigs));

		// get eigenvectors
		LAPACK_dstein(&nPts, hd, hod, nEigs, eigs, iblock, isplit, statesTemp, &nPts, work, iwork, ifail, &info);
		// copy eigenvectors to states
		*states = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*(*nEigs));
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

	void CrankNicolson::findInhomogeneousEigenStates(double* v, double* es, std::complex<double>* states, int nElec){
		std::complex<double>* lhs_d = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
		std::complex<double>* lhs_ld= (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
		std::complex<double>* lhs_ud= (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
		std::complex<double>* rhs   = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

		double* kls = (double*)sq_malloc(sizeof(double)*nElec);
		double* krs = (double*)sq_malloc(sizeof(double)*nElec);

		std::complex<double>* phaseAdvancement = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nElec);

		for(int i = 0; i < nElec; i++){
			//std::complex<double> phase = (1.0 - 0.5*PhysCon::im*es[i]*dt/PhysCon::hbar) / (1.0 + 0.5*PhysCon::im*es[i]*dt/PhysCon::hbar);
			phaseAdvancement[i] = phaseAdvanceFromEnergy(es[i], dt);

			// define inner system
			for(int k = 1; k < nPts-1; k++)
				lhs_d[k] = phaseAdvancement[i]*lhsDiag0 - rhsDiag0 + (1.0+phaseAdvancement[i])*potmul*v[k];
			for(int k = 0; k < nPts-2; k++){
				lhs_ld[k] = phaseAdvancement[i]*lhsOffDiag0 - rhsOffDiag;
				lhs_ud[k+1] = phaseAdvancement[i]*lhsOffDiag0 - rhsOffDiag;
			}

			// get wavenumbers on either side associated with energy
			// if energy is less than potential, wavenumber is set to decay constant
			//     this is encoded with a negative value
			try{ kls[i] = wavenumberFromEnergy(es[i], v[0], dx, dt, m_eff); }
			catch(const std::exception& e) { kls[i] = -wavenumberFromEnergy(-es[i], -v[0], dx, dt, m_eff); }
			try{ krs[i] = wavenumberFromEnergy(es[i], v[nPts-1], dx, dt, m_eff); }
			catch(const std::exception& e) { krs[i] = -wavenumberFromEnergy(-es[i], -v[nPts-1], dx, dt, m_eff); }
			
			// define boundaries of system
			lhs_d[0] = lbc->getSteadyLHSEle(phaseAdvancement[i], kls[i], v[0]);
			lhs_d[nPts-1] = rbc->getSteadyLHSEle(phaseAdvancement[i], krs[i], v[nPts-1]);
			lhs_ud[0] = lbc->getSteadyLHSAdjEle(phaseAdvancement[i], kls[i], v[0]);
			lhs_ld[nPts-2] = rbc->getSteadyLHSAdjEle(phaseAdvancement[i], krs[i], v[nPts-1]);

			// get inhomogeneous matrix
			std::fill_n(rhs, nPts, 0.0);
			rhs[0] = lbc->getSteadyRHS(phaseAdvancement[i], kls[i], v[0]);
			rhs[nPts-1] = rbc->getSteadyRHS(phaseAdvancement[i], krs[i], v[nPts-1]);

			// check if the sytem is inhomogeneous
			if(std::abs(rhs[0]) < 1e-10 && std::abs(rhs[nPts-1]) < 1e-10)
				throw std::runtime_error("System must be inhomogeneous to use findInhomogeneousEigenStates");
			//SOLVE
			int info, one=1;
			LAPACK_zgtsv(&nPts, &one, reinterpret_cast<dcomplex*>(lhs_ld), reinterpret_cast<dcomplex*>(lhs_d), reinterpret_cast<dcomplex*>(lhs_ud), reinterpret_cast<dcomplex*>(rhs), &nPts, &info);
		
			vtls::copyArray(nPts, rhs, &states[i*nPts]);
		}

		projectHistory(states, phaseAdvancement, phaseAdvancement, v, nElec);

		sq_free(lhs_d);
		sq_free(lhs_ld);
		sq_free(lhs_ud);
		sq_free(rhs);
		sq_free(kls);
		sq_free(krs);
		sq_free(phaseAdvancement);
	}

	double CrankNicolson::evaluateKineticEnergy(std::complex<double>* psi){
		std::complex<double>* temp = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-2));
		vtls::scaMulArray(nPts-2, PhysCon::hbar*PhysCon::hbar/(PhysCon::me*m_eff*dx*dx), &psi[1], temp);
		vtls::scaMulAddArrays(nPts-3, -0.5*PhysCon::hbar*PhysCon::hbar/(PhysCon::me*m_eff*dx*dx), &psi[2], &temp[0]);
		vtls::scaMulAddArrays(nPts-3, -0.5*PhysCon::hbar*PhysCon::hbar/(PhysCon::me*m_eff*dx*dx), &psi[1], &temp[1]);

		sq_free(temp);

		return std::real(vtlsInt::rSumMulConj(nPts-2, &psi[1], temp, 1.0) / vtls::getNorm(nPts-2, &psi[1], 1.0));
	}

}