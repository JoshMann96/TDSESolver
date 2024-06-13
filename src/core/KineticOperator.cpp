#include "KineticOperator.h"
#include "fftw3.h"
#include <omp.h>
#include <fftw3.h>
#include <lapack.h>

#define MULTIELEC_FFTW_POLICY FFTW_PATIENT

namespace KineticOperators {

	std::mutex fftw_plan_mutex;

	GenDisp_PSM::~GenDisp_PSM(){
		if (osKineticPhase)
			fftw_free(osKineticPhase);
		if (osPotentialPhase)
			fftw_free(osPotentialPhase);
		if (opMat)
			fftw_free(opMat);
		if (osKineticEnergy)
			fftw_free(osKineticEnergy);
		if (temp1)
			fftw_free(temp1);
		if (temp2)
			fftw_free(temp2);

		fftw_plan_mutex.lock();
		if(fftwOneForward)
			fftw_destroy_plan(fftwOneForward);
		if(fftwOneBackward)
			fftw_destroy_plan(fftwOneBackward);
		if(fftwAllForward)
			fftw_destroy_plan(fftwAllForward);
		if(fftwAllBackward)
			fftw_destroy_plan(fftwAllBackward);
		fftw_plan_mutex.unlock();
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

			
			fftw_plan_mutex.lock();

			//initialize FFTW for performance, find best algo
			std::complex<double>* test = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nElec);
			if(fftwAllForward)
				fftw_destroy_plan(fftwAllForward); fftwAllForward=NULL;
			if(fftwAllBackward)
				fftw_destroy_plan(fftwAllBackward);fftwAllBackward=NULL;
			
			fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			fftwAllForward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_FORWARD, MULTIELEC_FFTW_POLICY);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_BACKWARD, MULTIELEC_FFTW_POLICY);
			fftw_free(test);

			fftw_plan_mutex.unlock();

			if (firstStepAll) {
				if (osPotentialPhase)
					fftw_free(osPotentialPhase); osPotentialPhase = nullptr;
				if (osKineticPhase)
					fftw_free(osKineticPhase); osKineticPhase = nullptr;

				//initialize phase multipliers
				osPotentialPhase = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
				osKineticPhase = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
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

			fftw_plan_mutex.lock();

			if (temp1)
				fftw_free(temp1); temp1 = nullptr;
			if (temp2)
				fftw_free(temp2); temp2 = nullptr;

			temp1 = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
			temp2 = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);

			if(fftwOneForward)
				fftw_destroy_plan(fftwOneForward); fftwOneForward=NULL;
			if(fftwOneBackward)
				fftw_destroy_plan(fftwOneBackward); fftwOneBackward=NULL;

			//fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			fftwOneForward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1), FFTW_FORWARD, FFTW_ESTIMATE);
			fftwOneBackward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2), FFTW_BACKWARD, FFTW_ESTIMATE);

			fftw_plan_mutex.unlock();

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
			if (nPts > 46340) {
				std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (GenDisp_PSM::calcOpMat)" << std::endl;
				throw -1;
			}

			initializeOneFFT();

			if(opMat)
				fftw_free(opMat); opMat = nullptr;

			opMat = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*(nPts*(nPts+1))/2);//new std::complex<double>[(nPts * (nPts+1))/2];

			std::complex<double>* kinDiags = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*nPts);//new std::complex<double>[nPts];

			vtls::copyArray(nPts, osKineticEnergy, kinDiags);
			//DftiComputeBackward(dftiHandleMat, kinDiags);

			executeOneFFTBackward(kinDiags);

			//good for row major
			/*for (int i = 0; i < nPts; i++)
				vtls::copyArray(nPts - i, kinDiags, &opMat[(i * (2 * nPts + 1 - i)) / 2]);*/
			for (int d = 0; d < nPts; d++) {
				std::complex<double> cv = kinDiags[d];
				for (int i = 0; i < nPts - d; i++)
					opMat[(i * i + (2 * d + 3) * i + d * (d + 1)) / 2] = cv;
			}

			fftw_free(kinDiags); kinDiags = nullptr; //delete[] kinDiags;

		}
	}

	void GenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>* states, int* nEigs) {
		if (nPts > 46340) {
			std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (GenDisp_PSM::findEigenStates)" << std::endl;
			throw -1;
		}

		calcOpMat();
		for (int i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];


		dcomplex * work = (dcomplex *)fftw_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)fftw_malloc(sizeof(double)*7*nPts);
		int * iwork3 = (int *)fftw_malloc(sizeof(int)*5*nPts);
		double * eigs = (double *)fftw_malloc(sizeof(double)*nPts);
		int * ifail = (int *)fftw_malloc(sizeof(int)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		int info;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0, 0, &prec, nEigs, eigs, reinterpret_cast<dcomplex *>(states), &nPts, work, work2, iwork3, ifail, &info);

		clearOpMat();

		nElec = nEigs[0];

		if (work)
			fftw_free(work); work = nullptr;
		if (work2)
			fftw_free(work2); work2 = nullptr;
		if (iwork3)
			fftw_free(iwork3); iwork3 = nullptr;
		if (eigs)
			fftw_free(eigs); eigs = nullptr;
		if (ifail)
			fftw_free(ifail); ifail = nullptr;
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
		std::complex<double>* osKineticEnergy = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*nPts);//new std::complex<double>[nPts];

		double dphs = PhysCon::hbar*PhysCon::hbar / (2.0 * PhysCon::me*m_eff) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[0] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dphs * (double)(i * i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			fftw_free(osKineticEnergy); osKineticEnergy = nullptr;
	}

	GenDisp_PSM_Series::GenDisp_PSM_Series(int nPts, double dx, double dt, int nPoly, double* polyCoeffs) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*nPts);

		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		osKineticEnergy[0] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dk * (double)(i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}
		vtls::polyEval(nPts, osKineticEnergy, nPoly, polyCoeffs, osKineticEnergy);

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			fftw_free(osKineticEnergy); osKineticEnergy = nullptr;
	}

	GenDisp_PSM_MathExpr::GenDisp_PSM_MathExpr(int nPts, double dx, double dt, std::string expr) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*nPts);
		double* ks = new double[nPts];

		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		ks[0] = 0.0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			ks[i] = dk * (double)(i);
			ks[nPts - i] = -dk * (double)(i);
		}
		vtls::evalMathExpr(nPts, "k", ks, expr, osKineticEnergy);

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			fftw_free(osKineticEnergy); osKineticEnergy = nullptr;
		if (ks)
			delete[] ks; ks = nullptr;
	}

	NonUnifGenDisp_PSM::~NonUnifGenDisp_PSM(){
		if (osKineticEnergy)
			fftw_free(osKineticEnergy);
		if (osPotentialPhase)
			fftw_free(osPotentialPhase);
		if (opMat)
			fftw_free(opMat);
		if (tempPsi)
			fftw_free(tempPsi);
		if (tempPsiCum)
			fftw_free(tempPsiCum);
		if (temp1)
			fftw_free(temp1);
		if (temp2)
			fftw_free(temp2);
		if (temp3)
			fftw_free(temp3);
		if (osKineticMask)
			fftw_free(osKineticMask);
		if (norms)
			delete[] norms;

		fftw_plan_mutex.lock();
		if(fftwOneForward)
			fftw_destroy_plan(fftwOneForward);
		if(fftwOneBackward)
			fftw_destroy_plan(fftwOneBackward);
		if(fftwAllForward)
			fftw_destroy_plan(fftwAllForward);
		if(fftwAllBackward)
			fftw_destroy_plan(fftwAllBackward);
		fftw_plan_mutex.unlock();
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
			
			fftw_plan_mutex.lock();

			std::complex<double>* test = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nElec);
			if(fftwAllForward)
				fftw_destroy_plan(fftwAllForward); fftwAllForward=NULL;
			if(fftwAllBackward)
				fftw_destroy_plan(fftwAllBackward); fftwAllBackward=NULL;

			fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;
			
			fftwAllForward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_FORWARD, MULTIELEC_FFTW_POLICY);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts, nElec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_BACKWARD, MULTIELEC_FFTW_POLICY);
			fftw_free(test);
			
			fftw_plan_mutex.unlock();

			if (osPotentialPhase)
				fftw_free(osPotentialPhase); osPotentialPhase = nullptr;
			if (tempPsi)
				fftw_free(tempPsi); tempPsi = nullptr;
			if (tempPsiCum)
				fftw_free(tempPsiCum); tempPsiCum = nullptr;
			if (norms)
				delete[] norms; norms = nullptr;

			osPotentialPhase = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
			tempPsi = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nElec * nDisp);
			tempPsiCum = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nElec);
			norms = new double[nElec];

			firstStepAll = 0;
		}
	}

	void NonUnifGenDisp_PSM::initializeOneFFT() {
		if (firstStepOne) {
			/*DftiCreateDescriptor(&dftiHandleKin, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleKin, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleKin);*/

			fftw_plan_mutex.lock();

			if (temp1)
				fftw_free(temp1); temp1 = nullptr;
			if (temp2)
				fftw_free(temp2); temp2 = nullptr;
			if (temp3)
				fftw_free(temp3); temp3 = nullptr;

			temp1 = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
			temp2 = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
			temp3 = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);

			//initialize FFTW for performance, find best algo
			if(fftwOneForward)
				fftw_destroy_plan(fftwOneForward); fftwOneForward=NULL;
			if(fftwOneBackward)
				fftw_destroy_plan(fftwOneBackward); fftwOneBackward=NULL;

			//fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			fftwOneForward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1), FFTW_FORWARD, FFTW_ESTIMATE);
			fftwOneBackward = fftw_plan_dft(1, &nPts, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp2), FFTW_BACKWARD, FFTW_ESTIMATE);

			fftw_plan_mutex.unlock();

			firstStepOne = 0;

			if (temp1)
				fftw_free(temp1); temp1 = nullptr;
			if (temp2)
				fftw_free(temp2); temp2 = nullptr;
			if (temp3)
				fftw_free(temp3); temp3 = nullptr;
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
				fftw_free(opMat); opMat = nullptr;
			opMat = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1)) / 2);
			std::fill_n(opMat, (nPts * (nPts + 1)) / 2, 0.0);

			std::complex<double>* kinDiags = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
			std::complex<double>* kinMat = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1)) / 2);
			std::complex<double>* temp = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1)) / 2);

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
				fftw_free(kinDiags); kinDiags = nullptr;
			if (kinMat)
				fftw_free(kinMat); kinMat = nullptr;
			if (temp)
				fftw_free(temp); temp = nullptr;

		}
	}

	void NonUnifGenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>* states, int* nEigs) {
		calcOpMat();
		for (int i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];

		dcomplex * work = (dcomplex *)fftw_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)fftw_malloc(sizeof(double)*7*nPts);
		int * iwork3 = (int *)fftw_malloc(sizeof(int)*5*nPts);
		double * eigs = (double *)fftw_malloc(sizeof(double)*nPts);
		int * ifail = (int *)fftw_malloc(sizeof(int)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		int info;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0, 0, &prec, nEigs, eigs, reinterpret_cast<dcomplex *>(states), &nPts, work, work2, iwork3, ifail, &info);

		clearOpMat();

		nElec = nEigs[0];

		if (work)
			fftw_free(work); work = nullptr;
		if (work2)
			fftw_free(work2); work2 = nullptr;
		if (iwork3)
			fftw_free(iwork3); iwork3 = nullptr;
		if (eigs)
			fftw_free(eigs); eigs = nullptr;
		if (ifail)
			fftw_free(ifail); ifail = nullptr;
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
		std::complex<double>* osKineticEnergy = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*nPts*2);
		double* mask = (double*)fftw_malloc(sizeof(double)*nPts*2);

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
			fftw_free(osKineticEnergy); osKineticEnergy = nullptr;
		if (mask)
			fftw_free(mask); mask = nullptr;
	}

	NonUnifGenDisp_PSM_MathExprBoundary::NonUnifGenDisp_PSM_MathExprBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, int nDisp, std::vector<std::string> exprs, double* transRates, int* transPoss) : NonUnifGenDisp_PSM(nPts, dx, dt, nDisp, expOrder, forceNormalization) {
		std::complex<double>* osKineticEnergy = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*nPts*nDisp);
		double* mask = new double[nPts * nDisp];
		std::fill_n(mask, nPts * nDisp, 1.0);
		double* ks = new double[nPts];

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
			fftw_free(osKineticEnergy); osKineticEnergy = nullptr;
		if (mask)
			delete[] mask; mask = nullptr;
	}

}