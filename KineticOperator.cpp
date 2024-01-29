#include "KineticOperator.h"

namespace KineticOperators {

	void GenDisp_PSM::stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeAllFFT(nelec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//Perform FFT, apply full momentum-space phase contribution, invert FFT
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);

		//Apply half of phase contribution from potential, apply DFT normalization
#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osPotentialPhase[j];
			}
		}

	}

	void GenDisp_PSM::stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeAllFFT(nelec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//Perform FFT, apply full momentum-space phase contribution, invert FFT
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);
	}

	void GenDisp_PSM::stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * spatialDamp[i];

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}
	}

	void GenDisp_PSM::initializeAllFFT(int nelec) {
		if (firstStepAll || GenDisp_PSM::nelec != nelec) {
			GenDisp_PSM::nelec = nelec;
			/*DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nelec);
			DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
			DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
			DftiCommitDescriptor(dftiHandle);*/

			//initialize FFTW for performance, find best algo
			std::complex<double>* test = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec);
			if(fftwAllForward)
				fftw_destroy_plan(fftwAllForward); fftwAllForward=NULL;
			if(fftwAllBackward)
				fftw_destroy_plan(fftwAllBackward);fftwAllBackward=NULL;
			
			fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

			fftwAllForward = fftw_plan_many_dft(1, &nPts, nelec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_FORWARD, FFTW_PATIENT);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts, nelec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_BACKWARD, FFTW_PATIENT);
			fftw_free(test);

			if (firstStepAll) {
				if (osPotentialPhase)
					fftw_free(osPotentialPhase); osPotentialPhase = NULL;
				if (osKineticPhase)
					fftw_free(osKineticPhase); osKineticPhase = NULL;

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

			if (temp1)
				fftw_free(temp1); temp1 = NULL;
			if (temp2)
				fftw_free(temp2); temp2 = NULL;

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

			firstStepOne = 0;
		}
	}

	void GenDisp_PSM::executeAllFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void GenDisp_PSM::executeAllFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
#pragma omp parallel for
		for(int i = 0; i < nelec; i++)
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

			initializeOneFFT();

			if(opMat)
				LAPACKE_free(opMat); opMat = NULL;

			opMat = (std::complex<double>*)LAPACKE_malloc(sizeof(std::complex<double>)*(nPts*(nPts+1l))/2l);//new std::complex<double>[(nPts * (nPts+1))/2];

			std::complex<double>* kinDiags = (std::complex<double>*)LAPACKE_malloc(sizeof(std::complex<double>)*nPts);//new std::complex<double>[nPts];

			vtls::copyArray(nPts, osKineticEnergy, kinDiags);

			executeOneFFTBackward(kinDiags);

			//good for row major
			/*for (int i = 0; i < nPts; i++)
				vtls::copyArray(nPts - i, kinDiags, &opMat[(i * (2 * nPts + 1 - i)) / 2]);*/
		#pragma omp parallel for
			for (long d = 0l; d < nPts; d++) {
				std::complex<double> cv = kinDiags[(int)d];
				for (long i = 0l; i < nPts - d; i++)
					opMat[(i * i + (2l * d + 3l) * i + d * (d + 1l)) / 2l] = cv;
			}

			LAPACKE_free(kinDiags); kinDiags = NULL; //delete[] kinDiags;

		}
	}

	void GenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>* states, int* nEigs) {

		calcOpMat();
		for (long i = 0l; i < nPts; i++)
			opMat[(i * (i + 3l)) / 2l] += v[(int)i];


		dcomplex * work = (dcomplex *)LAPACKE_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)LAPACKE_malloc(sizeof(double)*7*nPts);
		long * iwork3 = (long *)LAPACKE_malloc(sizeof(long)*5*nPts);
		double * eigs = (double *)LAPACKE_malloc(sizeof(double)*nPts);
		long * ifail = (long *)LAPACKE_malloc(sizeof(long)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		long info;

		long nEigs_l, nPts_l=nPts;

		std::cout << "BEGINNING EIGENSOLVING" << std::endl;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts_l, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0l, 0l, &prec, &nEigs_l, eigs, reinterpret_cast<dcomplex *>(states), &nPts_l, work, work2, iwork3, ifail, &info);

		std::cout << "DONE EIGENSOLVING" << std::endl;

		*nEigs = (int)nEigs_l;
		nelec = (int)nEigs_l;

		clearOpMat();

		if (work)
			LAPACKE_free(work); work = NULL;
		if (work2)
			LAPACKE_free(work2); work2 = NULL;
		if (iwork3)
			LAPACKE_free(iwork3); iwork3 = NULL;
		if (eigs)
			LAPACKE_free(eigs); eigs = NULL;
		if (ifail)
			LAPACKE_free(ifail); ifail = NULL;
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
			fftw_free(osKineticEnergy); osKineticEnergy = NULL;
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
			fftw_free(osKineticEnergy); osKineticEnergy = NULL;
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
			fftw_free(osKineticEnergy); osKineticEnergy = NULL;
		if (ks)
			delete[] ks; ks = NULL;
	}


	void NonUnifGenDisp_PSM::stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeAllFFT(nelec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for collapse(2)
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//get original norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (int i = 0; i < nelec; i++)
				norms[i] = vtls::getNorm(nPts, &targ[i*nPts], dx);

		//Apply each order of kinetic exponential
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);

		vtls::copyArray(nPts * nelec, targ, tempPsiCum);
		for (int o = 1; o <= expOrder; o++) {
			for (int d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
#pragma omp parallel for
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nelec]);

				//back to real space
				//DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nelec]);
				executeAllFFTBackward(&tempPsi[d * nPts * nelec]);
				//apply mask
#pragma omp parallel for
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nelec]);

				//back to recip space
				//DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nelec]);
				executeAllFFTForward(&tempPsi[d * nPts * nelec]);
				//apply rest of kinetic energy
#pragma omp parallel for
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
			}
			//reset cumulative psi to first contribution
			vtls::copyArray(nPts * nelec, tempPsi, tempPsiCum);

			//combine all the other new psi components
			for (int d = 1; d < nDisp; d++)
				vtls::addArrays(nPts * nelec, &tempPsi[d * nPts * nelec], tempPsiCum);

			//apply factor (becomes factorial with multiple applications)
			vtls::scaMulArray(nPts * nelec, (-PhysCon::im * dt / PhysCon::hbar) / (double)o, tempPsiCum);

			//and add contribution to result
			vtls::addArrays(nPts * nelec, tempPsiCum, targ);
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);
		//restore norm if requested, else apply DFT normalization
		if (forceNorm)
#pragma omp parallel for
			for (int i = 0; i < nelec; i++)
				vtls::setNorm(nPts, &targ[i * nPts], dx, norms[i]);

		//Apply half of phase contribution from potential
#pragma omp parallel for collapse(2)
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osPotentialPhase[j];
			}
		}
	}

	void NonUnifGenDisp_PSM::stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeAllFFT(nelec);

		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}

		//get original norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (int i = 0; i < nelec; i++)
				norms[i] = vtls::getNorm(nPts, &targ[i * nPts], dx);

		//Apply each order of kinetic exponential
		//DftiComputeForward(dftiHandle, targ);
		executeAllFFTForward(targ);
		vtls::copyArray(nPts * nelec, targ, tempPsiCum);
		for (int o = 1; o <= expOrder; o++) {
			for (int d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
				//back to real space
				//DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nelec]);
				executeAllFFTBackward(&tempPsi[d * nPts * nelec]);
				//apply mask
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
				//back to recip space
				//DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nelec]);
				executeAllFFTForward(&tempPsi[d * nPts * nelec]);
				//apply rest of kinetic energy
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
			}
			//reset cumulative psi to first contribution
			vtls::copyArray(nPts * nelec, tempPsi, tempPsiCum);
			//combine all the other new psi components
			for (int d = 1; d < nDisp; d++)
				vtls::addArrays(nPts * nelec, &tempPsi[d * nPts * nelec], tempPsiCum);
			//apply factor (becomes factorial with multiple applications)
			vtls::scaMulArray(nPts * nelec, (-PhysCon::im * dt / PhysCon::hbar) / (double)o, tempPsiCum);
			//and add contribution to result
			vtls::addArrays(nPts * nelec, tempPsiCum, targ);
		}

		//DftiComputeBackward(dftiHandle, targ);
		executeAllFFTBackward(targ);

		//restore norm if requested
		if (forceNorm)
#pragma omp parallel for
			for (int i = 0; i < nelec; i++)
				vtls::setNorm(nPts, &targ[i * nPts], dx, norms[i]);
	}

	void NonUnifGenDisp_PSM::stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
		for (int i = 0; i < nPts; i++)
			osPotentialPhase[i] = std::exp(vcnst * v[i]) * std::sqrt(spatialDamp[i]);

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
			}
		}
	}

	void NonUnifGenDisp_PSM::initializeAllFFT(int nelec) {
		if (firstStepAll || NonUnifGenDisp_PSM::nelec != nelec) {
			NonUnifGenDisp_PSM::nelec = nelec;
			/*DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nelec);
			DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
			DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
			DftiCommitDescriptor(dftiHandle);*/

			std::complex<double>* test = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec);
			if(fftwAllForward)
				fftw_destroy_plan(fftwAllForward); fftwAllForward=NULL;
			if(fftwAllBackward)
				fftw_destroy_plan(fftwAllBackward); fftwAllBackward=NULL;

			fftw_plan_with_nthreads(omp_get_max_threads());
			//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;
			
			fftwAllForward = fftw_plan_many_dft(1, &nPts, nelec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_FORWARD, FFTW_PATIENT);
			fftwAllBackward = fftw_plan_many_dft(1, &nPts, nelec, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, reinterpret_cast<fftw_complex*>(test), &nPts, 1, nPts, FFTW_BACKWARD, FFTW_PATIENT);
			fftw_free(test);

			if (osPotentialPhase)
				fftw_free(osPotentialPhase); osPotentialPhase = NULL;
			if (tempPsi)
				fftw_free(tempPsi); tempPsi = NULL;
			if (tempPsiCum)
				fftw_free(tempPsiCum); tempPsiCum = NULL;
			if (norms)
				fftw_free(norms); norms = NULL;

			osPotentialPhase = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts);
			tempPsi = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec * nDisp);
			tempPsiCum = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * nPts * nelec);
			norms = new double[nelec];

			firstStepAll = 0;
		}
	}

	void NonUnifGenDisp_PSM::initializeOneFFT() {
		if (firstStepOne) {
			/*DftiCreateDescriptor(&dftiHandleKin, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleKin, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleKin);*/

			if (temp1)
				fftw_free(temp1); temp1 = NULL;
			if (temp2)
				fftw_free(temp2); temp2 = NULL;
			if (temp3)
				fftw_free(temp3); temp3 = NULL;

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

			firstStepOne = 0;
		}
	}

	void NonUnifGenDisp_PSM::executeAllFFTForward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllForward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
	}
	
	void NonUnifGenDisp_PSM::executeAllFFTBackward(std::complex<double>* targ){
		fftw_execute_dft(fftwAllBackward, reinterpret_cast<fftw_complex*>(targ), reinterpret_cast<fftw_complex*>(targ));
#pragma omp parallel for
		for(int i = 0; i < nelec; i++)
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

			initializeOneFFT();
			if(opMat)
				LAPACKE_free(opMat); opMat = NULL;
			opMat = (std::complex<double>*) LAPACKE_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1l)) / 2l);
			std::fill_n(opMat, (nPts * (nPts + 1l)) / 2l, 0.0);

			std::complex<double>* kinDiags = (std::complex<double>*) LAPACKE_malloc(sizeof(std::complex<double>) * nPts);
			std::complex<double>* kinMat = (std::complex<double>*) LAPACKE_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1l)) / 2l);
			std::complex<double>* temp = (std::complex<double>*) LAPACKE_malloc(sizeof(std::complex<double>) * (nPts * (nPts + 1l)) / 2l);

			for (int d = 0; d < nDisp; d++) {
				vtls::copyArray(nPts, &osKineticEnergy[d*nPts], kinDiags);
				//DftiComputeBackward(dftiHandleMat, kinDiags);
				executeOneFFTBackward(kinDiags);

			#pragma omp parallel for
				for (long dk = 0; dk < nPts; dk++) {
					std::complex<double> cv = kinDiags[dk];
					for (long i = 0; i < nPts - dk; i++)
						kinMat[(i * i + (2l * dk + 3l) * i + dk * (dk + 1l)) / 2l] = cv;
				}

				vtls::mulTriagDiagTriag(nPts, kinMat, &osKineticMask[(long)d * nPts], temp);
				vtls::addArrays((nPts * (nPts + 1l)) / 2l, temp, opMat);
			}

			if (kinDiags)
				LAPACKE_free(kinDiags); kinDiags = NULL;
			if (kinMat)
				LAPACKE_free(kinMat); kinMat = NULL;
			if (temp)
				LAPACKE_free(temp); temp = NULL;

		}
	}

	void NonUnifGenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>* states, int* nEigs) {

		calcOpMat();
		for (long i = 0l; i < nPts; i++)
			opMat[(i * (i + 3l)) / 2l] += v[i];

		dcomplex * work = (dcomplex *)LAPACKE_malloc(sizeof(dcomplex)*2*nPts);
		double * work2 = (double *)LAPACKE_malloc(sizeof(double)*7*nPts);
		long * iwork3 = (long *)LAPACKE_malloc(sizeof(long)*5*nPts);
		double * eigs = (double *)LAPACKE_malloc(sizeof(double)*nPts);
	    long * ifail = (long *)LAPACKE_malloc(sizeof(long)*nPts);

		char cV = 'V', cU = 'U', cS = 'S';

		double prec = LAPACK_dlamch(&cS);//(2 * dlamch_(&cS));
		long info;

		long nEigs_l, nPts_l=nPts;

		LAPACK_zhpevx(&cV, &cV, &cU, &nPts_l, reinterpret_cast<dcomplex *>(opMat), &emin, &emax, 0l, 0l, &prec, &nEigs_l, eigs, reinterpret_cast<dcomplex *>(states), &nPts_l, work, work2, iwork3, ifail, &info);

		clearOpMat();

		*nEigs = (int)nEigs_l;
		nelec = (int)nEigs_l;

		if (work)
			LAPACKE_free(work); work = NULL;
		if (work2)
			LAPACKE_free(work2); work2 = NULL;
		if (iwork3)
			LAPACKE_free(iwork3); iwork3 = NULL;
		if (eigs)
			LAPACKE_free(eigs); eigs = NULL;
		if (ifail)
			LAPACKE_free(ifail); ifail = NULL;
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
			fftw_free(osKineticEnergy); osKineticEnergy = NULL;
		if (mask)
			fftw_free(mask); mask = NULL;
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
			fftw_free(osKineticEnergy); osKineticEnergy = NULL;
		if (mask)
			delete[] mask; mask = NULL;
	}

}