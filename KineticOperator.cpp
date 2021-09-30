#include "stdafx.h"

namespace KineticOperators {

	void GenDisp_PSM::stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeFFT(nelec);

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
		DftiComputeForward(dftiHandle, targ);

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		DftiComputeBackward(dftiHandle, targ);

		//Apply half of phase contribution from potential
#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osPotentialPhase[j];
			}
		}
	}

	void GenDisp_PSM::stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeFFT(nelec);

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
		DftiComputeForward(dftiHandle, targ);

#pragma omp parallel for
		for (int i = 0; i < nelec; i++) {
			for (int j = 0; j < nPts; j++) {
				targ[i * nPts + j] *= osKineticPhase[j];
			}
		}

		DftiComputeBackward(dftiHandle, targ);
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

	void GenDisp_PSM::initializeFFT(int nelec) {
		if (firstStep || GenDisp_PSM::nelec != nelec) {
			GenDisp_PSM::nelec = nelec;
			DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nelec);
			DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
			DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
			DftiCommitDescriptor(dftiHandle);

			if (firstStep) {
				if (osPotentialPhase)
					delete[] osPotentialPhase; osPotentialPhase = NULL;
				if (osKineticPhase)
					delete[] osKineticPhase; osKineticPhase = NULL;

				osPotentialPhase = new std::complex<double>[nPts];
				osKineticPhase = new std::complex<double>[nPts];
				for (int i = 0; i < nPts; i++)
					osKineticPhase[i] = std::exp(-PhysCon::im * dt / PhysCon::hbar * osKineticEnergy[i]);
			}

			firstStep = 0;
		}
	}

	void GenDisp_PSM::initializeMatFFT() {
		if (firstStepMat) {
			DftiCreateDescriptor(&dftiHandleMat, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleMat, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleMat);

			firstStepMat = 0;
		}
	}

	void GenDisp_PSM::initializeKinFFT() {
		if (firstStepKin) {
			DftiCreateDescriptor(&dftiHandleKin, DFTI_COMPLEX, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleKin, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleKin);

			if (temp1)
				delete[] temp1; temp1 = NULL;
			if (temp2)
				delete[] temp2; temp2 = NULL;

			temp1 = new std::complex<double>[nPts];
			temp2 = new std::complex<double>[nPts];

			firstStepKin = 0;
		}
	}

	void GenDisp_PSM::calcOpMat() {
		if(needMat){
			needMat = 0;
			if (nPts > 46340) {
				std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (GenDisp_PSM::calcOpMat)" << std::endl;
				throw NULL;
			}

			initializeMatFFT();

			opMat = new std::complex<double>[(nPts * (nPts+1))/2];
			std::complex<double>* kinDiags = new std::complex<double>[nPts];
			vtls::copyArray(nPts, osKineticEnergy, kinDiags);
			DftiComputeBackward(dftiHandleMat, kinDiags);
			//good for row major
			/*for (int i = 0; i < nPts; i++)
				vtls::copyArray(nPts - i, kinDiags, &opMat[(i * (2 * nPts + 1 - i)) / 2]);*/
			for (int d = 0; d < nPts; d++) {
				std::complex<double> cv = kinDiags[d];
				for (int i = 0; i < nPts - d; i++)
					opMat[(i * i + (2 * d + 3) * i + d * (d + 1)) / 2] = cv;
			}
			delete[] kinDiags;
		}
	}

	void GenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs) {
		if (nPts > 46340) {
			std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (GenDisp_PSM::findEigenStates)" << std::endl;
			throw NULL;
		}
		std::complex<double>* states = new std::complex<double>[nPts * nPts];
		double* eigs = new double[nPts];
		int* ifail = new int[nPts];

		calcOpMat();
		for (int i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];

		LAPACKE_zhpevx(LAPACK_COL_MAJOR, 'V', 'V', 'U', nPts, opMat, emin, emax, 0, 0, 2 * LAPACKE_dlamch('S'), nEigs, eigs, states, nPts, ifail);

		clearOpMat();

		nelec = nEigs[0];

		psi[0] = new std::complex<double>[nPts * nelec];

		vtls::copyArray(nPts * nelec, states, psi[0]);
		for (int i = 0; i < nelec; i++)
			vtls::normalizeSqrNorm(nPts, &psi[0][i * nPts], dx);

		if (eigs)
			delete[] eigs; eigs = NULL;
		if (states)
			delete[] states; states = NULL;
		if (ifail)
			delete[] ifail; ifail = NULL;

		std::cout << "GenDisp_PSM created " << nelec << " electrons." << std::endl;
	}

	double GenDisp_PSM::evaluateKineticEnergy(std::complex<double>* psi) {
		initializeMatFFT();

		DftiComputeForward(dftiHandleKin, psi, temp1);
		vtls::seqMulArrays(nPts, osKineticEnergy, temp1, temp2);
		for (int i = 0; i < nPts; i++)
			temp1[i] = std::conj(temp1[i]);

		return std::real(vtlsInt::rSumMul(nPts, temp1, temp2, 1.0) / vtls::getNorm(nPts, temp1, 1.0));
	}


	GenDisp_PSM_FreeElec::GenDisp_PSM_FreeElec(int nPts, double dx, double dt, double m_eff) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = new std::complex<double>[nPts];

		double dphs = PhysCon::hbar*PhysCon::hbar / (2.0 * PhysCon::me*m_eff) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2);
		osKineticEnergy[0] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dphs * (double)(i * i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			delete[] osKineticEnergy; osKineticEnergy = NULL;
	}

	GenDisp_PSM_Series::GenDisp_PSM_Series(int nPts, double dx, double dt, int nPoly, double* polyCoeffs) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = new std::complex<double>[nPts];

		double dk = 2.0 * PhysCon::pi / (nPts * dx);
		osKineticEnergy[0] = 0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticEnergy[i] = dk * (double)(i);
			osKineticEnergy[nPts - i] = osKineticEnergy[i];
		}
		vtls::polyEval(nPts, osKineticEnergy, nPoly, polyCoeffs, osKineticEnergy);

		GenDisp_PSM::set_osKineticEnergy(osKineticEnergy);
		if (osKineticEnergy)
			delete[] osKineticEnergy; osKineticEnergy = NULL;
	}

	GenDisp_PSM_MathExpr::GenDisp_PSM_MathExpr(int nPts, double dx, double dt, std::string expr) : GenDisp_PSM(nPts, dx, dt) {
		std::complex<double>* osKineticEnergy = new std::complex<double>[nPts];
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
			delete[] osKineticEnergy; osKineticEnergy = NULL;
		if (ks)
			delete[] ks; ks = NULL;
	}


	void NonUnifGenDisp_PSM::stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeFFT(nelec);

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
		DftiComputeForward(dftiHandle, targ);

		vtls::copyArray(nPts * nelec, targ, tempPsiCum);
		for (int o = 1; o <= expOrder; o++) {
			for (int d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
#pragma omp parallel for
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
				//back to real space
				DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nelec]);
				//apply mask
#pragma omp parallel for
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
				//back to recip space
				DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nelec]);
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

		DftiComputeBackward(dftiHandle, targ);

		//restore norm if requested
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
		initializeFFT(nelec);

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
		DftiComputeForward(dftiHandle, targ);
		vtls::copyArray(nPts * nelec, targ, tempPsiCum);
		for (int o = 1; o <= expOrder; o++) {
			for (int d = 0; d < nDisp; d++) {
				//apply half current dispersion kinetic energy
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticEnergy[d * nPts], &tempPsiCum[i * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
				//back to real space
				DftiComputeBackward(dftiHandle, &tempPsi[d * nPts * nelec]);
				//apply mask
				for (int i = 0; i < nelec; i++)
					vtls::seqMulArrays(nPts, &osKineticMask[d * nPts], &tempPsi[i * nPts + d * nPts * nelec]);
				//back to recip space
				DftiComputeForward(dftiHandle, &tempPsi[d * nPts * nelec]);
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

		DftiComputeBackward(dftiHandle, targ);

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

	void NonUnifGenDisp_PSM::initializeFFT(int nelec) {
		if (firstStep || NonUnifGenDisp_PSM::nelec != nelec) {
			NonUnifGenDisp_PSM::nelec = nelec;
			DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nelec);
			DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
			DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
			DftiCommitDescriptor(dftiHandle);

			if (osPotentialPhase)
				delete[] osPotentialPhase; osPotentialPhase = NULL;
			if (tempPsi)
				delete[] tempPsi; tempPsi = NULL;
			if (tempPsiCum)
				delete[] tempPsiCum; tempPsiCum = NULL;
			if (norms)
				delete[] norms; norms = NULL;

			osPotentialPhase = new std::complex<double>[nPts];
			tempPsi = new std::complex<double>[nPts * nelec * nDisp];
			tempPsiCum = new std::complex<double>[nPts * nelec];
			norms = new double[nelec];

			firstStep = 0;
		}
	}

	void NonUnifGenDisp_PSM::initializeMatFFT() {
		if (firstStepMat) {
			DftiCreateDescriptor(&dftiHandleMat, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleMat, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleMat);

			firstStepMat = 0;
		}
	}

	void NonUnifGenDisp_PSM::calcOpMat() {
		if (needMat) {
			needMat = 0;
			if (nPts > 46340) {
				std::cout << "Long datatype is required for grids of size nPts>46340. Rewrite this code (NonUnifGenDisp_PSM::calcOpMat)" << std::endl;
				throw NULL;
			}

			initializeMatFFT();

			opMat = new std::complex<double>[(nPts * (nPts + 1)) / 2];
			std::fill_n(opMat, (nPts * (nPts + 1)) / 2, 0.0);
			std::complex<double>* kinDiags = new std::complex<double>[nPts];
			std::complex<double>* kinMat = new std::complex<double>[(nPts * (nPts + 1)) / 2];
			std::complex<double>* temp = new std::complex<double>[(nPts * (nPts + 1)) / 2];
			for (int d = 0; d < nDisp; d++) {
				vtls::copyArray(nPts, &osKineticEnergy[d*nPts], kinDiags);
				DftiComputeBackward(dftiHandleMat, kinDiags);

				for (int dk = 0; dk < nPts; dk++) {
					std::complex<double> cv = kinDiags[dk];
					for (int i = 0; i < nPts - dk; i++)
						kinMat[(i * i + (2 * dk + 3) * i + dk * (dk + 1)) / 2] = cv;
				}

				vtls::mulTriagDiagTriag(nPts, kinMat, &osKineticMask[d * nPts], temp);
				vtls::addArrays((nPts * (nPts + 1)) / 2, temp, opMat);
			}

			if (kinDiags)
				delete[] kinDiags; kinDiags = NULL;
			if (kinMat)
				delete[] kinMat; kinMat = NULL;
			if (temp)
				delete[] temp; temp = NULL;
		}
	}

	void NonUnifGenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs) {
		std::complex<double>* states = new std::complex<double>[nPts * nPts];
		double* eigs = new double[nPts];
		int* ifail = new int[nPts];

		calcOpMat();
		for (int i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];

		LAPACKE_zhpevx(LAPACK_COL_MAJOR, 'V', 'V', 'U', nPts, opMat, emin, emax, 0, 0, 2 * LAPACKE_dlamch('S'), nEigs, eigs, states, nPts, ifail);

		clearOpMat();

		nelec = nEigs[0];

		*psi = new std::complex<double>[nPts * nelec];

		vtls::copyArray(nPts * nelec, states, psi[0]);
		for (int i = 0; i < nelec; i++)
			vtls::normalizeSqrNorm(nPts, &psi[0][i * nPts], dx);

		if (eigs)
			delete[] eigs; eigs = NULL;
		if (states)
			delete[] states; states = NULL;
		if (ifail)
			delete[] ifail; ifail = NULL;

		std::cout << "NonUnifGenDisp_PSM created " << nelec << " electrons." << std::endl;
	}

	void NonUnifGenDisp_PSM::initializeKinFFT() {
		if (firstStepKin) {
			DftiCreateDescriptor(&dftiHandleKin, DFTI_COMPLEX, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandleKin, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			DftiCommitDescriptor(dftiHandleKin);

			if (temp1)
				delete[] temp1; temp1 = NULL;
			if (temp2)
				delete[] temp2; temp2 = NULL;
			if (temp3)
				delete[] temp3; temp3 = NULL;

			temp1 = new std::complex<double>[nPts];
			temp2 = new std::complex<double>[nPts];
			temp3 = new std::complex<double>[nPts];

			firstStepKin = 0;
		}
	}

	double NonUnifGenDisp_PSM::evaluateKineticEnergy(std::complex<double>* psi) {
		initializeMatFFT();

		DftiComputeForward(dftiHandleKin, psi, temp1);
		std::fill_n(temp2, nPts, 0.0);
		for (int d = 0; d < nDisp; d++) {
			vtls::seqMulArrays(nPts, osKineticEnergy, temp1, temp3);
			vtls::addArrays(nPts, temp3, temp2);
		}
		for (int i = 0; i < nPts; i++)
			temp1[i] = std::conj(temp1[i]);

		return std::real(vtlsInt::rSumMul(nPts, temp1, temp2, 1.0) / vtls::getNorm(nPts, temp1, 1.0));
	}


	NonUnifGenDisp_PSM_EffMassBoundary::NonUnifGenDisp_PSM_EffMassBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, double meff_l, double meff_r, double transRate, int transPos) : NonUnifGenDisp_PSM(nPts, dx, dt, 2, expOrder, forceNormalization) {
		std::complex<double>* osKineticEnergy = new std::complex<double>[nPts*2];
		double* mask = new double[nPts * 2];

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

		for (int i = 0; i < nPts; i++) {
			mask[i] = 1.0 / (1.0 + std::exp(-dx * transRate * (i - transPos)));
			mask[i + nPts] = 1.0 - mask[i];
		}

		NonUnifGenDisp_PSM::set_osKineticEnergy(osKineticEnergy, mask);
		if (osKineticEnergy)
			delete[] osKineticEnergy; osKineticEnergy = NULL;
		if (mask)
			delete[] mask; mask = NULL;
	}

	NonUnifGenDisp_PSM_MathExprBoundary::NonUnifGenDisp_PSM_MathExprBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, int nDisp, std::vector<std::string> exprs, double* transRates, int* transPoss) : NonUnifGenDisp_PSM(nPts, dx, dt, nDisp, expOrder, forceNormalization) {
		std::complex<double>* osKineticEnergy = new std::complex<double>[nPts * nDisp];
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
			delete[] osKineticEnergy; osKineticEnergy = NULL;
		if (mask)
			delete[] mask; mask = NULL;
	}

}