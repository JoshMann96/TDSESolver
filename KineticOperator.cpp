#include "stdafx.h"

namespace KineticOperators {

	void GenDisp_PSM::stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) {
		initializeFFT();

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
		initializeFFT();

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

	void GenDisp_PSM::initializeFFT() {
		if (firstStep) {
			DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
			DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nelec);
			DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
			DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
			//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
			DftiCommitDescriptor(dftiHandle);

			if (osPotentialPhase)
				delete[] osPotentialPhase; osPotentialPhase = NULL;
			if (osKineticPhase)
				delete[] osKineticPhase; osKineticPhase = NULL;

			osPotentialPhase = new std::complex<double>[nPts];
			osKineticPhase = new std::complex<double>[nPts];
			for (int i = 0; i < nPts; i++)
				osKineticPhase[i] = std::exp(-PhysCon::im * dt / PhysCon::hbar * osKineticEnergy[i]);

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

	void GenDisp_PSM::clearOpMat() {
		if (opMat)
			delete[] opMat; opMat = NULL;
		needMat = 1;
	}

	void GenDisp_PSM::findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs) {
		std::complex<double>* states = new std::complex<double>[nPts * nPts];
		double* eigs = new double[nPts];
		int* ifail = new int[nPts];

		calcOpMat();
		for (int i = 0; i < nPts; i++)
			opMat[(i * (i + 3)) / 2] += v[i];
 
		LAPACKE_zhpevx(LAPACK_COL_MAJOR, 'V', 'V', 'U', nPts, opMat, emin, emax, 0, 0, 2 * LAPACKE_dlamch('S'), nEigs, eigs, states, nPts, ifail);

		clearOpMat();

		nelec = nEigs[0];

		vtlsPrnt::printArray(nelec, eigs);

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
}