#include "stdafx.h"



TDSEIterator1D::TDSEIterator1D(double dx, double dt, int simSize) {
	using namespace PhysCon;
	TDSEIterator1D::dt = dt;
	TDSEIterator1D::dx = dx;
	cnAl = hbar * dt / (2 * me) / (2 * (dx*dx)) * im;
	cnBe = dt / (2 * hbar) * im;
	eulAl = hbar * dt / (2 * me * dx * dx) * im;
	eulBe = -dt / hbar * im;
	nPts = simSize;
	timePhase = new std::complex<double>[nPts];
	std::fill_n(timePhase, nPts, 1.0);
	spatialDamp = new std::complex<double>[nPts];
	std::fill_n(spatialDamp, nPts, 1.0);

	h2al = std::pow(PhysCon::hbar, 4) / (4.0*PhysCon::me*PhysCon::me*dx*dx*dx*dx);
	h2be = std::pow(PhysCon::hbar, 2) / (2.0*PhysCon::me*dx*dx);
	h2dt = 4.0*std::pow(PhysCon::me, 2)*std::pow(dx / (PhysCon::hbar*PhysCon::pi), 4)/10000.0;
	//h2dt = dt * dt;
	//numThreads = std::fmax(std::fmin(threadMax, thread::hardware_concurrency()), 2);
	//if (omp_get_nested())
	//	numThreads = omp_get_num_threads();
	//else
	//	numThreads = 1;
	//pool = new ThreadPool(numThreads);
	//std::cout << "Number of Electron Threads: " << numThreads << std::endl;
}

void TDSEIterator1D::addTimePhase(std::complex<double>* arr) {
	for (int i = 0; i < nPts; i++)
		timePhase[i] += arr[i];
}

void TDSEIterator1D::addSpatialDamp(double* arr) {
	for (int i = 0; i < nPts; i++)
		spatialDamp[i] *= arr[i];
}

void TDSEIterator1D::initializeCN() {
	dl = new std::complex<double>[nPts - 1];
	d = new std::complex<double>[nPts];
	du = new std::complex<double>[nPts - 1];
	du2 = new std::complex<double>[nPts - 2];
	ipiv = new int[nPts];

	temp = new std::complex<double>[nPts];
	dl0 = new std::complex<double>[nPts - 1];
	d0 = new std::complex<double>[nPts];
	du0 = new std::complex<double>[nPts - 1];
}

void TDSEIterator1D::initializeCNPA_OBS() {
	pt = new int[64];
	iparm = new int[64];
	perm = new int[nPts];
	x = new std::complex<double>[nPts];
	pardisoinit(pt, &mtype, iparm);
	iparm[0] = 1;
	iparm[34] = 1;
	iparm[5] = 1;
	//iparm[26] = 1;
	a = new std::complex<double>[nPts * 3 - 2];
	ia = new int[nPts + 1];
	ja = new int[nPts * 3 - 2];
	ia[0] = 0;
	for (int i = 1; i < nPts; i++) {
		ia[i] = i * 3 - 1;
	}
	ia[nPts] = nPts * 3 - 2;
	ja[0] = 0;
	ja[1] = 1;
	for (int i = 1; i < nPts - 1; i++) {
		ja[i * 3 - 1] = i - 1;
		ja[i * 3] = i;
		ja[i * 3 + 1] = i + 1;
	}
	ja[nPts * 3 - 4] = nPts - 2;
	ja[nPts * 3 - 3] = nPts - 1;

	a[1] = -cnAl * timePhase[0];
	for (int i = 1; i < nPts - 1; i++) {
		a[3 * i - 1] = -cnAl * timePhase[i];
		a[3 * i + 1] = -cnAl * timePhase[i];
	}
	a[3 * nPts - 4] = -cnAl * timePhase[nPts - 1];
}

void TDSEIterator1D::stepEuler(std::complex<double> * psi0, double * v, std::complex<double> * targ) {
	//IMPLEMENT EULER
	/*if (d2psi == NULL) d2psi = new std::complex<double>[nPts];
	vtls::secondDifference(nPts, psi0, d2psi);
	vtls::scaMulArray(nPts, eulAl, d2psi);

	vtls::seqMulArrays(nPts, v, psi0, targ);
	vtls::scaMulArray(nPts, eulAl, targ);

	vtls::addArrays(nPts, d2psi, targ);
	vtls::addArrays(nPts, psi0, targ);*/

	targ[0] = (eulBe*v[0] - 2.0 * eulAl + 1.0)*psi0[0] + eulAl * psi0[1];
	for (int i = 1; i < nPts - 1; i++)
		targ[i] = eulAl * (psi0[i - 1] + psi0[i + 1]) + (eulBe*v[i] - 2.0 * eulAl + 1.0)*psi0[i];
	targ[nPts - 1] = eulAl * psi0[nPts - 2] + (eulBe*v[nPts - 1] - 2.0 * eulAl + 1.0)*psi0[nPts - 1];
}

void TDSEIterator1D::stepEulerImT(std::complex<double> * psi0, double * v, std::complex<double> * targ) {
	//IMPLEMENT EULER IMAGINARY TIME
	std::complex<double> nAl = -PhysCon::im*eulAl;
	std::complex<double> nBe = -PhysCon::im*eulBe;
	targ[0] = (nBe*v[0] - 2.0 * nAl + 1.0)*psi0[0] + nAl * psi0[1];
	for (int i = 1; i < nPts - 1; i++)
		targ[i] = nAl * (psi0[i - 1] + psi0[i + 1]) + (nBe*v[i] - 2.0 * nAl + 1.0)*psi0[i];
	targ[nPts - 1] = nAl * psi0[nPts - 2] + (nBe*v[nPts - 1] - 2.0 * nAl + 1.0)*psi0[nPts - 1];
}

void TDSEIterator1D::stepCN(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ) {
	// INTEL MKL APPROACH
#pragma omp parallel
	{
	#pragma omp sections
		{
		#pragma omp section
			{
				std::fill_n(du, nPts - 1, -cnAl);
				vtls::seqMulArrays(nPts - 1, timePhase, du);
			}
		#pragma omp section
			{
				std::fill_n(dl, nPts - 1, -cnAl);
				vtls::seqMulArrays(nPts - 1, &timePhase[1], dl);
			}
		}
		#pragma omp for
		for (int i = 0; i < nPts; i++) {
			d[i] = 1.0 + (2.0 * cnAl * timePhase[i] + cnBe * vf[i]);
		}
		std::complex<double> nAl = timePhase[0] * cnAl;
		targ[0] = (1.0 - 2.0 * nAl - cnBe * v0[0]) * psi0[0] + nAl * psi0[1];
		#pragma omp for private(nAl)
		for (int i = 1; i < nPts - 1; i++) {
			nAl = timePhase[i] * cnAl;
			targ[i] = (1.0 - 2.0 * nAl - cnBe * v0[i]) * psi0[i] + nAl * (psi0[i - 1] + psi0[i + 1]);
		}
		nAl = timePhase[nPts - 1] * cnAl;
		targ[nPts - 1] = (1.0 - 2.0 * nAl - cnBe * v0[nPts - 1]) * psi0[nPts - 1] + nAl * psi0[nPts - 2];
	}

	LAPACKE_zgttrf(nPts, dl, d, du, du2, ipiv);
	LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl, d, du, du2, ipiv, targ, nPts);
	//~18.5s for 1,000,000 points 100 times
	//~18s for 100,000 points 1000 times
}

void TDSEIterator1D::stepOS_U2TU(std::complex<double>* __restrict psi0, double* __restrict v0, std::complex<double>* targ, int nelec) {
		if (first_stepOS_U2TU) {
		DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
		DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nelec);
		DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
		DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
		//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
		DftiCommitDescriptor(dftiHandle);

		osKineticPhase = new std::complex<double>[nPts];
		osPotentialPhase = new std::complex<double>[nPts];
		//tempD = new double[nPts];
		std::complex<double> dphs = -PhysCon::im * PhysCon::hbar / (2.0 * PhysCon::me) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2) * dt;
		osKineticPhase[0] = 1.0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticPhase[i] = std::exp(dphs * (std::complex<double>)(i * i));
			osKineticPhase[nPts - i] = osKineticPhase[i];
		}
		/*if (nPts % 2)
			osKineticPhase[nPts / 2 + 1] = std::exp(dphs * (std::complex<double>)((nPts / 2 + 1) * (nPts / 2 + 1)));*/

		first_stepOS_U2TU = 0;
	}

	std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
	for (int i = 0; i < nPts; i++)
		osPotentialPhase[i] = std::exp(vcnst * v0[i]) * spatialDamp[i];

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
	/*
	//Move over potential for out-of-place evaluation
	//cblas_dcopy(nPts, v0, 1, tempD, 1);
	for (int i = 0; i < nPts; i++)
		tempD[i] = v0[i];
	//Multiply dt/hbar/2
	//cblas_dscal(nPts, -dt / PhysCon::hbar / 2.0, tempD, 1);
	for (int i = 0; i < nPts; i++)
		tempD[i] *= -dt / PhysCon::hbar / 2.0;
	//Get complex exponential
	//vzCIS(nPts, tempD, osPotentialPhase);
	for (int i = 0; i < nPts; i++)
		osPotentialPhase[i] = std::cos(tempD[i]) + PhysCon::im * std::sin(tempD[i]);
	//Multiply by spatial damper
	//vzMul(nPts, spatialDamp, osPotentialPhase, osPotentialPhase);
	for (int i = 0; i < nPts; i++)
		osPotentialPhase[i] *= spatialDamp[i];
	//Apply half of phase contribution from potential

	for (int i = 0; i < nelec; i++)
		for (int j = 0; j < nPts; j++)
			targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
		//vzMul(nPts, osPotentialPhase, &psi0[i * nPts], &targ[i * nPts]);

	//Perform FFT, apply full momentum-space phase contribution, invert FFT
	DftiComputeForward(dftiHandle, targ);

	for (int i = 0; i < nelec; i++)
		for (int j = 0; j < nPts; j++)
			targ[i * nPts + j] *= osKineticPhase[j];
		//vzMul(nPts, osKineticPhase, &targ[i * nPts], &targ[i * nPts]);
	DftiComputeBackward(dftiHandle, targ);

	//Apply half of phase contribution from potential

	for (int i = 0; i < nelec; i++)
		for (int j = 0; j < nPts; j++)
			targ[i * nPts + j] *= osPotentialPhase[j];
		//vzMul(nPts, osPotentialPhase, &targ[i * nPts], &targ[i * nPts]);

		*/
}

void TDSEIterator1D::stepOS_UW2T(std::complex<double>* __restrict psi0, double* __restrict v0, std::complex<double>* targ, int nelec) {
	if (first_stepOS_U2TU) {
		DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nPts);
		DftiSetValue(dftiHandle, DFTI_NUMBER_OF_TRANSFORMS, nelec);
		DftiSetValue(dftiHandle, DFTI_INPUT_DISTANCE, nPts);
		DftiSetValue(dftiHandle, DFTI_BACKWARD_SCALE, 1.0 / nPts);
		//DftiSetValue(dftiHandle, DFTI_THREAD_LIMIT, numThreads);
		DftiCommitDescriptor(dftiHandle);

		osKineticPhase = new std::complex<double>[nPts];
		osPotentialPhase = new std::complex<double>[nPts];
		//tempD = new double[nPts];
		std::complex<double> dphs = -PhysCon::im * PhysCon::hbar / (2.0 * PhysCon::me) * std::pow(2.0 * PhysCon::pi / ((nPts)*dx), 2) * dt;
		osKineticPhase[0] = 1.0;
		for (int i = 1; i < nPts / 2 + 1; i++) {
			osKineticPhase[i] = std::exp(dphs * (std::complex<double>)(i * i));
			osKineticPhase[nPts - i] = osKineticPhase[i];
		}
		/*if (nPts % 2)
			osKineticPhase[nPts / 2 + 1] = std::exp(dphs * (std::complex<double>)((nPts / 2 + 1) * (nPts / 2 + 1)));*/

		first_stepOS_U2TU = 0;
	}

	std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
	for (int i = 0; i < nPts; i++)
		osPotentialPhase[i] = std::exp(vcnst * v0[i]) * spatialDamp[i];

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

void TDSEIterator1D::stepOS_UW(std::complex<double>* __restrict psi0, double* __restrict v0, std::complex<double>* targ, int nelec) {
	std::complex<double> vcnst = -PhysCon::im * dt / PhysCon::hbar / 2.0;
#pragma omp parallel for
	for (int i = 0; i < nPts; i++)
		osPotentialPhase[i] = std::exp(vcnst * v0[i]) * spatialDamp[i];

#pragma omp parallel for
	for (int i = 0; i < nelec; i++) {
		for (int j = 0; j < nPts; j++) {
			targ[i * nPts + j] = psi0[i * nPts + j] * osPotentialPhase[j];
		}
	}
}


void TDSEIterator1D::multiStepCN(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ, int nelec) {
	// INTEL MKL APPROACH
	//auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
	{
	#pragma omp sections
		{
		#pragma omp section
			{
				std::fill_n(du, nPts - 1, -cnAl);
				vtls::seqMulArrays(nPts - 1, timePhase, du);
			}
		#pragma omp section
			{
				std::fill_n(dl, nPts - 1, -cnAl);
				vtls::seqMulArrays(nPts - 1, &timePhase[1], dl);
			}
		}
	#pragma omp for
		for (int i = 0; i < nPts; i++) {
			d[i] = 1.0 + (2.0 * cnAl * timePhase[i] + cnBe * vf[i]);
		}
	#pragma omp for
		for (int j = 0; j < nelec; j++) {
			std::complex<double> nAl = timePhase[0] * cnAl;
			targ[j * nPts] = (1.0 - 2.0 * nAl - cnBe * v0[0]) * psi0[j * nPts] + nAl * psi0[j * nPts + 1];
			for (int i = 1; i < nPts - 1; i++) {
				nAl = timePhase[i] * cnAl;
				targ[j * nPts + i] = (1.0 - 2.0 * nAl - cnBe * v0[i]) * psi0[j * nPts + i] + nAl * (psi0[j * nPts + i - 1] + psi0[j * nPts + i + 1]);
			}
			nAl = timePhase[nPts - 1] * cnAl;
			targ[j * nPts + nPts - 1] = (1.0 - 2.0 * nAl - cnBe * v0[nPts - 1]) * psi0[j * nPts + nPts - 1] + nAl * psi0[j * nPts + nPts - 2];
		}

		LAPACKE_zgttrf(nPts, dl, d, du, du2, ipiv);
	#pragma omp for
		for (int i = 0; i < nelec; i++)
			LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl, d, du, du2, ipiv, &targ[i * nPts], nPts);

	}
		/*int sid;
		for (int j = 0; j < nelec; j++) {
			sid = j * nPts;
			nAl = timePhase[0] * cnAl;
			targ[sid] = (1.0 - 2.0*nAl - cnBe * v0[0])*psi0[sid] + nAl * psi0[sid+1];
			for (int i = 1; i < nPts - 1; i++) {
				nAl = timePhase[i] * cnAl;
				targ[sid+i] = (1.0 - 2.0*nAl - cnBe * v0[i])*psi0[sid+i] + nAl * (psi0[sid+i - 1] + psi0[sid+i + 1]);
			}
			nAl = timePhase[nPts - 1] * cnAl;
			targ[sid+nPts - 1] = (1.0 - 2.0*nAl - cnBe * v0[nPts - 1])*psi0[sid+nPts - 1] + nAl * psi0[sid+nPts - 2];
		}

		for (int j = 0; j < nelec; j++) {
			res.emplace_back(
				pool->enqueue([this, j, targ, v0, psi0]
			{
				int sid = j * nPts;
				std::complex<double> nAl = timePhase[0] * cnAl;
				targ[sid] = (1.0 - 2.0 * nAl - cnBe * v0[0]) * psi0[sid] + nAl * psi0[sid + 1];
				for (int i = 1; i < nPts - 1; i++) {
					nAl = timePhase[i] * cnAl;
					targ[sid + i] = (1.0 - 2.0 * nAl - cnBe * v0[i]) * psi0[sid + i] + nAl * (psi0[sid + i - 1] + psi0[sid + i + 1]);
				}
				nAl = timePhase[nPts - 1] * cnAl;
				targ[sid + nPts - 1] = (1.0 - 2.0 * nAl - cnBe * v0[nPts - 1]) * psi0[sid + nPts - 1] + nAl * psi0[sid + nPts - 2];
				return 0;
			}
				)
			);
		}


		for (auto&& result : res)
			result.get();
		res.clear();*/

		/*auto stop = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

		std::cout << duration.count() << "/";*/

		//LAPACKE_zgtsv(LAPACK_COL_MAJOR, nPts, nelec, dl, d, du, targ, nPts);

		//start = std::chrono::high_resolution_clock::now();

		//LAPACKE_zgttrf(nPts, dl, d, du, du2, ipiv);

		/*stop = std::chrono::high_resolution_clock::now();

		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

		std::cout << duration.count() << "/";

		start = std::chrono::high_resolution_clock::now();*/


	/*int numPer = nelec / (numThreads - 1);

	for (int i = 0; i < numThreads-1; i++) {
		res.emplace_back(
			pool->enqueue([i, this, numPer, targ] {
				int ret = LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, numPer, dl, d, du, du2, ipiv, &targ[i * nPts * numPer], nPts);
				return ret;
			})
		);
	}
	res.emplace_back(
		pool->enqueue([this, numPer, targ, nelec] {
			int ret = LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, nelec%(numThreads-1), dl, d, du, du2, ipiv, &targ[(numThreads-1) * nPts * numPer], nPts);
			return ret;
		})
	);

	for (auto&& result : res)
		result.get();
	res.clear();*/

	/*stop = std::chrono::high_resolution_clock::now();

	duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

	std::cout << duration.count() << std::endl;*/

	/*for (int i = 0; i < nelec / 16; i++) {
		for (int j = 0; j < 16; j++)
			fs[j] = std::async(r1, this, &targ[(i * 16 + j)*nPts]);
		for (int j = 0; j < 16; j++)
			fs[j].get();
	}
	for (int j = 0; j < nelec%16; j++)
		fs[j] = std::async(r1, this, &targ[(nelec-nelec%16 + j)*nPts]);
	for (int j = 0; j < nelec % 16; j++)
		fs[j].get();*/

	//start = std::chrono::high_resolution_clock::now();

	//LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, nelec, dl, d, du, du2, ipiv, targ, nPts);

	//stop = std::chrono::high_resolution_clock::now();

	//duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

	//std::cout << duration.count() << std::endl;
	//~18.5s for 1,000,000 points 100 times
	//~18s for 100,000 points 1000 times
}

void TDSEIterator1D::stepCNREF(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ) {
	// INTEL MKL APPROACH
	std::complex<double> tp;

	std::fill_n(du, nPts - 1, -cnAl);
	vtls::seqMulArrays(nPts - 1, timePhase, du);
	std::fill_n(dl, nPts - 1, -cnAl);
	for (int i = 0; i < nPts - 1; i++)
		dl[i] *= timePhase[i + 1];
	for (int i = 0; i < nPts; i++) {
		d[i] = 1.0 + (2.0*cnAl * timePhase[i] + cnBe * vf[i]);
	}
	tp = timePhase[0];
	std::complex<double> nAl = tp * cnAl;

	vtls::copyArray(nPts, d, d0);
	vtls::copyArray(nPts - 1, dl, dl0);
	vtls::copyArray(nPts - 1, du, du0);

	temp[0] = (1.0 - 2.0*nAl - cnBe * v0[0])*psi0[0] + nAl * psi0[1];
	for (int i = 1; i < nPts - 1; i++) {
		tp = timePhase[i];
		nAl = tp * cnAl;
		temp[i] = (1.0 - 2.0*nAl - cnBe * v0[i])*psi0[i] + nAl * (psi0[i - 1] + psi0[i + 1]);
	}
	tp = timePhase[nPts - 1];
	nAl = tp * cnAl;
	temp[nPts - 1] = (1.0 - 2.0*nAl - cnBe * v0[nPts - 1])*psi0[nPts - 1] + nAl * psi0[nPts - 2];
	LAPACKE_zgttrf(nPts, dl, d, du, du2, ipiv);
	LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl, d, du, du2, ipiv, temp, nPts);

	double ferr, berr;

	LAPACKE_zgtrfs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl0, d0, du0, dl, d, du, du2, ipiv, temp, nPts, targ, nPts, &ferr, &berr);
	//~18.5s for 1,000,000 points 100 times
	//~18s for 100,000 points 1000 times
}

void TDSEIterator1D::stepCN(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ, double ndt) {
	// INTEL MKL APPROACH
	std::complex<double> tp;

	std::fill_n(du, nPts - 1, -cnAl*ndt/dt);
	vtls::seqMulArrays(nPts - 1, timePhase, du);
	std::fill_n(dl, nPts - 1, -cnAl * ndt / dt);
	for (int i = 0; i < nPts - 1; i++)
		dl[i] *= timePhase[i + 1];
	for (int i = 0; i < nPts; i++) {
		d[i] = 1.0 + (2.0*cnAl*ndt / dt * timePhase[i] + cnBe * ndt / dt * vf[i]);
	}
	tp = timePhase[0];
	std::complex<double> nAl = tp * cnAl*ndt / dt;
	targ[0] = (1.0 - 2.0*nAl - cnBe * ndt / dt * v0[0])*psi0[0] + nAl * psi0[1];
	for (int i = 1; i < nPts - 1; i++) {
		tp = timePhase[i];
		nAl = tp * cnAl*ndt / dt;
		targ[i] = (1.0 - 2.0*nAl - cnBe * ndt / dt * v0[i])*psi0[i] + nAl * (psi0[i - 1] + psi0[i + 1]);
	}
	tp = timePhase[nPts - 1];
	nAl = tp * cnAl*ndt / dt;
	targ[nPts - 1] = (1.0 - 2.0*nAl - cnBe * ndt / dt * v0[nPts - 1])*psi0[nPts - 1] + nAl * psi0[nPts - 2];
	LAPACKE_zgttrf(nPts, dl, d, du, du2, ipiv);
	LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl, d, du, du2, ipiv, targ, nPts);
	//~18.5s for 1,000,000 points 100 times
	//~18s for 100,000 points 1000 times
}

void TDSEIterator1D::stepCNREF(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ, double ndt) {
	// INTEL MKL APPROACH
	std::complex<double> tp;

	std::fill_n(du, nPts - 1, -cnAl * ndt / dt);
	vtls::seqMulArrays(nPts - 1, timePhase, du);
	std::fill_n(dl, nPts - 1, -cnAl * ndt / dt);
	for (int i = 0; i < nPts - 1; i++)
		dl[i] *= timePhase[i + 1];
	for (int i = 0; i < nPts; i++) {
		d[i] = 1.0 + (2.0*cnAl*ndt / dt * timePhase[i] + cnBe * ndt / dt * vf[i]);
	}
	tp = timePhase[0];
	std::complex<double> nAl = tp * cnAl*ndt / dt;

	vtls::copyArray(nPts, d, d0);
	vtls::copyArray(nPts - 1, dl, dl0);
	vtls::copyArray(nPts - 1, du, du0);

	temp[0] = (1.0 - 2.0*nAl - cnBe * ndt / dt * v0[0])*psi0[0] + nAl * psi0[1];
	for (int i = 1; i < nPts - 1; i++) {
		tp = timePhase[i];
		nAl = tp * cnAl*ndt / dt;
		temp[i] = (1.0 - 2.0*nAl - cnBe * ndt / dt * v0[i])*psi0[i] + nAl * (psi0[i - 1] + psi0[i + 1]);
	}
	tp = timePhase[nPts - 1];
	nAl = tp * cnAl*ndt / dt;
	temp[nPts - 1] = (1.0 - 2.0*nAl - cnBe * ndt / dt * v0[nPts - 1])*psi0[nPts - 1] + nAl * psi0[nPts - 2];
	LAPACKE_zgttrf(nPts, dl, d, du, du2, ipiv);
	LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl, d, du, du2, ipiv, temp, nPts);

	double ferr, berr;

	LAPACKE_zgtrfs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl0, d0, du0, dl, d, du, du2, ipiv, temp, nPts, targ, nPts, &ferr, &berr);
	//~18.5s for 1,000,000 points 100 times
	//~18s for 100,000 points 1000 times
}

void TDSEIterator1D::stepH2Euler(std::complex<double> *__restrict psi0, double *__restrict v0, std::complex<double> * targ){

	targ[0] = psi0[0] + (-(6.0*h2al + 4 * h2be * v0[0] + std::pow(v0[0], 2) - h2be * (v0[1] - 2.0*v0[0]))*psi0[0] +
		(4.0*h2al + 2.0*h2be*v0[0])*(psi0[1]) -
		h2al * (psi0[2]))*h2dt;

	targ[1] = psi0[1] + (-(6.0*h2al + 4 * h2be * v0[1] + std::pow(v0[1], 2) - h2be * (v0[2] + v0[0] - 2.0*v0[1]))*psi0[1] +
		(4.0*h2al + 2.0*h2be*v0[1])*(psi0[0] + psi0[2]) -
		h2al * (psi0[3]))*h2dt;

	for (int i = 2; i < nPts - 2; i++) {
		targ[i] = psi0[i] +(-(6.0*h2al + 4 * h2be * v0[i] + std::pow(v0[i], 2) - h2be * (v0[i + 1] + v0[i - 1] - 2.0*v0[i]))*psi0[i] +
			(4.0*h2al + 2.0*h2be*v0[i])*(psi0[i - 1] + psi0[i + 1]) -
			h2al * (psi0[i - 2] + psi0[i + 2]))*h2dt;
	}

	targ[nPts-2] = psi0[nPts-2] + (-(6.0*h2al + 4 * h2be * v0[nPts-2] + std::pow(v0[nPts-2], 2) - h2be * (v0[nPts-1] + v0[nPts-3] - 2.0*v0[nPts-2]))*psi0[nPts-2] +
		(4.0*h2al + 2.0*h2be*v0[nPts-2])*(psi0[nPts-3] + psi0[nPts-1]) -
		h2al * (psi0[nPts-4]))*h2dt;

	targ[nPts-1] = psi0[nPts-1] + (-(6.0*h2al + 4 * h2be * v0[nPts-1] + std::pow(v0[nPts-1], 2) - h2be * (v0[nPts-2] - 2.0*v0[nPts-1]))*psi0[nPts-1] +
		(4.0*h2al + 2.0*h2be*v0[nPts-1])*(psi0[nPts-2]) -
		h2al * (psi0[nPts-3]))*h2dt;
}


void TDSEIterator1D::stepCNImT(std::complex<double> *__restrict psi0, double *__restrict v, std::complex<double> * targ, double rate) {
	//IMPLEMENT CRANK_NICHOLSON IMAGINARY TIME, ASSUME CONSTANT POTENTIAL
	double nAl = std::imag(cnAl)*rate;
	double nBe = std::imag(cnBe)*rate;
	std::fill_n(du, nPts - 1, -nAl);
	std::fill_n(dl, nPts - 1, -nAl);
	for (int i = 0; i < nPts; i++) {
		d[i] = 1.0 + (2.0*nAl + nBe * v[i]);
	}
	targ[0] = (1.0 - 2.0*nAl - nBe * v[0])*psi0[0] + nAl * psi0[1];
	for (int i = 1; i < nPts - 1; i++) {
		targ[i] = (1.0 - 2.0*nAl - nBe * v[i])*psi0[i] + nAl * (psi0[i - 1] + psi0[i + 1]);
	}
	targ[nPts - 1] = (1.0 - 2.0*nAl - nBe * v[nPts - 1])*psi0[nPts - 1] + nAl * psi0[nPts - 2];

	LAPACKE_zgttrf(nPts, dl, d, du, du2, ipiv);
	LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', nPts, 1, dl, d, du, du2, ipiv, targ, nPts);
}

void TDSEIterator1D::stepCNPA_OBS(std::complex<double> *__restrict psi0, double *__restrict v0, double *__restrict vf, std::complex<double> * targ) {

	a[0] = 1.0 + (2.0*cnAl * timePhase[0] + cnBe * vf[0]);
	for (int i = 1; i < nPts - 1; i++)
		a[3 * i] = 1.0 + (2.0*cnAl * timePhase[i] + cnBe * vf[i]);
	a[3 * nPts - 3] = 1.0 + (2.0*cnAl * timePhase[nPts - 1] + cnBe * vf[nPts - 1]);

	std::complex<double> nAl;
	for (int i = 1; i < nPts - 1; i++) {
		nAl = timePhase[i] * cnAl;
		targ[i] = (1.0 - 2.0*nAl - cnBe * v0[i])*psi0[i] + nAl * (psi0[i - 1] + psi0[i + 1]);
	}

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nPts, a, ia, ja, perm, &nrhs, iparm, &msglvl, targ, x, &error);
}

void TDSEIterator1D::stepCNSL(std::complex<double> * psi0, double * v0, double * vf, std::complex<double> * targ) {
	//IMPLEMENT CRANK-NICHOLSON
	// STANDARD LIBRARY APPROACH
	//create d vector
	cnd[0] = (1.0 - 2.0*cnAl - cnBe * v0[0])*psi0[0] + cnAl * psi0[1];
	for (int i = 0; i < nPts - 1; i++) {
		cnd[i] = (1.0 - 2.0*cnAl - cnBe * v0[i])*psi0[i] + cnAl * (psi0[i - 1] + psi0[i + 1]);
	}
	cnd[nPts - 1] = (1.0 - 2.0*cnAl - cnBe * v0[nPts - 1])*psi0[nPts - 1] + cnAl * psi0[nPts - 2];

	//create b vector
	for (int i = 0; i < nPts; i++) cnb[i] = 1.0 + 2.0 * cnAl*vf[i];

	//create c' vector
	cncp[0] = -cnAl / cnb[0];
	for (int i = 0; i < nPts - 1; i++) {
		cncp[i] = -cnAl / (cnb[i] + cnAl * cncp[i]);
	}

	//create d' vector
	cndp[0] = cnd[0] / cnb[0];
	for (int i = 0; i < nPts - 1; i++) {
		cndp[i] = (cnd[i] + cnAl * cndp[i - 1]) / (cnb[i] + cnAl * cncp[i - 1]);
	}

	//get new wave function
	targ[nPts - 1] = cndp[nPts - 1];
	for (int i = nPts - 1; i >= 0; i--) {
		targ[i] = cndp[i] - cncp[i] * targ[i + 1];
	}
	//~15.25s for 1,000,000 points 100 times
	//~14.4s for 100,000 points 1000 times
	//but possible loss of precision due to
	//recursive nature of Thomas algorithm
}