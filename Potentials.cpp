#include "stdafx.h"

namespace Potentials {
	//Calculates a prefactor such that prefactor[i]*psi2[i*nelec] provides the 3-D electron density
	//	as a function of x
	/*void calcFermiBoxDimensionalityConversion(int nelec, int nPts, double dx, double ef, std::complex<double> * psi, Potential * totPot, double* prefactor) {
		if (nelec < 1)
			throw("calcFermiBoxDimensionalityConversion: Number of electrons is not finite! Failed to initialize.");
		//prefactor = new double[nelec];
		std::complex<double>* scratch1 = new std::complex<double>[nPts];
		std::complex<double>* scratch2 = new std::complex<double>[nPts];
		double* v = new double[nPts];
		totPot->getV(0, v);
		double minE, maxE;
		for (int i = 0; i < nelec; i++) {
			vtls::secondDerivative(nPts, &psi[i * nPts], scratch1, dx);
			vtls::scaMulArray(nPts, -PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me), scratch1);
			vtls::seqMulArrays(nPts, v, &psi[i * nPts], scratch2);
			vtls::addArrays(nPts, scratch2, scratch1);
			for (int j = 0; j < nPts; j++)
				scratch2[j] = std::conj(psi[i * nPts + j]);
			prefactor[i] = std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));
			if (i == 0) { minE = prefactor[i]; maxE = prefactor[i]; }
			if (prefactor[i] < minE) minE = prefactor[i];
			if (prefactor[i] > maxE) maxE = prefactor[i];
		}
		//Find center of Fermi slab, and set new 0-energy accoridngly, and convert to Fermi energy difference
		double bottom = (maxE + minE) / 2 - ef / 2;
		for (int i = 0; i < nelec; i++)
			prefactor[i] = ef - (prefactor[i] - bottom);
		//Get factored prefactor
		double fact = 0;
		for (int i = 0; i < nelec; i++)
			fact += prefactor[i];
		fact = 2 / (3 * PhysCon::pi) * PhysCon::me * ef / (PhysCon::hbar * PhysCon::hbar) * nelec / fact;
		//Combine factored prefactor with energy difference
		for (int i = 0; i < nelec; i++)
			prefactor[i] *= fact;

		delete[] scratch1, scratch2, v;

		//psi2 = new double[nelec * nPts];
		//if (nelec < 1)
		//	throw("SurfaceSpaceCharge: Number of electrons is not finite! Failed to initialize.");
		//prefactor = new double[nelec];
		//std::complex<double>* scratch1 = new std::complex<double>[nPts];
		//std::complex<double>* scratch2 = new std::complex<double>[nPts];
		//double* v = new double[nPts];
		//totPot->getV(0, v);
		//double minE, maxE;
		//for (int i = 0; i < nelec; i++) {
		//	vtls::secondDerivative(nPts, &psi[i*nPts], scratch1, dx);
		//	vtls::scaMulArray(nPts, -PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me), scratch1);
		//	vtls::seqMulArrays(nPts, v, &psi[i*nPts], scratch2);
		//	vtls::addArrays(nPts, scratch2, scratch1);
		//	for (int j = 0; j < nPts; j++)
		//		scratch2[j] = std::conj(psi[i*nPts+j]);
		//	prefactor[i] = std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));
		//	if (i == 0) { minE = prefactor[i]; maxE = prefactor[i]; }
		//	if (prefactor[i] < minE) minE = prefactor[i];
		//	if (prefactor[i] > maxE) maxE = prefactor[i];
		//}
		////Find center of Fermi slab, and set new 0-energy accoridngly, and convert to Fermi energy difference
		//double bottom = (maxE + minE) / 2 - ef / 2;
		//for (int i = 0; i < nelec; i++)
		//	prefactor[i] = ef - (prefactor[i] - bottom);
		////Get factored prefactor
		//double fact = 0;
		//for (int i = 0; i < nelec; i++)
		//	fact += prefactor[i];
		//fact = 2 / (3 * PhysCon::pi) * PhysCon::me * ef / (PhysCon::hbar * PhysCon::hbar) * nelec / fact * PhysCon::qe * PhysCon::qe / PhysCon::e0 / 2;
		////Combine factored prefactor with energy difference
		//for (int i = 0; i < nelec; i++)
		//	prefactor[i] *= fact;
	}*/

	FilePotential::FilePotential(int nPts, double * x, double offset, const char * fil, int refPoint) {
		FilePotential::nPts = nPts;
		v = new double[nPts];
		std::fstream ifil = std::fstream(fil, std::ios::in | std::ios::binary);
		int nRep;
		ifil.read(reinterpret_cast<char*>(&nRep), sizeof(int));

		double * fx = new double[nRep];
		double * fv = new double[nRep];
		ifil.read(reinterpret_cast<char*>(fx), sizeof(double)*nRep);
		ifil.read(reinterpret_cast<char*>(fv), sizeof(double)*nRep);
		for (int i = 0; i < nRep; i++)
			fx[i] += offset;
		vtls::linearInterpolate(nRep, fx, fv, nPts, x, v);
		vtls::scaAddArray(nPts, -v[refPoint], v);
		ifil.close();
		delete[] fx;
		delete[] fv;
	}

	void FilePotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void FilePotential::getV(std::complex<double> * psi, double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	int FilePotential::isDynamic() {
		return 0;
	}

	BiasFieldPotential::BiasFieldPotential(int nPts, double * x, double tstart, double tbuf, double xmin, double xmax, double xmin_buf, double xmax_buf, double fieldStrength, int refPoint) {
		BiasFieldPotential::nPts = nPts;
		BiasFieldPotential::tstart = tstart;
		BiasFieldPotential::tbuf = tbuf;
		double dx = (x[nPts - 1] - x[0]) / nPts;

		v = new double[nPts];
		double cx;
		for (int i = 0; i < nPts; i++) {
			cx = x[i];
			if (cx < xmin)
				v[i] = 0;
			else if (cx > xmin && cx < xmax) {
				if (cx > xmin + xmin_buf && cx < xmax - xmax_buf)
					v[i] = v[i-1] + dx*fieldStrength*PhysCon::qe;
				else if (cx < xmin + xmin_buf) {
					double k = (cx - xmin) / xmin_buf;
					v[i] = v[i - 1] + dx * fieldStrength*PhysCon::qe * (
						924.0*std::pow(k, 13) -
						6006.0*std::pow(k, 12) +
						16380.0*std::pow(k, 11) -
						24024.0*std::pow(k, 10) +
						20020.0*std::pow(k, 9) -
						9009.0*std::pow(k, 8) +
						1716.0*std::pow(k, 7)
						);
				}
				else if(cx < xmax) {
					double k = (xmax-cx) / xmax_buf;
					v[i] = v[i - 1] + dx * fieldStrength*PhysCon::qe * (
						924.0*std::pow(k, 13) -
						6006.0*std::pow(k, 12) +
						16380.0*std::pow(k, 11) -
						24024.0*std::pow(k, 10) +
						20020.0*std::pow(k, 9) -
						9009.0*std::pow(k, 8) +
						1716.0*std::pow(k, 7)
						);
				}
			}
			else
				v[i] = v[i-1];
		}
		vtls::scaAddArray(nPts, -v[refPoint], v);
	}

	void BiasFieldPotential::getV(double t, double * targ) {
		if (t > tstart + tbuf)
			vtls::copyArray(nPts, v, targ);
		else if (t > tstart) {
			double k = (t-tstart) / tbuf;
			vtls::scaMulArray(nPts,
				(
					924.0*std::pow(k, 13) -
					6006.0*std::pow(k, 12) +
					16380.0*std::pow(k, 11) -
					24024.0*std::pow(k, 10) +
					20020.0*std::pow(k, 9) -
					9009.0*std::pow(k, 8) +
					1716.0*std::pow(k, 7)
					),
				v, targ);
		}
		else
			std::fill_n(targ, nPts, 0.0);
	}

	void BiasFieldPotential::getV(std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	int BiasFieldPotential::isDynamic() {
		return 1;
	}


	CoulombPotential::CoulombPotential(int nPts, double * x, double ne, double chargePos, double minX, double maxX, int refPoint) {
		v = new double[nPts];
		double k = PhysCon::qe*PhysCon::qe / (4.0*PhysCon::pi*PhysCon::e0);
		double dx = (x[nPts - 1] - x[0]) / nPts;
		for (int i = 0; i < nPts; i++) {
			if (x[i] < minX)
				v[i] = -k / std::abs(minX - chargePos);
			else if (x[i] > maxX)
				v[i] = -k / std::abs(maxX - chargePos);
			else
				v[i] = -k / std::max(std::abs(x[i] - chargePos), dx);
		}
		double ref = v[refPoint];
		for (int i = 0; i < nPts; i++)
			v[i] -= ref;
	}

	void CoulombPotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void CoulombPotential::getV(std::complex<double> * psi, double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	int CoulombPotential::isDynamic() {
		return 0;
	}

	FiniteBox::FiniteBox(int nPts, double * x, double left, double right, double vin, int refPoint) {
		FiniteBox::nPts = nPts;
		v = new double[nPts];
		for (int i = 0; i < nPts; i++) {
			if (x[i] > left && x[i] < right)
				v[i] = vin;
			else
				v[i] = 0;
		}
		vtls::scaAddArray(nPts, -v[refPoint], v);
	}

	void FiniteBox::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void FiniteBox::getV(std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	int FiniteBox::isDynamic() {
		return 0;
	}



	JelliumPotential::JelliumPotential(int nPts, double * x, double center, double ef, double w, int refPoint) {
		JelliumPotential::nPts = nPts;
		v = new double[nPts];
		double nEf = ef / PhysCon::auE_ry;
		double nW = w / PhysCon::auE_ry;
		double v0 = nEf + nW;
		double nBulk = 1.0 / (3.0*std::pow(PhysCon::pi, 2))*std::pow(nEf, 2);
		double rs = std::cbrt(3.0 / (4.0*PhysCon::pi*nBulk));
		double zim = -0.2*rs + 1.25;
		double kf = std::sqrt(nEf);
		double b = kf;
		double aA = -1.0+2.0*v0/b;//4.0 * v0 / b - 1.0;
		double bB = v0/aA;//v0 / (4.0 * v0 / b - 1.0);
		for (int i = 0; i < nPts; i++) {
			double xc = (x[i] - center) / PhysCon::a0;
			if (xc < zim) {
				v[i] = -v0 / (aA*std::exp(bB*(xc-zim)) + 1.0) * PhysCon::auE_ry;
			}
			else {
				if (xc == 0)
					if (i >= 2)
						v[i] = 2.0 * v[i - 1] - v[i - 2];
					else if (i == 1)
						v[i] = v[i - 1];
					else
						v[i] = 0;
				else
					v[i] = -1.0 / (2.0 * (xc - zim)) * (1.0 - std::exp(-b * (xc - zim)))*PhysCon::auE_ry;
				/*-(1.0 - std::exp(-b * xc)) /
					(4.0*xc)*
					(-v0 / (aA*std::exp(bB*zim) + 1.0)) /
					(-(1.0 - std::exp(-b * zim)) /
					(4.0*zim))*
					PhysCon::auE;*/
			}
		}
		vtls::scaAddArray(nPts, -v[refPoint], v);
	}

	void JelliumPotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void JelliumPotential::getV(std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	int JelliumPotential::isDynamic() {
		return 0;
	}


	JelliumPotentialBacked::JelliumPotentialBacked(int nPts, double * x, double center, double ef, double w, double backStart, double backWidth, int refPoint) {
		JelliumPotentialBacked::nPts = nPts;
		v = new double[nPts];
		double nEf = ef / PhysCon::auE_ry;
		double nW = w / PhysCon::auE_ry;
		double v0 = nEf + nW;
		double nBulk = 1.0 / (3.0 * std::pow(PhysCon::pi, 2)) * std::pow(nEf, 2);
		double rs = std::cbrt(3.0 / (4.0 * PhysCon::pi * nBulk));
		double zim = -0.2 * rs + 1.25;
		double kf = std::sqrt(nEf);
		double b = kf;
		double aA = -1.0 + 2.0 * v0 / b;//4.0 * v0 / b - 1.0;
		double bB = v0 / aA;//v0 / (4.0 * v0 / b - 1.0);
		for (int i = 0; i < nPts; i++) {
			double xc = (x[i] - center) / PhysCon::a0;
			if (xc < zim) {
				v[i] = -v0 / (aA * std::exp(bB * (xc - zim)) + 1.0) * PhysCon::auE_ry;
			}
			else {
				if (xc == 0)
					if (i >= 2)
						v[i] = 2.0 * v[i - 1] - v[i - 2];
					else if (i == 1)
						v[i] = v[i - 1];
					else
						v[i] = 0;
				else
					v[i] = -1.0 / (2.0 * (xc - zim)) * (1.0 - std::exp(-b * (xc - zim))) * PhysCon::auE_ry;
				/*-(1.0 - std::exp(-b * xc)) /
					(4.0*xc)*
					(-v0 / (aA*std::exp(bB*zim) + 1.0)) /
					(-(1.0 - std::exp(-b * zim)) /
					(4.0*zim))*
					PhysCon::auE;*/
			}
		}
		vtls::scaAddArray(nPts, -v[refPoint], v);
		for (int i = 0; x[i] < backStart + backWidth / 2 && i < nPts; i++) {
			double k;
			if (x[i] > backStart - backWidth / 2) {
				k = (x[i] - backStart + backWidth / 2) / backWidth;
				v[i] *= (
					924.0*std::pow(k, 13) -
					6006.0*std::pow(k, 12) +
					16380.0*std::pow(k, 11) -
					24024.0*std::pow(k, 10) +
					20020.0*std::pow(k, 9) -
					9009.0*std::pow(k, 8) +
					1716.0*std::pow(k, 7));
			}
			else
				v[i] = 0.0;
		}

	}

	void JelliumPotentialBacked::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void JelliumPotentialBacked::getV(std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	int JelliumPotentialBacked::isDynamic() {
		return 0;
	}


	ElectricFieldProfileToPotential::ElectricFieldProfileToPotential(int nPts, ElectricFieldProfiles::ElectricFieldProfile * fieldProfile, double dx, double phase, double tmax, double lam, Envelopes::Envelope * env, int refPoint) {
		ElectricFieldProfileToPotential::tmax = tmax;
		ElectricFieldProfileToPotential::env = env;
		ElectricFieldProfileToPotential::nPts = nPts;
		ElectricFieldProfileToPotential::phase = phase;
		std::complex<double> * fieldMask = fieldProfile->getProfile();
		potMask = new std::complex<double>[nPts];
		vtlsInt::cumIntTrapz(nPts, fieldMask, -dx * PhysCon::qe, potMask);
		vtls::scaAddArray(nPts, -potMask[refPoint], potMask);
		w = PhysCon::c / lam * 2.0*PhysCon::pi;
	}

	void ElectricFieldProfileToPotential::getV(double t, double * targ) {
		vtls::scaMulArrayRe(nPts, std::exp(PhysCon::im*(w*(t-tmax)+phase))*PhysCon::im*env->getValue(t), potMask, targ);
	}

	void ElectricFieldProfileToPotential::getV(std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	int ElectricFieldProfileToPotential::isDynamic() {
		return 1;
	}


	ShieldedAtomicPotential::ShieldedAtomicPotential(int nPts, double * x, double center, double latticeSpacing, double zProtons, double decayConst) {
		ShieldedAtomicPotential::nPts = nPts;
		v = new double[nPts];
		using namespace PhysCon;
		for (int i = 0; i < nPts; i++)
			v[i] = -zProtons * qe*qe / (2 * e0*latticeSpacing*latticeSpacing / decayConst)*std::exp(-std::abs(x[i] - center) / decayConst);
	}

	void ShieldedAtomicPotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void ShieldedAtomicPotential::getV(std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	int ShieldedAtomicPotential::isDynamic() {
		return 0;
	}


	WaveFunctionSelfPotentialJellPotMask::WaveFunctionSelfPotentialJellPotMask(int nPts, double * x, double dx, double strength, double otherDimensionDistance, double center, double ef, double w, int refPoint) {
		WaveFunctionSelfPotentialJellPotMask::nPts = nPts;
		WaveFunctionSelfPotentialJellPotMask::dx = dx;
		WaveFunctionSelfPotentialJellPotMask::strength = strength;
		WaveFunctionSelfPotentialJellPotMask::odimDist = otherDimensionDistance;
		WaveFunctionSelfPotentialJellPotMask::refPoint = refPoint;
		psi2 = new double[nPts];
		mask = new double[nPts];
		stepFunc = new double[nPts];
		add = new double[nPts];
		for (int i = 0; i < nPts; i++)
			add[i] = 0.0;
		int mid = nPts / 2;
		for (int i = 0; i < nPts; i++)
			mask[i] = PhysCon::qe*PhysCon::qe*strength*dx / std::sqrt(std::pow((i - mid)*dx, 2) + odimDist * odimDist);
		double b = std::sqrt(ef * 2.0 / PhysCon::auE_ha);
		double aA = 4.0*(ef + w) / PhysCon::auE_ha / b - 1.0;
		double bB = (ef + w) / PhysCon::auE_ha / aA;
		for (int i = 0; i < nPts; i++)
			stepFunc[i] = 1.0 / (aA*std::exp(-bB * (x[i] - center) / PhysCon::a0) + 1.0);

		initializeConvolution();
	}

	void WaveFunctionSelfPotentialJellPotMask::initializeConvolution() {
		status = VSL_STATUS_OK;
		task_ptr = &task;

		vsldConvNewTaskX1D(task_ptr, VSL_CONV_MODE_FFT, nPts, nPts, nPts, mask, 1);
		int nPts2 = nPts / 2;
		vslConvSetStart(task, &nPts2);
		int one = 1;
		vslConvSetDecimation(task, &one);
	}

	void WaveFunctionSelfPotentialJellPotMask::negateGroundEffects(std::complex<double> * psi) {
		double * temp = new double[nPts];
		getV(psi, 0.0, temp);
		for (int i = 0; i < nPts; i++)
			add[i] = -temp[i] + add[i];
		delete temp;
	}

	void WaveFunctionSelfPotentialJellPotMask::getV(double t, double * targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void WaveFunctionSelfPotentialJellPotMask::getV(std::complex<double> * psi, double t, double * targ) {
		vtls::normSqr(nPts, psi, psi2);
		vtls::seqMulArrays(nPts, stepFunc, psi2);

		status = vsldConvExecX1D(task, psi2, 1, targ, 1);

		if (status != VSL_STATUS_OK) {
			std::cout << "ERROR: convolution: bad status=" << status << std::endl;
			std::cin.ignore();
			exit(1);
		}

		vtls::addArrays(nPts, add, targ);
		vtls::scaAddArray(nPts, -targ[refPoint], targ);

		//vtls::mkl_ddcon(mask, 1, psi2, 1, targ, 1, nPts, nPts, nPts / 2, nPts, 1);
		/*
		for (int i = 0; i < nPts; i++) {
		targ[i] = 0.0;
		for (int j = std::max(i-5, 0); j < std::min(j+5, nPts); j++)
		targ[i] += std::norm(psi[j])*PhysCon::qe*PhysCon::k*strength*dx / std::sqrt(std::pow((i - j)*dx, 2) + odimDist * odimDist);
		}
		for (int i = 0; i < nPts; i++) {
		targ[i] -= targ[refPoint];
		}*/
	}

	int WaveFunctionSelfPotentialJellPotMask::isDynamic() {
		return 2;
	}

	/*
	WaveFunctionSelfPotential::WaveFunctionSelfPotential(int nPts, double dx, double strength, double otherDimensionDistance, int refPoint) {
		WaveFunctionSelfPotential::nPts = nPts;
		WaveFunctionSelfPotential::dx = dx;
		WaveFunctionSelfPotential::strength = strength;
		WaveFunctionSelfPotential::odimDist = otherDimensionDistance;
		WaveFunctionSelfPotential::refPoint = refPoint;
		psi2 = new double[nPts];
		mask = new double[nPts*2];
		add = new double[nPts];
		std::fill_n(add, nPts, 0.0);
		int mid = nPts;
		double oconst = PhysCon::qe*PhysCon::qe*strength*dx / (PhysCon::pi*4.0*PhysCon::e0);
		for (int i = 0; i < nPts*2; i++) {
			if (i != mid)
				mask[i] = oconst / std::sqrt(std::pow((i - mid)*dx, 2) + odimDist * odimDist);
			else
				mask[i] = oconst / std::sqrt(dx*dx + odimDist * odimDist);
		}

		initializeConvolution();
	}

	void WaveFunctionSelfPotential::initializeConvolution() {
		status = VSL_STATUS_OK;
		task_ptr = &task;

		vsldConvNewTaskX1D(task_ptr, VSL_CONV_MODE_FFT, nPts*2, nPts, nPts, mask, 1);
		int nPts2 = nPts;
		vslConvSetStart(task, &nPts2);
		int one = 1;
		vslConvSetDecimation(task, &one);
	}

	void WaveFunctionSelfPotential::negateGroundEffects(std::complex<double> * psi) {
		double * temp = new double[nPts];
		getV(psi, 0.0, temp);
		for (int i = 0; i < nPts; i++)
			add[i] = -temp[i] + add[i];
		delete temp;
	}

	void WaveFunctionSelfPotential::getV(double t, double * targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void WaveFunctionSelfPotential::getV(std::complex<double> * psi, double t, double * targ) {
		vtls::normSqr(nPts, psi, psi2);

		status = vsldConvExecX1D(task, psi2, 1, targ, 1);

		if (status != VSL_STATUS_OK) {
			std::cout << "ERROR: Space charge convolution: bad status=" << status << std::endl;
			std::cin.ignore();
			exit(1);
		}

		vtls::addArrays(nPts, add, targ);
		vtls::scaAddArray(nPts, -targ[refPoint], targ);
	}

	int WaveFunctionSelfPotential::isDynamic() {
		return 2;
	}
	*/


	SurfaceSpaceCharge::SurfaceSpaceCharge(int nPts, double dx, double ef, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int refPoint) {
		SurfaceSpaceCharge::nPts = nPts;
		SurfaceSpaceCharge::dx = dx;
		SurfaceSpaceCharge::ef = ef;
		SurfaceSpaceCharge::refPoint = refPoint;
		SurfaceSpaceCharge::nelecPtr = nelec;
		SurfaceSpaceCharge::totPot = totPot;
		SurfaceSpaceCharge::wght = wght;
		SurfaceSpaceCharge::dens = dens;
		if (posMin < 0)
			SurfaceSpaceCharge::posMin = 0;
		else
			SurfaceSpaceCharge::posMin = posMin;
		if (posMax > nPts - 1)
			SurfaceSpaceCharge::posMax = nPts - 1;
		else
			SurfaceSpaceCharge::posMax = posMax;

		fldTot = new double[nPts];
		origPot = new double[nPts];
		rho = new double[nPts];
		psi0 = new std::complex<double>[nPts];
	}

	void SurfaceSpaceCharge::negateGroundEffects(std::complex<double>* psi) {
		if (first) doFirst(psi);
		calcPot(psi, origPot);
	}

	void SurfaceSpaceCharge::getV(double t, double* targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void SurfaceSpaceCharge::getV(std::complex<double>* psi, double t, double* targ) {
		if (first) {
			doFirst(psi);
			std::cout << "SurfaceSpaceCharge potential should be negated after initial state is found." << std::endl;
		}
		calcPot(psi, targ);

		double ref = targ[refPoint] - origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	int SurfaceSpaceCharge::isDynamic() {
		return 2;
	}

	void SurfaceSpaceCharge::doFirst(std::complex<double>* psi) {
		first = 0;
		nelec = nelecPtr[0];
		psi2 = new double[nPts * nelec];
		prefactor = new double[nelec];
		double* energies = new double[nelec];
		double* v0w = new double[nPts];
		totPot->getV(0, v0w);
		WfcToRho::calcEnergies(nelec, nPts, dx, psi, v0w, energies);
		wght->calcWeights(nelec, energies, prefactor);
		delete[] energies, v0w;
		//calcFermiBoxDimensionalityConversion(nelec, nPts, dx, ef, psi, totPot, prefactor);
		//psi2 = new double[nelec * nPts];
		//if (nelec < 1)
		//	throw("SurfaceSpaceCharge: Number of electrons is not finite! Failed to initialize.");
		//prefactor = new double[nelec];
		//std::complex<double>* scratch1 = new std::complex<double>[nPts];
		//std::complex<double>* scratch2 = new std::complex<double>[nPts];
		//double* v = new double[nPts];
		//totPot->getV(0, v);
		//double minE, maxE;
		//for (int i = 0; i < nelec; i++) {
		//	vtls::secondDerivative(nPts, &psi[i*nPts], scratch1, dx);
		//	vtls::scaMulArray(nPts, -PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me), scratch1);
		//	vtls::seqMulArrays(nPts, v, &psi[i*nPts], scratch2);
		//	vtls::addArrays(nPts, scratch2, scratch1);
		//	for (int j = 0; j < nPts; j++)
		//		scratch2[j] = std::conj(psi[i*nPts+j]);
		//	prefactor[i] = std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));
		//	if (i == 0) { minE = prefactor[i]; maxE = prefactor[i]; }
		//	if (prefactor[i] < minE) minE = prefactor[i];
		//	if (prefactor[i] > maxE) maxE = prefactor[i];
		//}
		////Find center of Fermi slab, and set new 0-energy accoridngly, and convert to Fermi energy difference
		//double bottom = (maxE + minE) / 2 - ef / 2;
		//for (int i = 0; i < nelec; i++)
		//	prefactor[i] = ef - (prefactor[i] - bottom);
		////Get factored prefactor
		//double fact = 0;
		//for (int i = 0; i < nelec; i++)
		//	fact += prefactor[i];
		//fact = 2 / (3 * PhysCon::pi) * PhysCon::me * ef / (PhysCon::hbar * PhysCon::hbar) * nelec / fact * PhysCon::qe * PhysCon::qe / PhysCon::e0 / 2;
		////Combine factored prefactor with energy difference
		//for (int i = 0; i < nelec; i++)
		//	prefactor[i] *= fact;
	}

	void SurfaceSpaceCharge::calcPot(std::complex<double>* psi, double* targ) {
		if (first) doFirst(psi);
		dens->calcRho(nPts, nelec, dx, prefactor, psi, rho);
		for (int i = 0; i < posMin; i++)
			rho[i] = 0;
		for (int i = posMax; i < nPts; i++)
			rho[i] = 0;

		// NEW CALC, STAGGERED FIELD
		vtlsInt::cumIntRectLeft(nPts, rho, dx, fldTot);
		vtlsInt::cumIntRectRight(nPts, fldTot, -dx * std::pow(PhysCon::qe, 2) / PhysCon::e0, targ);

		/* OLD CALC, NON-STAGGERED FIELD
		//Calculate ~field
		vtlsInt::cumIntTrapz(nPts, rho, dx, fldTot);
		//Calculate potential
		vtlsInt::cumIntTrapz(nPts, fldTot, -dx * std::pow(PhysCon::qe, 2) / PhysCon::e0, targ);
		*/

		//if (first) doFirst(psi);
		////Calculate field due to each electron
		//std::fill_n(fldTot, nPts, 0);
		//vtls::normSqr(nPts * nelec, psi, psi2);
		//for (int i = 0; i < nelec; i++) {
		//	vtlsInt::cumIntTrapz(nPts, &psi2[i*nPts], dx * prefactor[i], fldInd);
		//	vtls::addArrays(nPts, fldInd, fldTot);
		//}
		//vtls::scaAddArray(nPts, -fldTot[nPts - 1]/2.0, fldTot);
		////Calculate potential due to field
		//vtlsInt::cumIntTrapz(nPts, fldTot, -dx * std::pow(PhysCon::qe, 2) / PhysCon::e0, targ);
	}


	FullCylindricalSpaceCharge::FullCylindricalSpaceCharge(int nPts, double * x, double dx, double ef, double r, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int surfPos, int refPoint) {
		FullCylindricalSpaceCharge::nPts = nPts;
		FullCylindricalSpaceCharge::dx = dx;
		FullCylindricalSpaceCharge::ef = ef;
		FullCylindricalSpaceCharge::refPoint = refPoint;
		FullCylindricalSpaceCharge::nelecPtr = nelec;
		FullCylindricalSpaceCharge::totPot = totPot;
		FullCylindricalSpaceCharge::r = r;
		FullCylindricalSpaceCharge::wght = wght;
		FullCylindricalSpaceCharge::dens = dens;
		if (posMin < 0)
			FullCylindricalSpaceCharge::posMin = 0;
		else
			FullCylindricalSpaceCharge::posMin = posMin;
		if(posMax > nPts-1)
			FullCylindricalSpaceCharge::posMax = nPts-1;
		else
			FullCylindricalSpaceCharge::posMax = posMax;
		if (surfPos > nPts - 1)
			FullCylindricalSpaceCharge::surfPos = nPts - 1;
		else
			FullCylindricalSpaceCharge::surfPos = surfPos;
		potTemp = new double[nPts];
		origPot = new double[nPts];
		rho = new double[nPts];
		lrxr = new double[nPts];
		for (int i = 0; i < nPts; i++)
			if (x[i] - x[surfPos] <= -r)
				lrxr[i] = 0;
			else
				lrxr[i] = std::log((r + x[i] - x[surfPos]) / r);
	}

	void FullCylindricalSpaceCharge::negateGroundEffects(std::complex<double>* psi) {
		if (first) doFirst(psi);
		calcPot(psi, origPot);
	}

	void FullCylindricalSpaceCharge::getV(double t, double* targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void FullCylindricalSpaceCharge::getV(std::complex<double>* psi, double t, double* targ) {
		if (first) {
			doFirst(psi);
			std::cout << "FullCylindricalSpaceCharge potential should be negated after initial state is found." << std::endl;
		}
		calcPot(psi, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	int FullCylindricalSpaceCharge::isDynamic() {
		return 2;
	}

	void FullCylindricalSpaceCharge::doFirst(std::complex<double>* psi) {
		first = 0;
		nelec = nelecPtr[0];
		psi2 = new double[nPts * nelec];
		prefactor = new double[nelec];
		double* energies = new double[nelec];
		double* v0w = new double[nPts];
		totPot->getV(0, v0w);
		WfcToRho::calcEnergies(nelec, nPts, dx, psi, v0w, energies);
		wght->calcWeights(nelec, energies, prefactor);
		delete[] energies, v0w;
		//calcFermiBoxDimensionalityConversion(nelec, nPts, dx, ef, psi, totPot, prefactor);
	}

	void FullCylindricalSpaceCharge::calcPot(std::complex<double>* psi, double* targ) {
		if (first) doFirst(psi);
		//Caclulate 3-D density, only for within bounds
		dens->calcRho(nPts, nelec, dx, prefactor, psi, rho);
		for (int i = 0; i < posMin; i++)
			rho[i] = 0;
		for (int i = posMax; i < nPts; i++)
			rho[i] = 0;
		//Calculate first integral
		std::fill_n(targ, nPts, 0);
		vtlsInt::cumIntTrapz(nPts, rho, dx, targ);
		vtls::seqMulArrays(nPts, lrxr, targ);
		//Calculate second integral
		std::fill_n(potTemp, nPts, 0);
		vtls::seqMulArrays(nPts, lrxr, rho);
		vtlsInt::cumIntTrapz(nPts, rho, -dx, potTemp);
		vtls::addArrays(nPts, potTemp, targ);
		//Apply final constants
		vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe * r / PhysCon::e0, targ);
	}


	LinearBulkCylindricalFieldSpaceCharge::LinearBulkCylindricalFieldSpaceCharge(int nPts, double* x, double dx, double ef, double rad, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int surfPos, int refPoint) {
		LinearBulkCylindricalFieldSpaceCharge::nPts = nPts;
		LinearBulkCylindricalFieldSpaceCharge::dx = dx;
		LinearBulkCylindricalFieldSpaceCharge::ef = ef;
		LinearBulkCylindricalFieldSpaceCharge::refPoint = refPoint;
		LinearBulkCylindricalFieldSpaceCharge::nelecPtr = nelec;
		LinearBulkCylindricalFieldSpaceCharge::totPot = totPot;
		LinearBulkCylindricalFieldSpaceCharge::rad = rad;
		LinearBulkCylindricalFieldSpaceCharge::wght = wght;
		LinearBulkCylindricalFieldSpaceCharge::dens = dens;
		if (posMin < 0)
			LinearBulkCylindricalFieldSpaceCharge::posMin = 0;
		else
			LinearBulkCylindricalFieldSpaceCharge::posMin = posMin;
		if (posMax > nPts - 1)
			LinearBulkCylindricalFieldSpaceCharge::posMax = nPts - 1;
		else
			LinearBulkCylindricalFieldSpaceCharge::posMax = posMax;
		if (surfPos > nPts - 1)
			LinearBulkCylindricalFieldSpaceCharge::surfPos = nPts - 1;
		else
			LinearBulkCylindricalFieldSpaceCharge::surfPos = surfPos;
		potTemp = new double[nPts];
		genTemp = new double[nPts];
		origPot = new double[nPts];
		rho = new double[nPts];
		lrxr = new double[nPts];

		for (int i = 0; i < nPts; i++)
			if (x[i] - x[surfPos] <= -rad)
				lrxr[i] = 0;
			else
				lrxr[i] = std::log((rad + x[i] - x[surfPos]) / rad);
	}

	void LinearBulkCylindricalFieldSpaceCharge::negateGroundEffects(std::complex<double>* psi) {
		if (first) doFirst(psi);
		calcPot(psi, origPot);
	}

	void LinearBulkCylindricalFieldSpaceCharge::getV(double t, double* targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void LinearBulkCylindricalFieldSpaceCharge::getV(std::complex<double>* psi, double t, double* targ) {
		if (first) {
			doFirst(psi);
			std::cout << "LinearBulkCylindricalFieldSpaceCharge potential should be negated after initial state is found." << std::endl;
		}
		calcPot(psi, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	int LinearBulkCylindricalFieldSpaceCharge::isDynamic() {
		return 2;
	}

	void LinearBulkCylindricalFieldSpaceCharge::doFirst(std::complex<double>* psi) {
		first = 0;
		nelec = nelecPtr[0];
		psi2 = new double[nPts * nelec];
		prefactor = new double[nelec];
		double* energies = new double[nelec];
		double* v0w = new double[nPts];
		totPot->getV(0, v0w);
		WfcToRho::calcEnergies(nelec, nPts, dx, psi, v0w, energies);
		wght->calcWeights(nelec, energies, prefactor);
		delete[] energies, v0w;
		//calcFermiBoxDimensionalityConversion(nelec, nPts, dx, ef, psi, totPot, prefactor);
	}

	void LinearBulkCylindricalFieldSpaceCharge::calcPot(std::complex<double>* psi, double* targ) {
		//auto t1 = std::chrono::system_clock::now();
		if (first) doFirst(psi);
		//Caclulate 3-D density, only for within bounds
		dens->calcRho(nPts, nelec, dx, prefactor, psi, rho);
		for (int i = 0; i < posMin; i++)
			rho[i] = 0;
		for (int i = posMax; i < nPts; i++)
			rho[i] = 0;
		//auto t2 = std::chrono::system_clock::now();
		std::fill_n(targ, nPts, 0);

		//ORIGINAL ATTEMPT
		//CALCULATE FIELDS FROM FIELD ELECTRONS
		//Calculate first integral
		vtlsInt::cumIntTrapz(nPts-surfPos, &rho[surfPos], dx*rad, &targ[surfPos]);
		vtls::seqMulArrays(nPts-surfPos, &lrxr[surfPos], &targ[surfPos]);
		//Calculate second integral
		std::fill_n(potTemp, nPts, 0);
		vtls::seqMulArrays(nPts-surfPos, &lrxr[surfPos], &rho[surfPos], &genTemp[surfPos]);
		vtlsInt::cumIntTrapz(nPts-surfPos, &genTemp[surfPos], -dx*rad, &potTemp[surfPos]);
		vtls::addArrays(nPts-surfPos, &potTemp[surfPos], &targ[surfPos]);

		//CALCULATE FIELDS FROM BULK ELECTRONS
		//Calculate cumulative charge integral within metal

		////surfPos -> surfPos+1
		vtlsInt::cumIntTrapz(surfPos+1, rho, dx, genTemp);
		//Calculate potential outside the metal (as it has an extra constant)

		////surfPos-1 -> surfPos
		vtls::scaMulArray(nPts-surfPos, rad * genTemp[surfPos], &lrxr[surfPos], &potTemp[surfPos]);
		//Calculate potential inside the metal

		////added step to append to potTemp
		vtlsInt::cumIntTrapz(surfPos, genTemp, dx, potTemp);
		potTemp[surfPos] += 0.5 * dx * (genTemp[surfPos - 1] + genTemp[surfPos]);

		//Match potentials at surface, sans the expected difference due to field.
		////genTemp[surfPos-1] -> genTemp[surfPos]
		if(surfPos < nPts && surfPos > 0){
			double diff = potTemp[surfPos] - potTemp[surfPos - 1] - dx * genTemp[surfPos];
			vtls::scaAddArray(nPts - surfPos, -diff, &potTemp[surfPos]);
		}
		//Recombine two potentials
		vtls::addArrays(nPts, potTemp, targ);

		//Apply final constants
		vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe / PhysCon::e0, targ);
		//auto t3 = std::chrono::system_clock::now();

		//std::chrono::duration<double> dur1 = t2 - t1;
		//std::chrono::duration<double> dur2 = t3 - t2;
		//std::cout << dur1.count() << " " << dur2.count() << std::endl;

		//CALCULATING WITHOUT WEIRD POINTER STUFF
		////Fields from field electrons
		////Calculate first integral
		//vtlsInt::cumIntTrapz(nPts, rho, dx * rad, temp1);
		//vtls::seqMulArrays(nPts, lrxr, temp1);
		//vtls::scaAddArray(nPts, -temp1[surfPos], temp1);
		//vtls::addArrays(nPts - surfPos, &temp1[surfPos], &targ[surfPos]);
		////Calculate second integral
		//vtls::seqMulArrays(nPts, lrxr, rho, temp1);
		//vtlsInt::cumIntTrapz(nPts, temp1, -dx * rad, temp2);
		//vtls::scaAddArray(nPts, -temp2[surfPos], temp2);
		//vtls::addArrays(nPts - surfPos, &temp2[surfPos], &targ[surfPos]);

		////Calculate fields from bulk electrons
		////Cumulative charge integral
		//vtlsInt::cumIntTrapz(nPts, rho, dx, temp1);
		////Potential outside metal, simple function
		//vtls::scaMulArray(nPts, rad * temp1[surfPos - 1], lrxr, temp2);
		//vtls::scaAddArray(nPts, -temp2[surfPos], temp2);
		//vtls::addArrays(nPts - surfPos, &temp2[surfPos], &targ[surfPos]);
		////Potential inside metal
		//vtlsInt::cumIntTrapz(nPts, temp1, dx, temp2);
		//vtls::scaAddArray(nPts, -temp2[surfPos], temp2);
		//vtls::addArrays(surfPos, temp2, targ);
	}


	DielectricBulkCylindricalFieldSpaceCharge::DielectricBulkCylindricalFieldSpaceCharge(int nPts, double* x, double dx, double dt, double ef, double rad, double wellWidth, double dampRate, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int surfPos, int refPoint) {
		DielectricBulkCylindricalFieldSpaceCharge::nPts = nPts;
		DielectricBulkCylindricalFieldSpaceCharge::dx = dx;
		DielectricBulkCylindricalFieldSpaceCharge::dt = dt;
		DielectricBulkCylindricalFieldSpaceCharge::ef = ef;
		DielectricBulkCylindricalFieldSpaceCharge::refPoint = refPoint;
		DielectricBulkCylindricalFieldSpaceCharge::nelecPtr = nelec;
		DielectricBulkCylindricalFieldSpaceCharge::totPot = totPot;
		DielectricBulkCylindricalFieldSpaceCharge::rad = rad;
		DielectricBulkCylindricalFieldSpaceCharge::wght = wght;
		DielectricBulkCylindricalFieldSpaceCharge::dens = dens;
		DielectricBulkCylindricalFieldSpaceCharge::wellWidth = wellWidth;
		DielectricBulkCylindricalFieldSpaceCharge::dampRate = dampRate;
		DielectricBulkCylindricalFieldSpaceCharge::xs = x;
		if (posMin < 0)
			DielectricBulkCylindricalFieldSpaceCharge::posMin = 0;
		else
			DielectricBulkCylindricalFieldSpaceCharge::posMin = posMin;
		if (posMax > nPts - 1)
			DielectricBulkCylindricalFieldSpaceCharge::posMax = nPts - 1;
		else
			DielectricBulkCylindricalFieldSpaceCharge::posMax = posMax;
		if (surfPos > nPts - 1)
			DielectricBulkCylindricalFieldSpaceCharge::surfPos = nPts - 1;
		else
			DielectricBulkCylindricalFieldSpaceCharge::surfPos = surfPos;
		potTemp = new double[nPts];
		genTemp = new double[nPts];
		origPot = new double[nPts];
		rho = new double[nPts];
		lrxr = new double[nPts];
		xfsf = new double[nPts];
		xfsf2 = new double[nPts];

		for (int i = 0; i < nPts; i++) {
			xfsf[i] = xs[i] - xs[surfPos];
			xfsf2[i] = xfsf[i] * xfsf[i];
		}

		for (int i = 0; i < nPts; i++)
			if (x[i] - x[surfPos] <= -rad)
				lrxr[i] = 0;
			else
				lrxr[i] = std::log((rad + x[i] - x[surfPos]) / rad);

		//pol = 0; lpol = 0;
	}

	void DielectricBulkCylindricalFieldSpaceCharge::negateGroundEffects(std::complex<double>* psi) {
		if (first) doFirst(psi);
		calcPot(psi, origPot);
	}

	void DielectricBulkCylindricalFieldSpaceCharge::getV(double t, double* targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void DielectricBulkCylindricalFieldSpaceCharge::getV(std::complex<double>* psi, double t, double* targ) {
		if (first) {
			doFirst(psi);
			std::cout << "LinearBulkCylindricalFieldSpaceCharge potential should be negated after initial state is found." << std::endl;
		}
		calcPot(psi, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	int DielectricBulkCylindricalFieldSpaceCharge::isDynamic() {
		return 2;
	}

	void DielectricBulkCylindricalFieldSpaceCharge::doFirst(std::complex<double>* psi) {
		first = 0;
		nelec = nelecPtr[0];
		psi2 = new double[nPts * nelec];
		prefactor = new double[nelec];
		double* energies = new double[nelec];
		double* v0w = new double[nPts];
		totPot->getV(0, v0w);
		WfcToRho::calcEnergies(nelec, nPts, dx, psi, v0w, energies);
		wght->calcWeights(nelec, energies, prefactor);
		delete[] energies, v0w;

		dens->calcRho(nPts, nelec, dx, prefactor, psi, rho);
		//pol = vtlsInt::trapzMul(nPts, xs, rho, dx);

		if (wellWidth <= 0) {
			double x0 = vtlsInt::trapzMul(nPts, xs, rho, dx) / vtlsInt::trapz(nPts, rho, dx);
			double* xmx02 = new double[nPts];
			for (int i = 0; i < nPts; i++)
				xmx02[i] = std::pow(xs[i] - x0, 2);
			wellWidth = 2 * std::sqrt(3) * sqrt(vtlsInt::trapzMul(nPts, xmx02, rho, dx) / vtlsInt::trapz(nPts, rho, dx));
			delete[] xmx02;
		}
	}

	void DielectricBulkCylindricalFieldSpaceCharge::calcPot(std::complex<double>* psi, double* targ) {
		//lpol = pol;
		//auto t1 = std::chrono::system_clock::now();
		if (first) doFirst(psi);
		//Caclulate 3-D density, only for within bounds
		dens->calcRho(nPts, nelec, dx, prefactor, psi, rho);
		for (int i = 0; i < posMin; i++)
			rho[i] = 0;
		for (int i = posMax; i < nPts; i++)
			rho[i] = 0;
		//auto t2 = std::chrono::system_clock::now();
		std::fill_n(targ, nPts, 0);

		//ORIGINAL ATTEMPT
		//CALCULATE FIELDS FROM FIELD ELECTRONS
		//Calculate first integral
		vtlsInt::cumIntTrapz(nPts - surfPos, &rho[surfPos], dx * rad, &targ[surfPos]);
		vtls::seqMulArrays(nPts - surfPos, &lrxr[surfPos], &targ[surfPos]);
		std::cout << targ[surfPos + 1] - targ[surfPos];
		//Calculate second integral
		std::fill_n(potTemp, nPts, 0);
		vtls::seqMulArrays(nPts - surfPos, &lrxr[surfPos], &rho[surfPos], &genTemp[surfPos]);
		vtlsInt::cumIntTrapz(nPts - surfPos, &genTemp[surfPos], -dx * rad, &potTemp[surfPos]);
		vtls::addArrays(nPts - surfPos, &potTemp[surfPos], &targ[surfPos]);
		std::cout << " | " << potTemp[surfPos + 1] - potTemp[surfPos] << std::endl;
		//see if slight amount of field is consistent with charge at surface
		//try printing out the density between surfPos and surfPos+1

		//CALCULATE FIELDS FROM BULK ELECTRONS
		//Calculate charge integral within metal
		double totBulk = vtlsInt::trapz(surfPos, rho, dx);
		//vtlsInt::cumIntTrapz(surfPos, rho, dx, genTemp);
		//Calculate potential outside the metal (as it has an extra constant)
		vtls::scaMulArray(nPts - surfPos, rad * totBulk , &lrxr[surfPos], &potTemp[surfPos]);
		//Calculate potential inside the metal
		//Calculate polarization
		//pol = vtlsInt::trapzMul(nPts, xs, rho, dx);
		//Then find total potential
		//vtls::scaMulArray(surfPos, (-pol + dampRate/dt*(pol - lpol)) / wellWidth, xfsf, genTemp);
		vtls::scaMulArray(surfPos, totBulk, xfsf, genTemp);
		vtls::scaMulArray(surfPos, totBulk / wellWidth / 2.0, xfsf2, potTemp);
		vtls::addArrays(surfPos, genTemp, potTemp);
		//vtlsInt::cumIntTrapz(surfPos, genTemp, dx, potTemp);
		//Match potentials at surface, sans the expected difference due to field.
		if (surfPos < nPts && surfPos > 0) {
			double diff = potTemp[surfPos] - potTemp[surfPos - 1] - dx * totBulk;
			vtls::scaAddArray(nPts - surfPos, -diff, &potTemp[surfPos]);
		}
		//Recombine two potentials
		////vtls::addArrays(nPts, potTemp, targ);

		//Apply final constants
		vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe / PhysCon::e0, targ);
	}


	LinearBulkCylSectionFieldSpaceCharge::LinearBulkCylSectionFieldSpaceCharge(int nPts, double* x, double dx, double ef, double rad, double theta0, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int surfPos, int refPoint) {
		LinearBulkCylSectionFieldSpaceCharge::nPts = nPts;
		LinearBulkCylSectionFieldSpaceCharge::dx = dx;
		LinearBulkCylSectionFieldSpaceCharge::ef = ef;
		LinearBulkCylSectionFieldSpaceCharge::rad = rad;
		LinearBulkCylSectionFieldSpaceCharge::nelecPtr = nelec;
		LinearBulkCylSectionFieldSpaceCharge::totPot = totPot;
		LinearBulkCylSectionFieldSpaceCharge::surfPos = surfPos;
		LinearBulkCylSectionFieldSpaceCharge::refPoint = refPoint;
		LinearBulkCylSectionFieldSpaceCharge::wght = wght;
		LinearBulkCylSectionFieldSpaceCharge::dens = dens;

		int nrPt = vtls::findValue(nPts, x, -rad+x[surfPos]);
		if (nrPt == -1)
			nrPt = 0;

		mMat = new double[nPts * nPts];
		std::fill_n(mMat, nPts * nPts, 0.0);

		double pref = dx * dx / (2.0 * PhysCon::pi);

		for (int j = nrPt; j < nPts; j++) {
			double sum = 0.0;
			double leftIntegrand = pref * (-2 * rad * std::sin(theta0) / (rad + x[j] - x[surfPos]) - rad * std::sin(2.0 * theta0) / std::pow(rad + x[j] - x[surfPos], 2) * (rad + x[nrPt] - x[surfPos]));
			double rightIntegrand;
			for (int i = nrPt; i < nPts-1; i++) {
				if (i >= surfPos || j >= surfPos)
					mMat[i * nPts + j] = sum;
				if (i + 1 == j)
					rightIntegrand = pref * (rad / (rad + x[i + 1] - x[surfPos])) * (theta0);
				else
					rightIntegrand = pref * (rad / (rad + x[i+1] - x[surfPos])) * (theta0 + 2.0 * std::atan((2.0 * rad + x[i+1] + x[j] - 2.0 * x[surfPos]) / (x[i+1] - x[j]) * std::tan(theta0 / 2.0)));
				sum += 0.5 * (leftIntegrand + rightIntegrand);
				leftIntegrand = rightIntegrand;
			}
			mMat[(nPts-1)*nPts + j] = sum;
		}

		origPot = new double[nPts];
		rho = new double[nPts];
		genTemp = new double[nPts];
		potTemp = new double[nPts];
	}

	void LinearBulkCylSectionFieldSpaceCharge::negateGroundEffects(std::complex<double>* psi) {
		if (first) doFirst(psi);
		calcPot(psi, origPot);
	}

	void LinearBulkCylSectionFieldSpaceCharge::getV(double t, double* targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void LinearBulkCylSectionFieldSpaceCharge::getV(std::complex<double>* psi, double t, double* targ) {
		if (first) {
			doFirst(psi);
			std::cout << "LinearBulkCylSectionFieldSpaceCharge potential should be negated after initial state is found." << std::endl;
		}
		calcPot(psi, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	int LinearBulkCylSectionFieldSpaceCharge::isDynamic() {
		return 2;
	}

	void LinearBulkCylSectionFieldSpaceCharge::doFirst(std::complex<double>* psi) {
		first = 0;
		nelec = nelecPtr[0];
		psi2 = new double[nPts * nelec];
		prefactor = new double[nelec];
		double* energies = new double[nelec];
		double* v0w = new double[nPts];
		totPot->getV(0, v0w);
		WfcToRho::calcEnergies(nelec, nPts, dx, psi, v0w, energies);
		wght->calcWeights(nelec, energies, prefactor);
		delete[] energies, v0w;
		//calcFermiBoxDimensionalityConversion(nelec, nPts, dx, ef, psi, totPot, prefactor);
	}

	void LinearBulkCylSectionFieldSpaceCharge::calcPot(std::complex<double>* psi, double* targ) {
		//auto t1 = std::chrono::system_clock::now();
		if (first) doFirst(psi);
		//Caclulate 3-D density
		dens->calcRho(nPts, nelec, dx, prefactor, psi, rho);
		//auto t2 = std::chrono::system_clock::now();
		std::fill_n(targ, nPts, 0);
		std::fill_n(potTemp, nPts, 0);

		//CALCULATE FIELDS FROM BULK ELECTRONS
		//Calculate cumulative charge integral within metal
		vtlsInt::cumIntTrapz(surfPos, rho, dx, genTemp);
		//Calculate potential inside the metal from metallic electrons
		vtlsInt::cumIntTrapz(surfPos, genTemp, dx, potTemp);
		//Match potentials at surface, sans the expected difference due to field.
		double surfFieldElec = genTemp[surfPos - 1];

		cblas_dgemv(CblasRowMajor, CblasNoTrans, nPts, nPts, 1.0, mMat, nPts, rho, 1, 0.0, genTemp, 1);
		//vtls::mxv_plain_openmp(nPts, nPts, mMat, rho, genTemp);

		double surfField = surfFieldElec - (genTemp[surfPos + 1] - genTemp[surfPos] + genTemp[surfPos-1] - genTemp[surfPos - 2]) / (2.0 * dx);

		vtls::addArrays(nPts, genTemp, potTemp);

		if (surfPos < nPts && surfPos > 0) {
			double diff = potTemp[surfPos] - potTemp[surfPos - 1] - dx * surfField;
			vtls::scaAddArray(nPts - surfPos, -diff, &potTemp[surfPos]);
		}
		//Recombine two potentials
		vtls::addArrays(nPts, potTemp, targ);

		//Apply final constants
		vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe / PhysCon::e0, targ);
		//auto t3 = std::chrono::system_clock::now();

		//std::chrono::duration<double> dur1 = t2 - t1;
		//std::chrono::duration<double> dur2 = t3 - t2;
		//std::cout << dur1.count() << " " << dur2.count() << std::endl;
	}


	CompositePotential::CompositePotential(int nPts, int numSPots, int numDPots, int numWPots, Potential ** staticPots, Potential ** dynamicPots, Potential ** waveFuncDependentPots) {
		CompositePotential::nPts = nPts;
		CompositePotential::numSPots = numSPots;
		CompositePotential::numDPots = numDPots;
		CompositePotential::numWPots = numWPots;
		CompositePotential::staticPots = staticPots;
		CompositePotential::dynamicPots = dynamicPots;
		CompositePotential::waveFuncDependentPots = waveFuncDependentPots;
		v0 = new double[nPts];
		if (numSPots != 0)
			staticPots[0]->getV(0.0, v0);
		else
			std::fill_n(v0, nPts, 0.0);
		nv = new double[nPts];
		for (int i = 1; i < numSPots; i++) {
			staticPots[i]->getV(0.0, nv);
			vtls::addArrays(nPts, nv, v0);
		}
	}

	void CompositePotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v0, targ);
		for (int i = 0; i < numDPots; i++) {
			dynamicPots[i]->getV(t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
	}

	void CompositePotential::getV(std::complex<double> * psi, double t, double * targ) {
		vtls::copyArray(nPts, v0, targ);
		for (int i = 0; i < numDPots; i++) {
			dynamicPots[i]->getV(t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
		for (int i = 0; i < numWPots; i++) {
			waveFuncDependentPots[i]->getV(psi, t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
	}

	int CompositePotential::isDynamic() {
		if (numWPots > 0)
			return 2;
		else if (numDPots > 0)
			return 1;
		else
			return 0;
	}


	PotentialManager::PotentialManager(int nPts) {
		PotentialManager::nPts = nPts;
	}

	void PotentialManager::addPotential(Potential * pot) {
		int typ = pot->isDynamic();
		if (typ == 2) {
			waveFuncDependentPots.push_back(pot);
		}
		else if (typ == 1) {
			dynamicPots.push_back(pot);
		}
		else if (typ == 0) {
			staticPots.push_back(pot);
		}
		else {
			throw "Potential is neither dynamic nor static.";
		}
	}

	void PotentialManager::finishAddingPotentials() {
		int ns = staticPots.size();
		int nd = dynamicPots.size();
		int nw = waveFuncDependentPots.size();
		Potential ** spots = new Potential*[ns];
		Potential ** dpots = new Potential*[nd];
		Potential ** wpots = new Potential*[nw];
		for (int i = 0; i < ns; i++)
			spots[i] = staticPots[i];
		for (int i = 0; i < nd; i++)
			dpots[i] = dynamicPots[i];
		for (int i = 0; i < nw; i++)
			wpots[i] = waveFuncDependentPots[i];
		pot = new CompositePotential(nPts, ns, nd, nw, spots, dpots, wpots);
	}

	void PotentialManager::getV(double t, double * targ) {
		if (pot)
			pot->getV(t, targ);
		else {
			throw "Must run finishAddingPotentials() function on PotentialManager before using getV().";
		}
	}

	void PotentialManager::getV(std::complex<double> * psi, double t, double * targ) {
		if (pot)
			pot->getV(psi, t, targ);
		else {
			throw "Must run finishAddingPotentials() function on PotentialManager before using getV().";
		}
	}

	int PotentialManager::isDynamic() {
		return 2;
	}
}

namespace ElectricFieldProfiles {
	ConstantFieldProfile::ConstantFieldProfile(int nPts, double * x, double eMax, double minX, double maxX) {
		fs = new std::complex<double>[nPts];
		for (int i = 0; i < nPts; i++) {
			if (x[i] > minX && x[i] < maxX)
				fs[i] = eMax;
			else
				fs[i] = 0.0;
		}
	}
	std::complex<double> * ConstantFieldProfile::getProfile() {
		return fs;
	}

	CylindricalToLinearProfile::CylindricalToLinearProfile(int nPts, double * x, double minX, double maxX, double r, double eMax, double enhFact) {
		fs = new std::complex<double>[nPts];
		double xc = -enhFact * r + minX + std::sqrt((enhFact - 1.0)*r*(enhFact*r + maxX - minX));
		double xn;
		for (int i = 0; i < nPts; i++) {
			xn = x[i];
			if (xn < minX || xn > maxX)
				fs[i] = 0.0;
			else if (xn > minX && xn < xc)
				fs[i] = eMax / enhFact * ((enhFact - 1.0)*r / (xn - minX + r) + 1.0);
			else
				fs[i] = eMax / enhFact * ((enhFact - 1.0)*r / (xc - minX + r) + 1.0) * (maxX - xn) / (maxX - xc);
		}
	}

	std::complex<double> * CylindricalToLinearProfile::getProfile() {
		return fs;
	}

	CylindricalToCutoffProfile::CylindricalToCutoffProfile(int nPts, double * x, double minX, double maxX, double r, double eMax, double enhFact, double decayLength) {
		fs = new std::complex<double>[nPts];
		double xn, k;
		for (int i = 0; i < nPts; i++) {
			xn = x[i];
			if (xn < minX || xn > maxX)
				fs[i] = 0.0;
			else if (xn > minX && xn < maxX-decayLength)
				fs[i] = eMax / enhFact * ((enhFact - 1.0)*r / (xn - minX + r) + 1.0);
			else {
				k = -(x[i] - maxX) / decayLength;
				fs[i] = eMax / enhFact * ((enhFact - 1.0)*r / (xn - minX + r) + 1.0) * 
					(924.0*std::pow(k, 13) -
					6006.0*std::pow(k, 12) +
					16380.0*std::pow(k, 11) -
					24024.0*std::pow(k, 10) +
					20020.0*std::pow(k, 9) -
					9009.0*std::pow(k, 8) +
					1716.0*std::pow(k, 7));
			}
		}
	}

	std::complex<double> * CylindricalToCutoffProfile::getProfile() {
		return fs;
	}

	InMetalFieldProfile::InMetalFieldProfile(int nPts, double * x, double minX, double maxX, double eMax, double lam, std::complex<double> er, double cond) {
		fs = new std::complex<double>[nPts];
		double xn;
		//double w = PhysCon::c / lam * PhysCon::pi*2.0;
		std::complex<double> k, kx, kz;
		//plasmonic response https://en.wikipedia.org/wiki/Surface_plasmon_polariton#Propagation_length_and_skin_depth
		kx = 2.0*PhysCon::pi / lam * std::sqrt(er / (er + 1.0));
		kz = std::sqrt(er*std::pow((2 * PhysCon::pi / lam), 2.0) - kx * kx);
		//skin depth, slow frequencies
		//k = w * std::sqrt(er*PhysCon::e0*PhysCon::mu0 / 2.0)*std::sqrt(std::sqrt(1.0 + std::pow(cond / (er*PhysCon::e0*w), 2)) + 1.0);
		//k += PhysCon::im*w * std::sqrt(er*PhysCon::e0*PhysCon::mu0 / 2.0)*std::sqrt(std::sqrt(1.0 + std::pow(cond / (er*PhysCon::e0*w), 2)) - 1.0);
		for (int i = 0; i < nPts; i++) {
			xn = x[i];
			if (xn < minX || xn > maxX)
				fs[i] = 0.0;
			else {
				//plasmonic response
				fs[i] = -eMax * kx / kz * std::exp(PhysCon::im*(kz*std::abs(xn-maxX)));
				
				//skin depth, slow frequencies
				//fs[i] = 1.0 / er * eMax * std::exp(PhysCon::im*(k*(maxX - xn)));
				
			}
		}
	}

	std::complex<double> * InMetalFieldProfile::getProfile() {
		return fs;
	}

	FileFieldProfile::FileFieldProfile(int nPts, double * x, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, const char * fil) {
		fs = new std::complex<double>[nPts];
		double * tre = new double[nPts];
		double * tim = new double[nPts];
		std::fstream ifil = std::fstream(fil, std::ios::in | std::ios::binary);
		int nRep;
		ifil.read(reinterpret_cast<char*>(&nRep), sizeof(int));
		double * fx = new double[nRep];
		double * fre = new double[nRep];
		double * fim = new double[nRep];
		ifil.read(reinterpret_cast<char*>(fx), sizeof(double)*nRep);
		ifil.read(reinterpret_cast<char*>(fre), sizeof(double)*nRep);
		ifil.read(reinterpret_cast<char*>(fim), sizeof(double)*nRep);
		for (int i = 0; i < nRep; i++)
			fx[i] += offset;
		vtls::linearInterpolate(nRep, fx, fre, nPts, x, tre);
		vtls::linearInterpolate(nRep, fx, fim, nPts, x, tim);
		for (int i = 0; i < nPts; i++)
			fs[i] = (tre[i] + PhysCon::im*tim[i])*emax;
		double k;
		for (int i = 0; i < nPts; i++) {
			if (x[i] > leftDecayPos && x[i] < leftDecayPos + decayLength) {
				k = (x[i] - leftDecayPos) / decayLength;
				fs[i] *=
				924.0*std::pow(k, 13) -
					6006.0*std::pow(k, 12) +
					16380.0*std::pow(k, 11) -
					24024.0*std::pow(k, 10) +
					20020.0*std::pow(k, 9) -
					9009.0*std::pow(k, 8) +
					1716.0*std::pow(k, 7);
			}
			else if (x[i] < rightDecayPos && x[i] > rightDecayPos - decayLength) {
				k = -(x[i] - rightDecayPos) / decayLength;
				fs[i] *=
					924.0*std::pow(k, 13) -
					6006.0*std::pow(k, 12) +
					16380.0*std::pow(k, 11) -
					24024.0*std::pow(k, 10) +
					20020.0*std::pow(k, 9) -
					9009.0*std::pow(k, 8) +
					1716.0*std::pow(k, 7);
			}
			else if (x[i] < leftDecayPos || x[i] > rightDecayPos)
				fs[i] = 0;
		}
		ifil.close();
		delete[] tre;
		delete[] tim;
		delete[] fx;
		delete[] fre;
		delete[] fim;
	}

	std::complex<double> * FileFieldProfile::getProfile() {
		return fs;
	}
}

namespace Envelopes {
	GaussianEnvelope::GaussianEnvelope(double tau, double tmax) {
		GaussianEnvelope::tau = tau;
		GaussianEnvelope::tmax = tmax;
	}

	double GaussianEnvelope::getValue(double t) {
		return std::exp(-std::pow((t - tmax) / tau, 2) * 2.0 * std::log(2.0));
	}

	SmoothedInitialGaussianEnvelope::SmoothedInitialGaussianEnvelope(double tau, double tmax, double bufferTime) {
		SmoothedInitialGaussianEnvelope::tau = tau;
		SmoothedInitialGaussianEnvelope::tmax = tmax;
		SmoothedInitialGaussianEnvelope::buf = bufferTime;
	}

	double SmoothedInitialGaussianEnvelope::getValue(double t) {
		if (t < buf) {
			double k = t / buf;
			return std::exp(-std::pow((t - tmax) / tau, 2) * 2.0 * std::log(2.0)) * (
				924.0*std::pow(k, 13) -
				6006.0*std::pow(k, 12) +
				16380.0*std::pow(k, 11) -
				24024.0*std::pow(k, 10) +
				20020.0*std::pow(k, 9) -
				9009.0*std::pow(k, 8) +
				1716.0*std::pow(k, 7));
		} //NOTE: This is the 7th order smooth function -- IT WORKS! At higher orders we run into issues of floating point error.
		else {
			//return std::exp(-std::pow((t - tmax) / tau, 2) / 2.0);
			return std::pow(2.0, -2.0 * std::log(2.0)*std::pow((t - tmax) / tau, 2));
		}
	}
}