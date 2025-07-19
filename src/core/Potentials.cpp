#include "Potentials.h"
#include "Densities.h"
#include <bits/types/FILE.h>
#include <stdexcept>
#include "MathTools.h"
#include "PhysCon.h"
#include "blas.h"

namespace Potentials {
	FilePotential::FilePotential(size_t nPts, const double * x, double offset, const std::string fil, size_t refPoint) {
		FilePotential::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
		std::fstream ifil = std::fstream(fil, std::ios::in | std::ios::binary);
		int nRep;
		ifil.read(reinterpret_cast<char*>(&nRep), sizeof(int));

		double * fx = (double*) sq_malloc(sizeof(double)*nRep);
		double * fv = (double*) sq_malloc(sizeof(double)*nRep);
		ifil.read(reinterpret_cast<char*>(fx), sizeof(double)*nRep);
		ifil.read(reinterpret_cast<char*>(fv), sizeof(double)*nRep);
		for (size_t i = 0; i < nRep; i++)
			fx[i] += offset;
		vtls::linearInterpolate(nRep, fx, fv, nPts, x, v);
		vtls::scaAddArray(nPts, -v[refPoint], v);
		ifil.close();
		sq_free(fx);
		sq_free(fv);
	}

	FilePotential::~FilePotential(){
		sq_free(v);
	}

	void FilePotential::getVBare(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void FilePotential::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		vtls::copyArray(nPts, v, targ);
	}

	BiasFieldPotential::BiasFieldPotential(size_t nPts, const double * x, double tstart, double tbuf, double xmin, double xmax, double xmin_buf, double xmax_buf, double fieldStrength, size_t refPoint) {
		BiasFieldPotential::nPts = nPts;
		BiasFieldPotential::tstart = tstart;
		BiasFieldPotential::tbuf = tbuf;
		double dx = (x[nPts - 1] - x[0]) / nPts;

		v = (double*) sq_malloc(sizeof(double)*nPts);
		double cx;
		for (size_t i = 0; i < nPts; i++) {
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

	BiasFieldPotential::~BiasFieldPotential(){
		sq_free(v);
	}

	void BiasFieldPotential::getVBare(double t, double * targ) {
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

	void BiasFieldPotential::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVBare(t, targ);
	}

	CoulombPotential::CoulombPotential(size_t nPts, const double * x, double ne, double chargePos, double minX, double maxX, size_t refPoint) {
		v = (double*) sq_malloc(sizeof(double)*nPts);
		double k = PhysCon::qe*PhysCon::qe / (4.0*PhysCon::pi*PhysCon::e0);
		double dx = (x[nPts - 1] - x[0]) / nPts;
		for (size_t i = 0; i < nPts; i++) {
			if (x[i] < minX)
				v[i] = -k / std::abs(minX - chargePos);
			else if (x[i] > maxX)
				v[i] = -k / std::abs(maxX - chargePos);
			else
				v[i] = -k / std::max(std::abs(x[i] - chargePos), dx);
		}
		double ref = v[refPoint];
		for (size_t i = 0; i < nPts; i++)
			v[i] -= ref;
	}

	CoulombPotential::~CoulombPotential(){
		sq_free(v);
	}

	void CoulombPotential::getVBare(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void CoulombPotential::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		vtls::copyArray(nPts, v, targ);
	}

	FiniteBox::FiniteBox(size_t nPts, const double* x, double left, double right, double vin, size_t refPoint) {
		FiniteBox::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
		for (size_t i = 0; i < nPts; i++) {
			if (x[i] > left && x[i] < right)
				v[i] = vin;
			else
				v[i] = 0;
		}
		vtls::scaAddArray(nPts, -v[refPoint], v);
	}

	FiniteBox::~FiniteBox(){
		sq_free(v);
	}

	void FiniteBox::getVBare(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void FiniteBox::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVBare(t, targ);
	}


	JelliumPotential::JelliumPotential(size_t nPts, const double* x, double center, double ef, double w, size_t refPoint) {
		JelliumPotential::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
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
		for (size_t i = 0; i < nPts; i++) {
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

	JelliumPotential::~JelliumPotential(){
		sq_free(v);
	}

	void JelliumPotential::getVBare(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void JelliumPotential::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVBare(t, targ);
	}

	JelliumPotentialBacked::JelliumPotentialBacked(size_t nPts, const double* x, double center, double ef, double w, double backStart, double backWidth, size_t refPoint) {
		JelliumPotentialBacked::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
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
		for (size_t i = 0; i < nPts; i++) {
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
		for (size_t i = 0; x[i] < backStart + backWidth / 2 && i < nPts; i++) {
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

	JelliumPotentialBacked::~JelliumPotentialBacked(){
		sq_free(v);
	}

	void JelliumPotentialBacked::getVBare(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void JelliumPotentialBacked::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVBare(t, targ);
	}

	ElectricFieldProfileToPotential::ElectricFieldProfileToPotential(size_t nPts, ElectricFieldProfiles::ElectricFieldProfile * fieldProfile, double dx, double phase, double tmax, double lam, Envelopes::Envelope * env, size_t refPoint) {
		ElectricFieldProfileToPotential::tmax = tmax;
		ElectricFieldProfileToPotential::env = env;
		ElectricFieldProfileToPotential::nPts = nPts;
		ElectricFieldProfileToPotential::phase = phase;
		std::complex<double> * fieldMask = fieldProfile->getProfile();
		potMask = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
		vtlsInt::cumIntTrapz(nPts, fieldMask, -dx * PhysCon::qe, potMask);
		vtls::scaAddArray(nPts, -potMask[refPoint], potMask);
		w = PhysCon::c / lam * 2.0*PhysCon::pi;
	}

	ElectricFieldProfileToPotential::~ElectricFieldProfileToPotential(){
		sq_free(potMask);
	}

	void ElectricFieldProfileToPotential::getVBare(double t, double * targ) {
		vtls::scaMulArrayRe(nPts, std::exp(PhysCon::im*(w*(t-tmax)+phase))*env->getValue(t), potMask, targ);
	}

	void ElectricFieldProfileToPotential::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVBare(t, targ);
	}

	ShieldedAtomicPotential::ShieldedAtomicPotential(size_t nPts, const double* x, double center, double latticeSpacing, double zProtons, double decayLength) {
		ShieldedAtomicPotential::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
		using namespace PhysCon;
		for (size_t i = 0; i < nPts; i++)
			v[i] = -zProtons * qe*qe / (2 * e0*latticeSpacing*latticeSpacing / decayLength)*std::exp(-std::abs(x[i] - center) / decayLength);
	}

	ShieldedAtomicPotential::~ShieldedAtomicPotential(){
		sq_free(v);
	}

	void ShieldedAtomicPotential::getVBare(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void ShieldedAtomicPotential::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVBare(t, targ);
	}

	CurrentIntegrator::CurrentIntegrator(size_t nPts, double dx, size_t evalPoint, int side, const size_t* nElec, double * const * weights) :
		nPts(nPts), dx(dx), evalPoint(evalPoint), nElec(nElec), weights(weights), integratedFlux(0.0), tPrev(0.0), side(side) {};
	
	void CurrentIntegrator::integrate(const std::complex<double>* psi, double t) {
		size_t pt0;
		switch(side) {
		case 0: // central derivative
			for (size_t i = 0; i < *nElec; i++){
				pt0 = i * nPts + evalPoint;
				integratedFlux += (t - tPrev) * PhysCon::hbar / PhysCon::me * std::imag(std::conj(psi[pt0]) * \
					(psi[pt0 + 1] - psi[pt0 - 1]) / (2.0 * dx)) * (*weights)[i];
			}
			break;
		case 1: // right-side derivative
			for (size_t i = 0; i < *nElec; i++){
				pt0 = i * nPts + evalPoint;
				//integratedFlux += (t - tPrev) * PhysCon::hbar / PhysCon::me * std::imag(std::conj(psi[pt0]) * \
				//	(-psi[pt0 + 2] + 4.0*psi[pt0+1] - 3.0*psi[pt0]) / (2.0*dx)) * (*weights)[i];
				integratedFlux += (t - tPrev) * PhysCon::hbar / PhysCon::me * std::imag(std::conj(psi[pt0]) * \
					(psi[pt0 + 1] - psi[pt0]) / (dx)) * (*weights)[i];
			}
			break;
		case -1: // left-side derivative
			for (size_t i = 0; i < *nElec; i++){
				pt0 = i * nPts + evalPoint;
				//integratedFlux += (t - tPrev) * PhysCon::hbar / PhysCon::me * std::imag(std::conj(psi[pt0]) * \
				//	(3.0*psi[pt0] - 4.0*psi[pt0-1] + psi[pt0-2]) / (2.0*dx)) * (*weights)[i];
				integratedFlux += (t - tPrev) * PhysCon::hbar / PhysCon::me * std::imag(std::conj(psi[pt0]) * \
					(psi[pt0] - psi[pt0 - 1]) / (dx)) * (*weights)[i];
			}
			break;
		}

		tPrev = t;
	}

	CylindricalImageCharge::CylindricalImageCharge(size_t nPts, const double* x, double dx, double ef, double w, double rad, size_t surfPos, 
		const size_t* nElec, double * const * weights, const double* rho0, size_t posMin, size_t posMax, size_t refPoint) :
	 	nPts(nPts), dx(dx), ef(ef), w(w), rad(rad), refPoint(refPoint), x(x),
		posMin(posMin < 0 ? 0 : posMin),
		posMax(posMax > nPts - 1 ? nPts - 1 : posMax),
		surfPos(std::clamp(surfPos, (size_t)0, nPts - 1))
	{
		potTemp = (double*) sq_malloc(sizeof(double)*nPts);
		genTemp = (double*) sq_malloc(sizeof(double)*nPts);
		origPot = (double*) sq_malloc(sizeof(double)*nPts);
		myRho = (double*) sq_malloc(sizeof(double)*nPts);
		lrxr = (double*) sq_malloc(sizeof(double)*nPts);
		nsMask = (double*) sq_malloc(sizeof(double)*nPts);
		dethin = (double*) sq_malloc(sizeof(double)*nPts);

		curInt = new CurrentIntegrator(nPts, dx, posMax, -1, nElec, weights);

		for (size_t i = 0; i < nPts; i++)
			if (x[i] - x[surfPos] <= -rad)
				lrxr[i] = 0;
			else
				lrxr[i] = std::log((rad + x[i] - x[surfPos]) / rad);

		for (size_t i = 0; i < surfPos; i++)
			nsMask[i] = 0.0;
		for (size_t i = surfPos; i < nPts; i++)
			nsMask[i] = 1.0 - std::exp(-2 * std::sqrt(2.0 * PhysCon::me * w) / PhysCon::hbar * (i - surfPos) * dx);

		size_t temp1,temp2;
		Densities::CylindricalDensity::calcThinning(nPts, x[surfPos]-rad, rad, x[0], dx, dethin, &temp1, &temp2);
		for (size_t i = 0; i < nPts; i++)
			dethin[i] = 1.0/dethin[i];

		if (rho0 != nullptr)
			calcPot(rho0, nullptr, 0.0, origPot);
		else
			std::fill_n(origPot, nPts, 0.0);
	}

	CylindricalImageCharge::~CylindricalImageCharge(){
		sq_free(potTemp);
		sq_free(genTemp);
		sq_free(origPot);
		sq_free(myRho);
		sq_free(lrxr);
		sq_free(nsMask);
		sq_free(dethin);

		delete curInt;
	}

	void CylindricalImageCharge::getVBare(double t, double* targ) {
		for (size_t i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void CylindricalImageCharge::getVVirtual(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		calcPot(rho, psi, t, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (size_t i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	void CylindricalImageCharge::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVVirtual(rho, cur, psi, t, targ);
		curInt->integrate(psi, t);
	}

	void CylindricalImageCharge::calcPot(const double* rho, const std::complex<double>* psi, double cur_t, double* targ) {
		std::fill_n(targ, nPts, 0);

		vtls::seqMulArrays(nPts, dethin, rho, myRho);

		//CALCULATE FIELDS FROM VACUUM ELECTRONS
		//Calculate first integral
		vtlsInt::cumIntTrapz(nPts - surfPos, &myRho[surfPos], dx * rad, &targ[surfPos]);
		//Add in image charge
		vtls::scaAddArray(nPts - surfPos, -targ[nPts - 1] - curInt->getIntegratedFlux() * rad, &targ[surfPos]);
		vtls::seqMulArrays(nPts - surfPos, &lrxr[surfPos], &targ[surfPos]);
		//Calculate second integral
		std::fill_n(potTemp, nPts, 0);
		vtls::seqMulArrays(nPts - surfPos, &lrxr[surfPos], &myRho[surfPos], &genTemp[surfPos]);
		vtlsInt::cumIntTrapz(nPts - surfPos, &genTemp[surfPos], -dx * rad, &potTemp[surfPos]);
		vtls::addArrays(nPts - surfPos, &potTemp[surfPos], &targ[surfPos]);

		//Apply near-surface mask
		vtls::seqMulArrays(nPts, nsMask, targ);
		//Apply final constants
		vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe / PhysCon::e0, targ);
	}

	PlanarToCylindricalHartree::PlanarToCylindricalHartree(int mimicOpenSystem, int ghostCharge, size_t nPts, double dx, double rad, size_t surfPos,
		const size_t* nElec, double * const * weights, const double* rho0, size_t posMin, size_t posMax, size_t refPoint)  : 
		nPts(nPts), dx(dx), rad(rad), refPoint(refPoint),
		posMin(posMin < 0 ? 0 : posMin),
		posMax(posMax > nPts - 1 ? nPts - 1 : posMax),
		originalCharge(0.0),
		surfPos(std::clamp(surfPos, (size_t)0, nPts - 1)),
		mimicOpenSystem(mimicOpenSystem),
		ghostCharge(ghostCharge)
	{
		if(mimicOpenSystem != 0){
			curIntOpen = new CurrentIntegrator(nPts, dx, 
				mimicOpenSystem > 0 ? posMax : posMin,
				mimicOpenSystem > 0 ? (posMax == nPts - 1 ? -1 : 0) : (posMin == 0 ? 1 : 0),
				nElec, weights);
			if (mimicOpenSystem < 0)
				this->mimicOpenSystem = -1;
			if (mimicOpenSystem > 0)
				this->mimicOpenSystem = 1;
		}

		if(ghostCharge != 0){
			ghostPos = ghostCharge > 0 ? posMax : posMin;
			curIntGhost = new CurrentIntegrator(nPts, dx, 
				ghostPos,
				ghostCharge > 0 ? (posMax == nPts - 1 ? -1 : 0) : (posMin == 0 ? 1 : 0),
				nElec, weights);
			if (ghostCharge < 0)
				this->ghostCharge = -1;
			if (ghostCharge > 0)
				this->ghostCharge = 1;
		}

		potTemp = (double*) sq_malloc(sizeof(double)*nPts);
		origPot = (double*) sq_malloc(sizeof(double)*nPts);
		fieldScaler = (double*) sq_malloc(sizeof(double)*nPts);
		dethin = (double*) sq_malloc(sizeof(double)*nPts);
		myRho = (double*) sq_malloc(sizeof(double)*nPts);

		std::fill_n(potTemp, nPts, 0.0);
		std::fill_n(origPot, nPts, 0.0);
		std::fill_n(myRho, nPts, 0.0);

		// Calculate field scaler (R/z in vacuum, 1 in material) (z evaluated half a grid step to the right)
		for(size_t i = 0; i < nPts; i++)
			fieldScaler[i] = i >= surfPos ? rad / (rad + ((i-surfPos)+0.5)*dx) : 1.0;
		
		// Calculate dethin (1 in material, z/R in vacuum)
		size_t temp1, temp2;
		Densities::CylindricalDensity::calcThinning(nPts, dx*surfPos-rad, rad, 0.0, dx, dethin, &temp1, &temp2);
		for (size_t i = 0; i < nPts; i++)
			dethin[i] = 1.0 / dethin[i];

		if (rho0 != nullptr){
			calcPot(rho0, nullptr, 0.0, origPot);
			vtls::seqMulArrays(nPts, dethin, rho0, myRho);
			originalCharge = totalCharge;
		}
		else{
			std::fill_n(origPot, nPts, 0.0);
			std::fill_n(myRho, nPts, 0.0);
			originalCharge = 0.0;
		}
	}

	PlanarToCylindricalHartree::~PlanarToCylindricalHartree(){
		sq_free(potTemp);
		sq_free(origPot);
		sq_free(fieldScaler);
		sq_free(dethin);
		sq_free(myRho);

		if(curIntOpen)
			delete curIntOpen;
		if(curIntGhost)
			delete curIntGhost;
	}

	void PlanarToCylindricalHartree::getVBare(double t, double* targ) {
		std::fill_n(targ, nPts, 0.0);
	}

	void PlanarToCylindricalHartree::getVVirtual(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		calcPot(rho, psi, t, targ); // evaluates totalCharge as part of calculation
		if(mimicOpenSystem != 0){
			double lossFraction = -(originalCharge - totalCharge - mimicOpenSystem * curIntOpen->getIntegratedFlux()) / originalCharge + 1.0; // charge that left sim is lost, scale origPot by appropriate amount
			double ref = targ[refPoint] - lossFraction * origPot[refPoint];

			for (size_t i = 0; i < nPts; i++)
				targ[i] -= lossFraction * origPot[i] + ref;
		}
		else{
			double ref = targ[refPoint] - origPot[refPoint];
			for (size_t i = 0; i < nPts; i++)
				targ[i] -= origPot[i] + ref;
		}
	}

	void PlanarToCylindricalHartree::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		getVVirtual(rho, cur, psi, t, targ);
		if(curIntOpen)
			curIntOpen->integrate(psi, t);
		if(curIntGhost)
			curIntGhost->integrate(psi, t);
	}



	void PlanarToCylindricalHartree::calcPot(const double* rho, const std::complex<double>* psi, double t, double* targ) {
		std::fill_n(targ, posMin, 0);

		vtls::seqMulArrays(posMax-posMin, &dethin[posMin], &rho[posMin], &myRho[posMin]);
		if(curIntGhost)
			myRho[ghostPos] += ghostCharge * curIntGhost->getIntegratedFlux() / dx / dethin[ghostPos]; // add delta from ghost charge
		vtlsInt::cumIntTrapzToRight(nPts-posMin, &myRho[posMin], dx, &potTemp[posMin]); // cumulative integral of rho
		totalCharge = potTemp[nPts-1];
		vtls::seqMulArrays(nPts-posMin, &fieldScaler[posMin], &potTemp[posMin]); // scale by field scaler for 1/r term
		vtlsInt::cumIntTrapzToLeft(nPts-posMin, &potTemp[posMin], dx * -PhysCon::qe * PhysCon::qe / PhysCon::e0, &targ[posMin]); // final integral for potential, times constants
		//std::fill_n(&targ[posMax], nPts-posMax, targ[posMax-1]); // fill in right side with last value (zero field implied)
	}


	PlanarHartree::PlanarHartree(size_t nPts, double dx, const double* rho0, size_t refPoint) :
		nPts(nPts), dx(dx), refPoint(refPoint){
		origPot = (double*) sq_malloc(sizeof(double)*nPts);
		temp = (double*) sq_malloc(sizeof(double)*nPts);

		if(rho0 != nullptr)
			calcPot(rho0, origPot);
		else
			std::fill_n(origPot, nPts, 0.0);
	}

	PlanarHartree::~PlanarHartree(){
		sq_free(origPot);
		sq_free(temp);
	}

	void PlanarHartree::getVBare(double t, double* targ) {
		std::fill_n(targ, nPts, 0.0);
	}

	void PlanarHartree::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		calcPot(rho, targ);

		double ref = targ[refPoint] - origPot[refPoint];
		for (size_t i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	void PlanarHartree::calcPot(const double* rho, double* targ) {
		vtlsInt::cumIntTrapzToRight(nPts, rho, -PhysCon::qe * PhysCon::qe / PhysCon::e0*dx, temp);
		vtlsInt::cumIntTrapzToLeft(nPts, temp, dx, targ);
	}


	MixedGeometryHartreeGhostCharge::MixedGeometryHartreeGhostCharge(size_t nPts, size_t minPos, size_t maxPos, int gcSide, double dx, double mRTheta, const double* hRad, const double* rho0, const double* j0, size_t refPoint, bool includeVectorPotential) :
		nPts(nPts), minPos(minPos), maxPos(maxPos), gcSide(gcSide), dx(dx), mRTheta(mRTheta), refPoint(refPoint), includeVectorPotential(includeVectorPotential) 
	{
		assert(gcSide == 1 || gcSide == -1); // ghost charge must be either 1 or -1 here
		assert(nPts <= LAPACK_INT_MAX);

		gcPos = gcSide > 0 ? maxPos - 2 : minPos + 2;

		vld = (double*) sq_malloc(sizeof(double)*(nPts-1));
		vd  = (double*) sq_malloc(sizeof(double)*nPts);
		vud = (double*) sq_malloc(sizeof(double)*(nPts-1));
		vud2= (double*) sq_malloc(sizeof(double)*(nPts-2));
		vipiv=(lapack_int*) sq_malloc(sizeof(lapack_int)*nPts);
		vrhs= (double*) sq_malloc(sizeof(double)*nPts);
		newV= (double*) sq_malloc(sizeof(double)*nPts);
		drho= (double*) sq_malloc(sizeof(double)*nPts);
		dcur= (double*) sq_malloc(sizeof(double)*nPts);

		// fill matrix elements
		for (size_t i = 0; i < nPts-2; i++){
			vld[i] 		= 1.0 - (hRad[i+2] - hRad[i])/(4.0*hRad[i+1]);
			vud[i+1] 	= 1.0 + (hRad[i+2] - hRad[i])/(4.0*hRad[i+1]);
			vd[i+1] 	= -2.0 - dx*dx*mRTheta * mRTheta / (hRad[i+1]*hRad[i+1]);
		}

		// BCs (homogeneous, Neumann for GC side, Dirichlet for other side)
		if(gcSide == 1){
			vd[0] = 1.0; vud[0] = 0.0; // Dirichlet on left boundary
			vd[nPts-1] = 1.0; vld[nPts-2] = -1.0; // Neumann on right boundary
		}
		else{
			vd[0] = -1.0; vud[0] = 1.0; // Neumann on left boundary
			vd[nPts-1] = 1.0; vld[nPts-2] = 0.0; // Dirichlet on right boundary
		}

		// perform factorization (only needs to be done once! :) )
		lapack_int info;
		lapack_int nPtsL = static_cast<lapack_int>(nPts);
		LAPACK_dgttrf(&nPtsL, vld, vd, vud, vud2, vipiv, &info);
		if(info != 0)
			throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrf for potential failed with info = " + std::to_string(info));

		if(includeVectorPotential){
			ald = (double*) sq_malloc(sizeof(double)*(nPts-1));
			ad  = (double*) sq_malloc(sizeof(double)*nPts);
			aud = (double*) sq_malloc(sizeof(double)*(nPts-1));
			aud2= (double*) sq_malloc(sizeof(double)*(nPts-2));
			aipiv=(lapack_int*) sq_malloc(sizeof(lapack_int)*nPts);
			arhs= (double*) sq_malloc(sizeof(double)*nPts);
			newA= (double*) sq_malloc(sizeof(double)*nPts);
			oldV= (double*) sq_malloc(sizeof(double)*nPts);
			oldA= (double*) sq_malloc(sizeof(double)*nPts);
			aTemp= (double*) sq_malloc(sizeof(double)*nPts);
			oldVTrans = (double*) sq_malloc(sizeof(double)*nPts);

			// fill matrix elements (within simulation box, they are the same, only differ by BCs)
			for (size_t i = 0; i < nPts-2; i++){
				ald[i] = vld[i];
				aud[i+1] = vud[i+1];
				ad[i+1] = vd[i+1];
			}

			// BCs (always Neumann, but RHS will vary)
			ad[0] = -1.0; aud[0] = 1.0;
			ad[nPts-1] = 1.0; ald[nPts-2] = -1.0;

			// factor
			LAPACK_dgttrf(&nPtsL, ald, ad, aud, aud2, aipiv, &info);
			if(info != 0)
				throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrf for current failed with info = " + std::to_string(info));
		}

		if(rho0 != nullptr){ // include offset potential
			this->rho0 = (double*) sq_malloc(sizeof(double)*nPts);
			vtls::copyArray(nPts, rho0, this->rho0);
		}
		if(j0 != nullptr && includeVectorPotential){
			this->j0 = (double*) sq_malloc(sizeof(double)*nPts);
			vtls::copyArray(nPts, j0, this->j0);
		}
	}

	void MixedGeometryHartreeGhostCharge::calcPot(const double* rho, const double* cur, double* targ, double t, bool virt) {
		double dt = t - t0;

		// calculate differences if intial provided
		if(rho0)
			vtls::scaMulAddArrays(nPts, -1.0, rho0, rho, drho);
		else
			vtls::copyArray(nPts, rho, drho);

		if(j0)
			vtls::scaMulAddArrays(nPts, -1.0, j0, cur, dcur);
		else
			vtls::copyArray(nPts, cur, dcur);

		// electrostatic potential (newV will ultimately contain the preliminary result)
		// rhs
		vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe * dx * dx / PhysCon::e0, drho, newV);
		// set charge outside of system to zero
		std::fill_n(newV, minPos, 0.0);
		std::fill_n(&newV[maxPos+1], nPts-maxPos-1, 0.0);
		// add ghost charge
		newV[gcPos] += (ghostCharge + gcSide * dcur[gcPos]*dt) * (-PhysCon::qe * PhysCon::qe * dx / PhysCon::e0);
		newV[gcPos + gcSide * 1] = 0.0;
		newV[gcPos + gcSide * 2] = 0.0;
		// homogeneous BCs
		newV[0] = 0.0;
		newV[nPts-1] = 0.0;

		// solve
		lapack_int info, one=1;
		lapack_int nPtsL = static_cast<lapack_int>(nPts);
		LAPACK_dgttrs("N", &nPtsL, &one, vld, vd, vud, vud2, vipiv, newV, &nPtsL, &info);
		if(info != 0)
			throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrs for potential failed with info = " + std::to_string(info));

		// vector potential
		if(includeVectorPotential){
			size_t gcEdge = (gcSide == 1 ? maxPos : minPos);
			size_t diriEdge = (gcSide == 1 ? minPos : maxPos);

			vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe * dx * dx * PhysCon::mu0, dcur, newA);
			// set current outside of system to zero
			std::fill_n(newA, minPos, 0.0);
			std::fill_n(&newA[maxPos+1], nPts-maxPos-1, 0.0);

			// set BCs
			newA[diriEdge] = 0.0;
			if(!first)
				newA[gcEdge] = -oldAbDiff + gcSide * 2.0*dx/dt/PhysCon::c/PhysCon::c * (newV[gcEdge] - oldVb);
			else
				newA[gcEdge] = 0.0;

			// solve
			LAPACK_dgttrs("N", &nPtsL, &one, ald, ad, aud, aud2, aipiv, newA, &nPtsL, &info);
			if(info != 0)
				throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrs for current failed with info = " + std::to_string(info));

			// record new values for next iteration before gauge transformation
			if(!virt){
				oldVb = newV[gcEdge];
				oldAbDiff = newA[gcEdge] - newA[gcEdge - gcSide];
			}

			// apply gauge transformation, removing vector potential
			if(!first){
				vtlsInt::cumIntTrapz(nPts, newA, 2.0*dx/dt, targ);
				vtlsInt::cumIntTrapz(nPts, oldA,-2.0*dx/dt, aTemp);
				vtls::addArrays(nPts, aTemp, targ);
				vtls::addArrays(nPts, newV, targ);
				vtls::addArrays(nPts, oldV, targ);
				vtls::scaMulAddArrays(nPts, -1.0, oldVTrans, targ);
			}
			else{ // first step, no time derivative, gauge transformation is trivial
				vtls::copyArray(nPts, newV, targ);
			}

			// record whole potentials if not virtual step
			if(!virt){
				vtls::copyArray(nPts, newV, oldV);
				vtls::copyArray(nPts, newA, oldA);
				vtls::copyArray(nPts, targ, oldVTrans);
			}
		}
		else{
			vtls::copyArray(nPts, newV, targ);
		}

		// offset by reference point
		double ref = targ[refPoint];
		vtls::scaAddArray(nPts, -ref, targ);
		
		if(!virt){
			ghostCharge += gcSide * dcur[gcPos] * dt;
			first = false;
			t0 = t;
		}

		//vtlsPrnt::printArray(nPts, targ);
	}

	MixedGeometryHartreeGhostCharge::~MixedGeometryHartreeGhostCharge(){
		sq_free(vld);
		sq_free(vd);
		sq_free(vud);
		sq_free(vud2);
		sq_free(vipiv);
		sq_free(vrhs);
		sq_free(newV);
		if(rho0)
			sq_free(rho0);
		if(j0)
			sq_free(j0);
		sq_free(drho);
		sq_free(dcur);
		if(includeVectorPotential){
			sq_free(ald);
			sq_free(ad);
			sq_free(aud);
			sq_free(aud2);
			sq_free(aipiv);
			sq_free(arhs);
			sq_free(newA);
			sq_free(oldV);
			sq_free(oldA);
			sq_free(aTemp);
			sq_free(oldVTrans);
		}
	}



	MixedGeometryHartreeShielded::MixedGeometryHartreeShielded(size_t nPts, size_t minPos, size_t maxPos, size_t surfPos, double maskLength, double shieldLength, int neumannSide, double dx, double mRTheta, const double* hRad, const double* rho0, const double* j0, size_t refPoint, bool includeVectorPotential) :
		nPts(nPts), minPos(minPos), maxPos(maxPos), dx(dx), mRTheta(mRTheta), refPoint(refPoint), includeVectorPotential(includeVectorPotential), neumSide(neumannSide)
	{
		assert(nPts <= LAPACK_INT_MAX);

		vld = (double*) sq_malloc(sizeof(double)*(nPts-1));
		vd  = (double*) sq_malloc(sizeof(double)*nPts);
		vud = (double*) sq_malloc(sizeof(double)*(nPts-1));
		vud2= (double*) sq_malloc(sizeof(double)*(nPts-2));
		vipiv=(lapack_int*) sq_malloc(sizeof(lapack_int)*nPts);
		vrhs= (double*) sq_malloc(sizeof(double)*nPts);
		newV= (double*) sq_malloc(sizeof(double)*nPts);
		drho= (double*) sq_malloc(sizeof(double)*nPts);
		dcur= (double*) sq_malloc(sizeof(double)*nPts);
		shieldProfile = (double*) sq_malloc(sizeof(double)*nPts);
		maskProfile = (double*) sq_malloc(sizeof(double)*nPts);

		// create mask profile, use sigmoid according to maskLength
		vtls::masks::biSigmoid(nPts, minPos, maxPos, maskLength/dx, maskProfile);

		// build shield profile, decay to left for negative shieldLength, to right for positive shieldLength
		if(isnan(shieldLength) || isinf(shieldLength) || shieldLength == 0.0)
			useShielding = false;
		else{
			std::fill_n(shieldProfile, nPts, 1.0);
			if(shieldLength < 0.0)
				for(size_t i = 0; i < surfPos; i++)
					shieldProfile[i] = std::exp((surfPos - i) * dx / shieldLength);
			else if(shieldLength > 0.0)
				for(size_t i = surfPos; i < nPts; i++)
					shieldProfile[i] = std::exp((surfPos - i) * dx / shieldLength);
			useShielding = true;
		}

		// fill matrix elements
		for (size_t i = 0; i < nPts-2; i++){
			vld[i] 		= 1.0 - (hRad[i+2] - hRad[i])/(4.0*hRad[i+1]);
			vud[i+1] 	= 1.0 + (hRad[i+2] - hRad[i])/(4.0*hRad[i+1]);
			vd[i+1] 	= -2.0 - dx*dx*mRTheta*mRTheta / (hRad[i+1]*hRad[i+1]);
		}

		if(neumSide > 0){
			diriEdge = 0;
			neumEdge = nPts - 1;
			neumSide = 1;
		}
		else{
			diriEdge = nPts - 1;
			neumEdge = 0;
			neumSide = -1;
		}

		// BCs (homogeneous, Neumann for unshielded side, Dirichlet for other side
		if(neumSide == 1){
			vd[0] = 1.0; vud[0] = 0.0; // Dirichlet on left boundary
			vd[nPts-1] = 1.0; vld[nPts-2] = -1.0; // Neumann on right boundary
		}
		else{
			vd[0] = -1.0; vud[0] = 1.0; // Neumann on left boundary
			vd[nPts-1] = 1.0; vld[nPts-2] = 0.0; // Dirichlet on right boundary
		}

		// perform factorization (only needs to be done once! :) )
		lapack_int info;
		lapack_int nPtsL = static_cast<lapack_int>(nPts);
		LAPACK_dgttrf(&nPtsL, vld, vd, vud, vud2, vipiv, &info);
		if(info != 0)
			throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrf for potential failed with info = " + std::to_string(info));

		if(includeVectorPotential){
			ald = (double*) sq_malloc(sizeof(double)*(nPts-1));
			ad  = (double*) sq_malloc(sizeof(double)*nPts);
			aud = (double*) sq_malloc(sizeof(double)*(nPts-1));
			aud2= (double*) sq_malloc(sizeof(double)*(nPts-2));
			aipiv=(lapack_int*) sq_malloc(sizeof(lapack_int)*nPts);
			arhs= (double*) sq_malloc(sizeof(double)*nPts);
			newA= (double*) sq_malloc(sizeof(double)*nPts);
			oldV= (double*) sq_malloc(sizeof(double)*nPts);
			oldA= (double*) sq_malloc(sizeof(double)*nPts);
			aTemp= (double*) sq_malloc(sizeof(double)*nPts);
			oldVTrans = (double*) sq_malloc(sizeof(double)*nPts);

			// fill matrix elements (within simulation box, they are the same, only differ by BCs)
			for (size_t i = 0; i < nPts-2; i++){
				ald[i] = vld[i];
				aud[i+1] = vud[i+1];
				ad[i+1] = vd[i+1];
			}

			// BCs (always Neumann, but RHS will vary)
			ad[0] = -1.0; aud[0] = 1.0;
			ad[nPts-1] = 1.0; ald[nPts-2] = -1.0;

			// factor
			LAPACK_dgttrf(&nPtsL, ald, ad, aud, aud2, aipiv, &info);
			if(info != 0)
				throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrf for current failed with info = " + std::to_string(info));
		}

		if(rho0 != nullptr){ // include offset potential
			this->rho0 = (double*) sq_malloc(sizeof(double)*nPts);
			vtls::copyArray(nPts, rho0, this->rho0);
		}
		if(j0 != nullptr && includeVectorPotential){
			this->j0 = (double*) sq_malloc(sizeof(double)*nPts);
			vtls::copyArray(nPts, j0, this->j0);
		}
	}

	void MixedGeometryHartreeShielded::calcPot(const double* rho, const double* cur, double* targ, double t, bool virt) {
		double dt = t - t0;

		// calculate differences if intial provided
		if(rho0)
			vtls::scaMulAddArrays(nPts, -1.0, rho0, rho, drho);
		else
			vtls::copyArray(nPts, rho, drho);

		if(j0 && includeVectorPotential)
			vtls::scaMulAddArrays(nPts, -1.0, j0, cur, dcur);
		else if(includeVectorPotential)
			vtls::copyArray(nPts, cur, dcur);

		// electrostatic potential (newV will ultimately contain the preliminary result)
		// rhs
		vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe * dx * dx / PhysCon::e0, drho, newV);
		// apply charge mask profile
		vtls::seqMulArrays(nPts, maskProfile, newV);
			/*
			// set charge outside of system to zero
			std::fill_n(newV, minPos, 0.0);
			std::fill_n(&newV[maxPos+1], nPts-maxPos-1, 0.0);
			*/
		// homogeneous BCs
		newV[0] = 0.0;
		newV[nPts-1] = 0.0;

		// solve
		lapack_int info, one=1;
		lapack_int nPtsL = static_cast<lapack_int>(nPts);
		LAPACK_dgttrs("N", &nPtsL, &one, vld, vd, vud, vud2, vipiv, newV, &nPtsL, &info);
		if(info != 0)
			throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrs for potential failed with info = " + std::to_string(info));

		// vector potential
		if(includeVectorPotential){

			vtls::scaMulArray(nPts, -PhysCon::qe * PhysCon::qe * dx * dx * PhysCon::mu0, dcur, newA);
			// apply mask profile
			vtls::seqMulArrays(nPts, maskProfile, newA);
			/*
			// set current outside of system to zero
			std::fill_n(newA, minPos, 0.0);
			std::fill_n(&newA[maxPos+1], nPts-maxPos-1, 0.0);
			*/

			// set BCs
			newA[diriEdge] = 0.0;
			if(!first)
				newA[neumEdge] = -oldAbDiff + neumSide * 2.0*dx/dt/PhysCon::c/PhysCon::c * (newV[neumEdge] - oldVb);
			else
				newA[neumEdge] = 0.0;

			// solve
			LAPACK_dgttrs("N", &nPtsL, &one, ald, ad, aud, aud2, aipiv, newA, &nPtsL, &info);
			if(info != 0)
				throw std::runtime_error("MixedGeometryHartree: LAPACK_dgttrs for current failed with info = " + std::to_string(info));

			// record new values for next iteration before gauge transformation
			if(!virt){
				oldVb = newV[neumEdge];
				oldAbDiff = newA[neumEdge] - newA[neumEdge - neumSide];
			}

			// apply gauge transformation, removing vector potential
			if(!first){
				vtlsInt::cumIntTrapz(nPts, newA, 2.0*dx/dt, targ);
				vtlsInt::cumIntTrapz(nPts, oldA,-2.0*dx/dt, aTemp);
				vtls::addArrays(nPts, aTemp, targ);
				vtls::addArrays(nPts, newV, targ);
				vtls::addArrays(nPts, oldV, targ);
				vtls::scaMulAddArrays(nPts, -1.0, oldVTrans, targ);
			}
			else{ // first step, no time derivative, gauge transformation is trivial
				vtls::copyArray(nPts, newV, targ);
			}

			// record whole potentials if not virtual step
			if(!virt){
				vtls::copyArray(nPts, newV, oldV);
				vtls::copyArray(nPts, newA, oldA);
				vtls::copyArray(nPts, targ, oldVTrans);
			}
		}
		else{
			vtls::copyArray(nPts, newV, targ);
		}

		// apply shield profile to derivative if needed
		if(useShielding){
			vtls::firstDerivative(nPts, targ, aTemp, 1.0);
			vtls::seqMulArrays(nPts, shieldProfile, aTemp);
			vtlsInt::cumIntTrapz(nPts, aTemp, 1.0, targ);
	}

		// offset by reference point
		double ref = targ[refPoint];
		vtls::scaAddArray(nPts, -ref, targ);
		
		if(!virt){
			first = false;
			t0 = t;
		}

		//vtlsPrnt::printArray(nPts, targ);
	}

	MixedGeometryHartreeShielded::~MixedGeometryHartreeShielded(){
		sq_free(vld);
		sq_free(vd);
		sq_free(vud);
		sq_free(vud2);
		sq_free(vipiv);
		sq_free(vrhs);
		sq_free(newV);
		sq_free(shieldProfile);
		sq_free(maskProfile);
		if(rho0)
			sq_free(rho0);
		if(j0)
			sq_free(j0);
		sq_free(drho);
		sq_free(dcur);
		if(includeVectorPotential){
			sq_free(ald);
			sq_free(ad);
			sq_free(aud);
			sq_free(aud2);
			sq_free(aipiv);
			sq_free(arhs);
			sq_free(newA);
			sq_free(oldV);
			sq_free(oldA);
			sq_free(aTemp);
			sq_free(oldVTrans);
		}
	}



	LDAFunctional::LDAFunctional(LDAFunctionalType typ, size_t nPts, double dx, const double* rho0, size_t refPoint)
	: typ(typ), nPts(nPts), dx(dx), refPoint(refPoint) {
		origPot = (double*) sq_malloc(sizeof(double)*nPts);
		std::fill_n(origPot, nPts, 0.0);
		rho = (double*) sq_malloc(sizeof(double)*nPts);

		if (rho0 != nullptr)
			calcPot(rho0, origPot);
		else
			std::fill_n(origPot, nPts, 0.0);
	}

	LDAFunctional::~LDAFunctional(){
		sq_free(origPot);
		sq_free(rho);
	};

	void LDAFunctional::getVBare(double t, double* targ) {
		std::fill_n(targ, nPts, 0.0);
	}

	void LDAFunctional::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		calcPot(rho, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (size_t i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	void LDAFunctional::calcPot(const double* rho, double* targ) {
		switch (typ) {
		case LDAFunctionalType::X_SLATER: // slater exchange
		{
			double coef = -std::pow(3.0/PhysCon::pi, 1.0/3) * PhysCon::auE_ha * PhysCon::a0; //convert linear density to a.u., then to energy in SI
			for (size_t i = 0; i < nPts; i++)
				targ[i] = coef * std::pow(rho[i], 1.0/3);
			break;
		}
		case LDAFunctionalType::C_PW: // PW correlation
		{
			double aa = 0.031091, al = 0.21370, be1 = 7.5957, be2 = 3.5876, be3 = 1.6382, be4 = 0.49294;
			
			double crs, crho, q0, q1, q1p, drsdrho;
			double smallRho = 1e-10;
			for(size_t i = 0; i < nPts; i++){
				crho = rho[i] * std::pow(PhysCon::a0,3);
				if(crho < smallRho){
					targ[i] = PhysCon::auE_ha * al/be4*std::pow(4.0*PhysCon::pi/3.0 * crho, 2.0/3);
				} else{
					crs = std::pow(0.75/(PhysCon::pi * crho), 1.0/3) / PhysCon::a0;
					drsdrho = - crs / (3.0 * crho);
					q0 = -2*aa*(1+al*crs);
					q1 = 2*aa*(be1*std::sqrt(crs) + be2*crs + be3*std::pow(crs, 3.0/2) + be4*std::pow(crs, 2));
					q1p= aa*(be1/std::sqrt(crs) + 2*be2 + 3*be3*std::sqrt(crs) + 4*be4*crs);

					targ[i] = PhysCon::auE_ha * (-2*aa*al*std::log(1.0+1.0/q1) - q0*q1p/(q1*(q1+1.0))) * drsdrho;
				}
			}
			break;
		}
		}

	}

	CompositePotential::CompositePotential(size_t nPts, size_t numSPots, size_t numDPots, Potential ** staticPots, Potential ** dynamicPots) :
		nPts(nPts), numSPots(numSPots), numDPots(numDPots), staticPots(staticPots), dynamicPots(dynamicPots)
	{
		// calculate static part of potential
		v0 = (double*) sq_malloc(sizeof(double)*nPts);
		if (numSPots != 0)
			staticPots[0]->getVBare(0.0, v0);
		else
			std::fill_n(v0, nPts, 0.0);
		nv = (double*) sq_malloc(sizeof(double)*nPts);
		for (size_t i = 1; i < numSPots; i++) {
			staticPots[i]->getVBare(0.0, nv);
			vtls::addArrays(nPts, nv, v0);
		}

		// calculate complexity
		for (size_t i = 0; i < numDPots; i++)
			myDepend |= dynamicPots[i]->getDependence();
	}

	CompositePotential::~CompositePotential(){
		sq_free(v0);
		sq_free(nv);
	}

	void CompositePotential::getVBare(double t, double * targ) {
		vtls::copyArray(nPts, v0, targ);
		for (size_t i = 0; i < numDPots; i++) {
			dynamicPots[i]->getVBare(t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
	}

	void CompositePotential::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		vtls::copyArray(nPts, v0, targ);
		for (size_t i = 0; i < numDPots; i++) {
			dynamicPots[i]->getV(rho, cur, psi, t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
	}

	void CompositePotential::getVVirtual(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		vtls::copyArray(nPts, v0, targ);
		for (size_t i = 0; i < numDPots; i++) {
			dynamicPots[i]->getVVirtual(rho, cur, psi, t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
	}

	Dependence CompositePotential::getDependence() const {
		return myDepend;
	}

	PotentialManager::PotentialManager(size_t nPts) {
		PotentialManager::nPts = nPts;
	}

	void PotentialManager::addPotential(Potential * pot) {
		compositeRefreshed = false;
		if (pot -> getDependence() == Dependence::NONE)
			staticPots.push_back(pot);
		else{
			dynamicPots.push_back(pot);
			myDepend |= pot->getDependence();
		}
	}

	void PotentialManager::refreshCompositePotential() {
		// copy vectors to pointer arrays
		size_t ns = staticPots.size();
		size_t nd = dynamicPots.size();

		if (spots)
			delete[] spots;
		if (dpots)
			delete[] dpots;

		spots = new Potential*[ns > 0 ? ns : 1];
		dpots = new Potential*[nd > 0 ? nd : 1];

		for (size_t i = 0; i < ns; i++)
			spots[i] = staticPots[i];
		for (size_t i = 0; i < nd; i++)
			dpots[i] = dynamicPots[i];
		
		// recreate composite potential
		if (pot)
			delete pot;
		pot = new CompositePotential(nPts, ns, nd, spots, dpots);

		compositeRefreshed = true;
	}

	void PotentialManager::getVBare(double t, double * targ) {
		if(!compositeRefreshed)
			refreshCompositePotential();
		pot->getVBare(t, targ);
	}

	void PotentialManager::getV(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		if(!compositeRefreshed)
			refreshCompositePotential();
		pot->getV(rho, cur, psi, t, targ);
	}

	void PotentialManager::getVVirtual(const double* rho, const double* cur, const std::complex<double>* psi, double t, double* targ) {
		if(!compositeRefreshed)
			refreshCompositePotential();
		pot->getVVirtual(rho, cur, psi, t, targ);
	}

	namespace ElectricFieldProfiles {
		ConstantFieldProfile::ConstantFieldProfile(size_t nPts, const double * x, double eMax, double minX, double maxX) : ElectricFieldProfile(nPts) {
			for (size_t i = 0; i < nPts; i++) {
				if (x[i] > minX && x[i] < maxX)
					fs[i] = eMax;
				else
					fs[i] = 0.0;
			}
		}

		CylindricalToLinearProfile::CylindricalToLinearProfile(size_t nPts, const double * x, double minX, double maxX, double r, double eMax, double enhFact) : ElectricFieldProfile(nPts) {
			double xc = -enhFact * r + minX + std::sqrt((enhFact - 1.0)*r*(enhFact*r + maxX - minX));
			double xn;
			for (size_t i = 0; i < nPts; i++) {
				xn = x[i];
				if (xn < minX || xn > maxX)
					fs[i] = 0.0;
				else if (xn > minX && xn < xc)
					fs[i] = eMax / enhFact * ((enhFact - 1.0)*r / (xn - minX + r) + 1.0);
				else
					fs[i] = eMax / enhFact * ((enhFact - 1.0)*r / (xc - minX + r) + 1.0) * (maxX - xn) / (maxX - xc);
			}
		}

		CylindricalToCutoffProfile::CylindricalToCutoffProfile(size_t nPts, const double * x, double minX, double maxX, double r, double eMax, double enhFact, double decayLength) : ElectricFieldProfile(nPts) {
			double xn, k;
			for (size_t i = 0; i < nPts; i++) {
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

		InternalPlasmonicFieldProfile::InternalPlasmonicFieldProfile(size_t nPts, const double * x, double minX, double maxX, double eMax, double lam, std::complex<double> er, double cond) : ElectricFieldProfile(nPts) {
			double xn;
			//double w = PhysCon::c / lam * PhysCon::pi*2.0;
			std::complex<double> k, kx, kz;
			//plasmonic response https://en.wikipedia.org/wiki/Surface_plasmon_polariton#Propagation_length_and_skin_depth
			kx = 2.0*PhysCon::pi / lam * std::sqrt(er / (er + 1.0));
			kz = std::sqrt(er*std::pow((2 * PhysCon::pi / lam), 2.0) - kx * kx);
			//skin depth, slow frequencies
			//k = w * std::sqrt(er*PhysCon::e0*PhysCon::mu0 / 2.0)*std::sqrt(std::sqrt(1.0 + std::pow(cond / (er*PhysCon::e0*w), 2)) + 1.0);
			//k += PhysCon::im*w * std::sqrt(er*PhysCon::e0*PhysCon::mu0 / 2.0)*std::sqrt(std::sqrt(1.0 + std::pow(cond / (er*PhysCon::e0*w), 2)) - 1.0);
			for (size_t i = 0; i < nPts; i++) {
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

		FileFieldProfile::FileFieldProfile(size_t nPts, const double * x, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, const std::string fil) : ElectricFieldProfile(nPts) {
			double * tre = (double*) sq_malloc(sizeof(double)*nPts);
			double * tim = (double*) sq_malloc(sizeof(double)*nPts);
			std::fstream ifil = std::fstream(fil, std::ios::in | std::ios::binary);
			int nRep;
			ifil.read(reinterpret_cast<char*>(&nRep), sizeof(int));
			double * fx = (double*) sq_malloc(sizeof(double)*nRep);
			double * fre = (double*) sq_malloc(sizeof(double)*nRep);
			double * fim = (double*) sq_malloc(sizeof(double)*nRep);
			ifil.read(reinterpret_cast<char*>(fx), sizeof(double)*nRep);
			ifil.read(reinterpret_cast<char*>(fre), sizeof(double)*nRep);
			ifil.read(reinterpret_cast<char*>(fim), sizeof(double)*nRep);
			for (size_t i = 0; i < nRep; i++)
				fx[i] += offset;
			vtls::linearInterpolate(nRep, fx, fre, nPts, x, tre);
			vtls::linearInterpolate(nRep, fx, fim, nPts, x, tim);
			for (size_t i = 0; i < nPts; i++)
				fs[i] = (tre[i] + PhysCon::im*tim[i])*emax;
			double k;
			for (size_t i = 0; i < nPts; i++) {
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
			sq_free(tre);
			sq_free(tim);
			sq_free(fx);
			sq_free(fre);
			sq_free(fim);
		}

		ExponentialToLinearProfile::ExponentialToLinearProfile(size_t nPts, const double* x, double minX, double maxX, double r, double eMax) : ElectricFieldProfile(nPts) {
			double xn;
			for (size_t i = 0; i < nPts; i++) {
				xn = x[i];
				if (xn < minX || xn > maxX)
					fs[i] = 0.0;
				else if (xn < maxX - r)
					fs[i] = eMax * std::exp(-(xn - minX) / r);
				else
					fs[i] = eMax / r * std::exp(-(maxX - r) / r) * (maxX - xn);
			}
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

		CosSquaredEnvelope::CosSquaredEnvelope(double tau, double tmax) {
			CosSquaredEnvelope::tau = tau;
			CosSquaredEnvelope::tmax = tmax;
		}

		double CosSquaredEnvelope::getValue(double t) {
			if (std::abs((t - tmax) * a_t / tau) < PhysCon::pi / 2.0)
				return std::pow(std::cos(a_t * (t - tmax) / tau), 2);
			else
				return 0.0;
		}
	}
}