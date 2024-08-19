#include "Potentials.h"
#include <bits/types/FILE.h>
#include <stdexcept>
#include "MathTools.h"
#include "PhysCon.h"
#include "blas_externs.h"

namespace Potentials {
	FilePotential::FilePotential(int nPts, double * x, double offset, const char * fil, int refPoint) {
		FilePotential::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
		std::fstream ifil = std::fstream(fil, std::ios::in | std::ios::binary);
		int nRep;
		ifil.read(reinterpret_cast<char*>(&nRep), sizeof(int));

		double * fx = (double*) sq_malloc(sizeof(double)*nRep);
		double * fv = (double*) sq_malloc(sizeof(double)*nRep);
		ifil.read(reinterpret_cast<char*>(fx), sizeof(double)*nRep);
		ifil.read(reinterpret_cast<char*>(fv), sizeof(double)*nRep);
		for (int i = 0; i < nRep; i++)
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

	void FilePotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void FilePotential::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	BiasFieldPotential::BiasFieldPotential(int nPts, double * x, double tstart, double tbuf, double xmin, double xmax, double xmin_buf, double xmax_buf, double fieldStrength, int refPoint) {
		BiasFieldPotential::nPts = nPts;
		BiasFieldPotential::tstart = tstart;
		BiasFieldPotential::tbuf = tbuf;
		double dx = (x[nPts - 1] - x[0]) / nPts;

		v = (double*) sq_malloc(sizeof(double)*nPts);
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

	BiasFieldPotential::~BiasFieldPotential(){
		sq_free(v);
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

	void BiasFieldPotential::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	CoulombPotential::CoulombPotential(int nPts, double * x, double ne, double chargePos, double minX, double maxX, int refPoint) {
		v = (double*) sq_malloc(sizeof(double)*nPts);
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

	CoulombPotential::~CoulombPotential(){
		sq_free(v);
	}

	void CoulombPotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void CoulombPotential::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	FiniteBox::FiniteBox(int nPts, double * x, double left, double right, double vin, int refPoint) {
		FiniteBox::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
		for (int i = 0; i < nPts; i++) {
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

	void FiniteBox::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void FiniteBox::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}


	JelliumPotential::JelliumPotential(int nPts, double * x, double center, double ef, double w, int refPoint) {
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

	JelliumPotential::~JelliumPotential(){
		sq_free(v);
	}

	void JelliumPotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void JelliumPotential::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	JelliumPotentialBacked::JelliumPotentialBacked(int nPts, double * x, double center, double ef, double w, double backStart, double backWidth, int refPoint) {
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

	JelliumPotentialBacked::~JelliumPotentialBacked(){
		sq_free(v);
	}

	void JelliumPotentialBacked::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void JelliumPotentialBacked::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	ElectricFieldProfileToPotential::ElectricFieldProfileToPotential(int nPts, ElectricFieldProfiles::ElectricFieldProfile * fieldProfile, double dx, double phase, double tmax, double lam, Envelopes::Envelope * env, int refPoint) {
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

	void ElectricFieldProfileToPotential::getV(double t, double * targ) {
		vtls::scaMulArrayRe(nPts, std::exp(PhysCon::im*(w*(t-tmax)+phase))*env->getValue(t), potMask, targ);
	}

	void ElectricFieldProfileToPotential::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	ShieldedAtomicPotential::ShieldedAtomicPotential(int nPts, double * x, double center, double latticeSpacing, double zProtons, double decayConst) {
		ShieldedAtomicPotential::nPts = nPts;
		v = (double*) sq_malloc(sizeof(double)*nPts);
		using namespace PhysCon;
		for (int i = 0; i < nPts; i++)
			v[i] = -zProtons * qe*qe / (2 * e0*latticeSpacing*latticeSpacing / decayConst)*std::exp(-std::abs(x[i] - center) / decayConst);
	}

	ShieldedAtomicPotential::~ShieldedAtomicPotential(){
		sq_free(v);
	}

	void ShieldedAtomicPotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v, targ);
	}

	void ShieldedAtomicPotential::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		getV(t, targ);
	}

	CurrentIntegrator::CurrentIntegrator(int nPts, double dx, int evalPoint,  int* nElec, double** weights) :
		nPts(nPts), dx(dx), evalPoint(evalPoint), nElec(nElec), weights(weights), integratedFlux(0.0), last_t(0.0) {};
	
	void CurrentIntegrator::integrate(std::complex<double>* psi, double t) {
		for (int i = 0; i < *nElec; i++)
			integratedFlux += (t - last_t) * PhysCon::hbar / PhysCon::me * std::imag(std::conj(psi[i * nPts + evalPoint]) * (psi[i * nPts + evalPoint + 1] - psi[i * nPts + evalPoint - 1]) / (2.0 * dx)) * (*weights)[i];
		last_t = t;
	}

	CylindricalImageCharge::CylindricalImageCharge(int nPts, double* x, double dx, double ef, double w, double rad, int* nElec,
	 double** weights, int posMin, int posMax, int surfPos, int refPoint) :
	 	nPts(nPts), dx(dx), ef(ef), w(w), rad(rad), refPoint(refPoint),
		posMin(posMin < 0 ? 0 : posMin),
		posMax(posMax > nPts - 1 ? nPts - 1 : posMax),
		surfPos(surfPos > nPts - 1 ? nPts - 1 : surfPos)
	  {
		potTemp = (double*) sq_malloc(sizeof(double)*nPts);
		genTemp = (double*) sq_malloc(sizeof(double)*nPts);
		origPot = (double*) sq_malloc(sizeof(double)*nPts);
		myRho = (double*) sq_malloc(sizeof(double)*nPts);
		lrxr = (double*) sq_malloc(sizeof(double)*nPts);
		nsMask = (double*) sq_malloc(sizeof(double)*nPts);
		dethin = (double*) sq_malloc(sizeof(double)*nPts);

		for (int i = 0; i < nPts; i++)
			if (x[i] - x[surfPos] <= -rad)
				lrxr[i] = 0;
			else
				lrxr[i] = std::log((rad + x[i] - x[surfPos]) / rad);

		for (int i = 0; i < surfPos; i++)
			nsMask[i] = 0.0;
		for (int i = surfPos; i < nPts; i++)
			nsMask[i] = 1.0 - std::exp(-2 * std::sqrt(2.0 * PhysCon::me * w) / PhysCon::hbar * (i - surfPos) * dx);

		for (int i = 0; i < nPts; i++)
			dethin[i] = i >= surfPos ? 1.0 + (i-surfPos)*dx/rad : 1.0;

		curInt = new CurrentIntegrator(nPts, dx, posMax, nElec, weights);
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
	
	void CylindricalImageCharge::negateGroundEffects(double* rho, std::complex<double>* psi) {
		calcPot(rho, psi, 0.0, origPot);
	}

	void CylindricalImageCharge::getV(double t, double* targ) {
		for (int i = 0; i < nPts; i++)
			targ[i] = 0.0;
	}

	void CylindricalImageCharge::getV(double* rho, std::complex<double>* psi, double t, double* targ) {
		calcPot(rho, psi, t, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	void CylindricalImageCharge::calcPot(double* rho, std::complex<double>* psi, double cur_t, double* targ) {
		curInt->integrate(psi, cur_t);

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

	PlanarToCylindricalHartree::PlanarToCylindricalHartree(int nPts, double* x, double dx, double rad, int* nElec, double** weights, int posMin, int posMax, int surfPos, int refPoint) : 
		nPts(nPts), dx(dx), rad(rad), refPoint(refPoint),
		posMin(posMin < 0 ? 0 : posMin),
		posMax(posMax > nPts - 1 ? nPts - 1 : posMax),
		surfPos(surfPos > nPts - 1 ? nPts - 1 : surfPos),
		originalCharge(0.0)
	{
		curInt = new CurrentIntegrator(nPts, dx, posMin, nElec, weights);

		potTemp = (double*) sq_malloc(sizeof(double)*nPts);
		origPot = (double*) sq_malloc(sizeof(double)*nPts);
		fieldScaler = (double*) sq_malloc(sizeof(double)*nPts);
		dethin = (double*) sq_malloc(sizeof(double)*nPts);
		myRho = (double*) sq_malloc(sizeof(double)*nPts);

		std::fill_n(potTemp, nPts, 0.0);
		std::fill_n(origPot, nPts, 0.0);

		// Calculate field scaler (R/z in vacuum, 1 in material) (z evaluated half a grid step to the right)
		for(int i = 0; i < nPts; i++)
			fieldScaler[i] = i >= surfPos ? rad / (rad + ((i-surfPos)+0.5)*dx) : 1.0;
		
		// Calculate dethin (1 in material, z/R in vacuum)
		for (int i = 0; i < nPts; i++)
			dethin[i] = i >= surfPos ? 1.0 + (i-surfPos)*dx/rad : 1.0;
	}

	PlanarToCylindricalHartree::~PlanarToCylindricalHartree(){
		sq_free(potTemp);
		sq_free(origPot);
		sq_free(fieldScaler);
		sq_free(dethin);
		sq_free(myRho);

		delete curInt;
	}

	void PlanarToCylindricalHartree::negateGroundEffects(double* rho, std::complex<double>* psi) {
		calcPot(rho, psi, 0.0, origPot);
		vtls::seqMulArrays(nPts, dethin, rho, myRho);
		originalCharge = vtlsInt::trapz(nPts, myRho, dx);
	}

	void PlanarToCylindricalHartree::getV(double t, double* targ) {
		std::fill_n(targ, nPts, 0.0);
	}

	void PlanarToCylindricalHartree::getV(double* rho, std::complex<double>* psi, double t, double* targ) {
		if(originalCharge == 0.0)
			throw std::runtime_error("PlanarToCylindricalHartree::getV called before original charge was set (call negateGroundEffects first or there was zero net charge, somehow)");

		calcPot(rho, psi, t, targ);
		double lossFraction = -curInt->getIntegratedFlux() / originalCharge + 1.0; // charge that moved to left is lost, scale origPot by appropriate amount
		double ref = targ[refPoint] - lossFraction * origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= lossFraction * origPot[i] + ref;
	}

	void PlanarToCylindricalHartree::calcPot(double* rho, std::complex<double>* psi, double t, double* targ) {
		curInt->integrate(psi, t);

		std::fill_n(targ, nPts, 0);

		vtls::seqMulArrays(posMax-posMin, &dethin[posMin], &rho[posMin], &myRho[posMin]);
		vtlsInt::cumIntTrapzToRight(posMax-posMin, &myRho[posMin], dx, potTemp); // cumulative integral of rho
		vtls::seqMulArrays(posMax-posMin, &potTemp[posMin], &fieldScaler[posMin], &potTemp[posMin]); // scale by field scaler for 1/r term
		vtlsInt::cumIntTrapzToLeft(posMax-posMin, &potTemp[posMin], dx * -PhysCon::qe * PhysCon::qe / PhysCon::e0, &targ[posMin]); // final integral for potential, times constants
		std::fill_n(&targ[posMax], nPts-posMax, targ[posMax-1]); // fill in right side with last value (zero field implied)
	}

	LDAFunctional::LDAFunctional(LDAFunctionalType typ, int nPts, double dx, int refPoint)
	: typ(typ), nPts(nPts), dx(dx), refPoint(refPoint) {
		origPot = (double*) sq_malloc(sizeof(double)*nPts);
		std::fill_n(origPot, nPts, 0.0);
		rho = (double*) sq_malloc(sizeof(double)*nPts);
	}

	LDAFunctional::~LDAFunctional(){
		sq_free(origPot);
		sq_free(rho);
	};

	void LDAFunctional::negateGroundEffects(double* rho, std::complex<double>* psi) {
		calcPot(rho, origPot);
	}

	void LDAFunctional::getV(double t, double* targ) {
		std::fill_n(targ, nPts, 0.0);
	}

	void LDAFunctional::getV(double* rho, std::complex<double>* psi, double t, double* targ) {
		calcPot(rho, targ);
		double ref = targ[refPoint] - origPot[refPoint];
		for (int i = 0; i < nPts; i++)
			targ[i] -= origPot[i] + ref;
	}

	void LDAFunctional::calcPot(double* rho, double* targ) {
		switch (typ) {
		case LDAFunctionalType::X_SLATER: // slater exchange
		{
			double coef = -std::pow(3.0/PhysCon::pi, 1.0/3) * PhysCon::auE_ha * PhysCon::a0; //convert linear density to a.u., then to energy in SI
			for (int i = 0; i < nPts; i++)
				targ[i] = coef * std::pow(rho[i], 1.0/3);
			break;
		}
		case LDAFunctionalType::C_PW: // PW correlation
		{
			double aa = 0.031091, al = 0.21370, be1 = 7.5957, be2 = 3.5876, be3 = 1.6382, be4 = 0.49294;
			
			double crs, crho, q0, q1, q1p, drsdrho;
			double smallRho = 1e-10;
			for(int i = 0; i < nPts; i++){
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

	CompositePotential::CompositePotential(int nPts, int numSPots, int numDPots, int numWPots, Potential ** staticPots, Potential ** dynamicPots, Potential ** waveFuncDependentPots) :
		nPts(nPts), numSPots(numSPots), numDPots(numDPots), numWPots(numWPots), staticPots(staticPots), dynamicPots(dynamicPots), waveFuncDependentPots(waveFuncDependentPots)
	{
		v0 = (double*) sq_malloc(sizeof(double)*nPts);
		if (numSPots != 0)
			staticPots[0]->getV(0.0, v0);
		else
			std::fill_n(v0, nPts, 0.0);
		nv = (double*) sq_malloc(sizeof(double)*nPts);
		for (int i = 1; i < numSPots; i++) {
			staticPots[i]->getV(0.0, nv);
			vtls::addArrays(nPts, nv, v0);
		}
	}

	CompositePotential::~CompositePotential(){
		sq_free(staticPots);
		sq_free(dynamicPots);
		sq_free(waveFuncDependentPots);

		sq_free(v0);
		sq_free(nv);
	}

	void CompositePotential::getV(double t, double * targ) {
		vtls::copyArray(nPts, v0, targ);
		for (int i = 0; i < numDPots; i++) {
			dynamicPots[i]->getV(t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
	}

	void CompositePotential::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		vtls::copyArray(nPts, v0, targ);
		for (int i = 0; i < numDPots; i++) {
			dynamicPots[i]->getV(t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
		for (int i = 0; i < numWPots; i++) {
			waveFuncDependentPots[i]->getV(rho, psi, t, nv);
			vtls::addArrays(nPts, nv, targ);
		}
	}

	PotentialComplexity CompositePotential::getComplexity() {
		if (numWPots > 0)
			return PotentialComplexity::WAVEFUNCTION_DEPENDENT;
		else if (numDPots > 0)
			return PotentialComplexity::DYNAMIC;
		else
			return PotentialComplexity::STATIC;
	}

	PotentialManager::PotentialManager(int nPts) {
		PotentialManager::nPts = nPts;
	}

	void PotentialManager::addPotential(Potential * pot) {
		switch(pot->getComplexity()){
			case PotentialComplexity::STATIC:
				staticPots.push_back(pot);
				break;
			case PotentialComplexity::DYNAMIC:
				dynamicPots.push_back(pot);
				if(myComplex == PotentialComplexity::STATIC)
					myComplex = PotentialComplexity::DYNAMIC;
				break;
			case PotentialComplexity::WAVEFUNCTION_DEPENDENT:
				waveFuncDependentPots.push_back(pot);
				myComplex = PotentialComplexity::WAVEFUNCTION_DEPENDENT;
				break;
			default:
				throw std::runtime_error("Potential does not fall into an enum PotentialComplexity (something is very wrong?)");
		}
	}

	void PotentialManager::finishAddingPotentials() {
		int ns = staticPots.size();
		int nd = dynamicPots.size();
		int nw = waveFuncDependentPots.size();

		if (spots)
			delete[] spots;
		if (dpots)
			delete[] dpots;
		if (wpots)
			delete[] wpots;

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
		else
			throw std::runtime_error("Must run finishAddingPotentials() function on PotentialManager before using getV().");
	}

	void PotentialManager::getV(double* rho, std::complex<double> * psi, double t, double * targ) {
		if (pot)
			pot->getV(rho, psi, t, targ);
		else
			throw std::runtime_error("Must run finishAddingPotentials() function on PotentialManager before using getV().");
	}

	namespace ElectricFieldProfiles {
		ConstantFieldProfile::ConstantFieldProfile(int nPts, double * x, double eMax, double minX, double maxX) {
			fs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
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
			fs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
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

		CylindricalToLinearProfile::~CylindricalToLinearProfile() {
			sq_free(fs);
		}

		std::complex<double> * CylindricalToLinearProfile::getProfile() {
			return fs;
		}


		CylindricalToCutoffProfile::CylindricalToCutoffProfile(int nPts, double * x, double minX, double maxX, double r, double eMax, double enhFact, double decayLength) {
			fs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
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

		CylindricalToCutoffProfile::~CylindricalToCutoffProfile() {
			sq_free(fs);
		}

		std::complex<double> * CylindricalToCutoffProfile::getProfile() {
			return fs;
		}


		InMetalFieldProfile::InMetalFieldProfile(int nPts, double * x, double minX, double maxX, double eMax, double lam, std::complex<double> er, double cond) {
			fs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
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
			fs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
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
			sq_free(tre);
			sq_free(tim);
			sq_free(fx);
			sq_free(fre);
			sq_free(fim);
		}

		FileFieldProfile::~FileFieldProfile() {
			sq_free(fs);
		}

		std::complex<double> * FileFieldProfile::getProfile() {
			return fs;
		}


		ExponentialToLinearProfile::ExponentialToLinearProfile(int nPts, double* x, double minX, double maxX, double r, double eMax) {
			fs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts);
			double xn;
			for (int i = 0; i < nPts; i++) {
				xn = x[i];
				if (xn < minX || xn > maxX)
					fs[i] = 0.0;
				else if (xn < maxX - r)
					fs[i] = eMax * std::exp(-(xn - minX) / r);
				else
					fs[i] = eMax / r * std::exp(-(maxX - r) / r) * (maxX - xn);
			}
		}

		ExponentialToLinearProfile::~ExponentialToLinearProfile() {
			sq_free(fs);
		}

		std::complex<double>* ExponentialToLinearProfile::getProfile() {
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