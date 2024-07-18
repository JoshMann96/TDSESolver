#include "WfcRhoTools.h"
#include "MathTools.h"
#include "fftw3.h"
#include <stdexcept>

namespace WfcToRho {
	void calcEnergies(int nElec, int nPts, double dx, std::complex<double>* psi, double* totPot, KineticOperators::KineticOperator* kin, double* energies) {
		if (nElec < 1) {
			std::cout << "calcEnergies: Number of electrons is not finite! Failed to initialize." << std::endl;
			throw -1;
		}
		double* rho = new double[nPts];
		for (int i = 0; i < nElec; i++) {
			vtls::normSqr(nPts, &psi[i * nPts], rho);
			energies[i] = vtlsInt::rSumMul(nPts, rho, totPot, dx) + kin->evaluateKineticEnergy(&psi[i*nPts]);
			//potential energy + kinetic energy
		}
		delete[] rho;
	}


	void FermiGasDistro::calcWeights(int nElec, double* energies, double* weights) {
		double minE = energies[0], maxE = energies[0];
		//Find min/max energy
		for (int i = 0; i < nElec; i++) {
			if (energies[i] > maxE)
				maxE = energies[i];
			if (energies[i] < minE)
				minE = energies[i];
		}
		//Find center of Fermi slab, and set new 0-energy accoridngly, and convert to Fermi energy difference
		double bottom = (maxE + minE) / 2 - ef / 2;
		for (int i = 0; i < nElec; i++)
			weights[i] = ef - (energies[i] - bottom);
		//Get factored prefactor
		double fact = 0;
		for (int i = 0; i < nElec; i++)
			fact += weights[i];
		fact = 2 / (3 * PhysCon::pi) * PhysCon::me * ef / (PhysCon::hbar * PhysCon::hbar) * nElec / fact;
		//Combine factored prefactor with energy difference
		for (int i = 0; i < nElec; i++)
			weights[i] *= fact;
	}


	//NOTE: fl here is the Fermi level (relative to vacuum) of the model system, not Fermi energy. Typically -W
	//		ef is the Fermi energy
	//		leff is the effective well size
	FromDOS::FromDOS(double fl, double ef, double leff, const char* fil) : ef(ef), leff(leff) {
		std::fstream ifil = std::fstream(fil, std::ios::in | std::ios::binary);

		int nes;
		ifil.read(reinterpret_cast<char*>(&nes), sizeof(int));
		double *fE = new double[nes];
		double *dos = new double[nes];
		double *dosI = new double[nes];
		ifil.read(reinterpret_cast<char*>(fE), sizeof(double) * nes);
		ifil.read(reinterpret_cast<char*>(dos), sizeof(double) * nes);
		ifil.close();
		double dE = (fE[nes - 1] - fE[0]) / (nes - 1.0);
		vtlsInt::cumIntTrapz(nes, dos, dE, dosI);

		dosISpline = boost::math::interpolators::cardinal_cubic_b_spline<double>(dosI, nes, fE[0]*PhysCon::qe - fl, dE*PhysCon::qe);
		
		delete[] dos;
		delete[] dosI;
		delete[] fE;
	}

	void FromDOS::calcWeights(int nElec, double* energies, double* weights) {
		//approximate effective width of well... will need to be reconsidered if using non square-ish wells
		//corrects for lost normalized density for larger wells (densities should be O(1))
		//double leff = PhysCon::hbar * 2.0 * PhysCon::pi * nElec / (2.0 * std::sqrt(2.0 * ef * PhysCon::me));

		if (nElec == 1) {
			std::cout << "FromDOS: Only single electron provided, weight set to total DOS from energies[0]-ef to energies[0]" << std::endl;
			weights[0] = leff * (dosISpline(energies[0]) - dosISpline(energies[0]-ef));
			return;
		}

		//sort energy
		int* idx = new int[nElec];
		for (int i = 0; i < nElec; i++)
			idx[i] = i;
		vtls::insertSort_idxs(nElec, energies, idx);

		//represented energies are half-way between adjacent energies
		//eg, if we have energy states E = 0, 1, 3, then the state of energy 1 represents energies 2 <- 0.5
		for (int i = 1; i < nElec - 1; i++) {
			weights[idx[i]] = leff * (dosISpline((energies[i] + energies[i + 1]) / 2.0) - dosISpline((energies[i] + energies[i - 1]) / 2.0));
		}

		//the bottom state represents halfway above and the same amount below
		//eg, if we have energy states E = 0, 1 then state with energy 0 represents energies 0.5 <- -0.5
		weights[idx[0]] = leff * (dosISpline((energies[0] + energies[1]) / 2.0) - dosISpline((3.0 * energies[0] - energies[1]) / 2.0));

		//the top state is similar to the bottom state
		weights[idx[nElec - 1]] = leff * (dosISpline((3.0 * energies[nElec - 1] - energies[nElec - 2]) / 2.0) - dosISpline((energies[nElec - 1] + energies[nElec - 2]) / 2.0));

		delete[] idx;
	}

	DirectDensity::~DirectDensity(){
		if(psi2)
			delete[] psi2;
	}

	void DirectDensity::calcRho(int nPts, int nElec, double dx, double* weights, std::complex<double>* psi, double* rho) {
		if (first) {
			if(psi2)
				delete[] psi2;
			psi2 = new double[nPts * nElec];
			first = 0;
		}

		std::fill_n(rho, nPts, 0);
		vtls::normSqr(nPts * nElec, psi, psi2);
		for (int i = 0; i < nElec; i++) {
			double pref = weights[i];
			for (int j = 0; j < nPts; j++)
				rho[j] += psi2[i * nPts + j] * pref;
		}
	}

	CylindricalDensity::CylindricalDensity(double center, double radius, double minX) : center(center), radius(radius), minX(minX) {};

	CylindricalDensity::~CylindricalDensity(){
		if (thinning)
			delete[] thinning;
	}

	void CylindricalDensity::doFirst(int nPts, double dx) {
		first = 0;
		if(baseDens == nullptr)
			throw std::runtime_error("CylindricalDensity: Base density calculator must be set with setBaseDens before use.");
		if(thinning)
			delete[] thinning;
		thinning = new double[nPts];

		std::fill_n(thinning, nPts, 1.0);

		if(radius > dx/2){
			startIndex = std::max(0, (int)std::floor((center + radius - minX) / dx));
			endIndex = nPts;
		}
		else if(radius < -dx/2){
			startIndex = 0;
			endIndex = std::min(nPts, (int)std::ceil((center + radius - minX) / dx));
		}
		else
			throw std::runtime_error("CylindricalDensity: Radius must be larger than half of grid size.");
		
		for(int i = startIndex; i < endIndex; i++)
			thinning[i] = radius / (i * dx - center - minX);
	}

	void CylindricalDensity::calcRho(int nPts, int nElec, double dx, double* weights, std::complex<double>* psi, double* rho) {
		if (first)
			doFirst(nPts, dx);

		baseDens->calcRho(nPts, nElec, dx, weights, psi, rho);

		vtls::seqMulArrays(endIndex-startIndex, &thinning[startIndex], &rho[startIndex]);
	}
	
	GaussianSmoothedDensity::~GaussianSmoothedDensity(){
		if(psi2)
			fftw_free(psi2);
		if(tempRho)
			fftw_free(tempRho);
		delete conv;
	}

	void GaussianSmoothedDensity::calcRho(int nPts, int nElec, double dx, double* weights, std::complex<double>* psi, double* rho) {
		if (first) {
			//Initialize variables
			if(psi2)
				fftw_free(psi2);
			if(tempRho)
				fftw_free(tempRho);
			psi2 = (double*)fftw_malloc(sizeof(double)*nPts*nElec);
			double* mask = (double*)fftw_malloc(sizeof(double)*nPts);
			tempRho = (double*)fftw_malloc(sizeof(double)*nPts);
			first = 0;

			//Initialize Gaussian mask (in k space)
			for (int i = 0; i < nPts / 2; i++) {
				mask[i] = 1.0 / (sig/dx * std::sqrt(2.0 * PhysCon::pi)) * std::exp(-0.5 / (sig * sig) * (i * i * dx * dx));
				mask[nPts - i - 1] = 1.0 / (sig/dx * std::sqrt(2.0 * PhysCon::pi)) * std::exp(-0.5 / (sig * sig) * ((i+1) * (i+1) * dx * dx));
			}
			if(nPts%2)
				mask[nPts/2] = 1.0 / (sig/dx * std::sqrt(2.0 * PhysCon::pi)) * std::exp(-0.5 / (sig * sig) * (nPts * nPts / 4.0 * dx * dx));
				
			//Initialize FFT for convolution
			conv = new vtls::MaskConvolver<double>(nPts, mask);

			delete[] mask;
		}

		std::fill_n(tempRho, nPts, 0);
		vtls::normSqr(nPts * nElec, psi, psi2);
		for (int i = 0; i < nElec; i++) {
			double pref = weights[i];
			for (int j = 0; j < nPts; j++)
				tempRho[j] += psi2[i * nPts + j] * pref;
		}

		conv->compute(tempRho, rho);

		for (int i = 0; i < nPts; i++)
			rho[i] = std::real(tempRho[i]);
	}
}