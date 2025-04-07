#include "WfcRhoTools.h"
#include "PhysCon.h"
#include <stdexcept>

namespace WfcToRho {
	void BoundFermiGas::calcWeights(int nElec, const double* energies, double* weights, NormalizationScheme norm) {
		if(norm == NormalizationScheme::UNNORMALIZED)
			throw std::runtime_error("WfcToRho::BoundFermiGas Cannot use UNNORMALIZED scheme with BoundFermiGas");

		double minE = vtls::min(nElec, energies);
		double maxE = vtls::max(nElec, energies);
		//Find center of Fermi slab, and set new 0-energy accoridngly, and convert to Fermi energy difference
		double bottom = (maxE + minE) / 2 - ef / 2;
		for (int i = 0; i < nElec; i++)
			weights[i] = ef - (energies[i] - bottom);
		//Get factored prefactor
		double fact = 0;
		for (int i = 0; i < nElec; i++)
			fact += weights[i];
		fact = 2.0 / (3.0 * PhysCon::pi) * PhysCon::me * ef / (PhysCon::hbar * PhysCon::hbar) * nElec / fact;
		//Combine factored prefactor with energy difference
		for (int i = 0; i < nElec; i++)
			weights[i] *= fact;
	}

	void SemiInfiniteFermiGas::calcWeights(int nElec, const double* energies, double* weights, NormalizationScheme norm) {
		if(norm == NormalizationScheme::NORMALIZED)
			throw std::runtime_error("WfcToRho::SemiInfiniteFermiGas expects a wavefunction which has a normalization dicated by boundary conditions, not the total norm.");

		double minE = vtls::min(nElec, energies);
		double maxE = vtls::max(nElec, energies);
		double bottom = maxE - ef;

		// get energy range boundaries relative to bottom
		double* energyRangeBoundaries = (double*) sq_malloc(sizeof(double) * (nElec+1));
		energyRangeBoundaries[0] = 0.0;
		for (int i = 0; i < nElec-1; i++){
			if(energies[i] > energies[i+1])
				throw std::runtime_error("WfcToRho::SemiInfiniteFermiGas::calcWeights: energies not sorted");
			energyRangeBoundaries[i + 1] = (energies[i] + energies[i+1])/2.0 - bottom;
		}
		energyRangeBoundaries[nElec] = ef;

		// Fermi gas (one direction only): N = sqrt(2)/(3 pi^2) (Ef m/hbar^2)^(3/2)
		// density wrt k \propto Ef-E
		// density between energies \propto (Ef-E1)^2-(Ef-E0)^2
		// each state has density = N * [(Ef-E1)^2-(Ef-E0)^2] / Ef^2
		double n0 = std::sqrt(2)/(3.0*PhysCon::pi*PhysCon::pi)*std::pow(ef * PhysCon::me/PhysCon::hbar/PhysCon::hbar, 1.5);
		for(int i = 0; i < nElec; i++)
			weights[i] = n0 * (std::pow(ef - energyRangeBoundaries[i], 2) - std::pow(ef - energyRangeBoundaries[i+1], 2)) / ef / ef;

		sq_free(energyRangeBoundaries);
	}

	//NOTE: fl here is the Fermi level (relative to vacuum) of the model system, not Fermi energy. Typically -W
	//		ef is the Fermi energy
	//		leff is the effective well size
	FromDOS::FromDOS(double fl, double ef, double leff, const char* fil) : ef(ef), leff(leff) {
		std::fstream ifil = std::fstream(fil, std::ios::in | std::ios::binary);

		int nes;
		ifil.read(reinterpret_cast<char*>(&nes), sizeof(int));
		double *fE = (double*) sq_malloc(sizeof(double)*nes);
		double *dos = (double*) sq_malloc(sizeof(double)*nes);
		double *dosI = (double*) sq_malloc(sizeof(double)*nes);
		ifil.read(reinterpret_cast<char*>(fE), sizeof(double) * nes);
		ifil.read(reinterpret_cast<char*>(dos), sizeof(double) * nes);
		ifil.close();
		double dE = (fE[nes - 1] - fE[0]) / (nes - 1.0);
		vtlsInt::cumIntTrapz(nes, dos, dE, dosI);

		dosISpline = boost::math::interpolators::cardinal_cubic_b_spline<double>(dosI, nes, fE[0]*PhysCon::qe - fl, dE*PhysCon::qe);
		
		sq_free(dos);
		sq_free(dosI);
		sq_free(fE);
	}

	void FromDOS::calcWeights(int nElec, const double* energies, double* weights, NormalizationScheme norm) {
		//approximate effective width of well... will need to be reconsidered if using non square-ish wells
		//corrects for lost normalized density for larger wells (densities should be O(1))
		//double leff = PhysCon::hbar * 2.0 * PhysCon::pi * nElec / (2.0 * std::sqrt(2.0 * ef * PhysCon::me));

		if(norm == NormalizationScheme::UNNORMALIZED)
			std::cerr << "Warning: FromDOS: UNNORMALIZED scheme is not recommended with FromDOS, but could work if you know what you're doing." << std::endl;

		if (nElec == 1) {
			std::cout << "FromDOS: Only single electron provided, weight set to total DOS from energies[0]-ef to energies[0]" << std::endl;
			weights[0] = leff * (dosISpline(energies[0]) - dosISpline(energies[0]-ef));
			return;
		}

		//sort energy
		int* idx = (int*) sq_malloc(sizeof(int)*nElec);
		double* sortedEnergies = (double*) sq_malloc(sizeof(double) * nElec);
		vtls::copyArray( nElec, energies, sortedEnergies);

		for (int i = 0; i < nElec; i++)
			idx[i] = i;
		vtls::insertSort_idxs(nElec, sortedEnergies, idx);

		//represented energies are half-way between adjacent energies
		//eg, if we have energy states E = 0, 1, 3, then the state of energy 1 represents energies 2 <- 0.5
		for (int i = 1; i < nElec - 1; i++) {
			weights[idx[i]] = leff * (dosISpline((sortedEnergies[i] + sortedEnergies[i + 1]) / 2.0) - dosISpline((sortedEnergies[i] + sortedEnergies[i - 1]) / 2.0));
		}

		//the bottom state represents halfway above and the same amount below
		//eg, if we have energy states E = 0, 1 then state with energy 0 represents energies 0.5 <- -0.5
		weights[idx[0]] = leff * (dosISpline((sortedEnergies[0] + sortedEnergies[1]) / 2.0) - dosISpline((3.0 * sortedEnergies[0] - sortedEnergies[1]) / 2.0));

		//the top state is similar to the bottom state
		weights[idx[nElec - 1]] = leff * (dosISpline((3.0 * sortedEnergies[nElec - 1] - sortedEnergies[nElec - 2]) / 2.0) - dosISpline((sortedEnergies[nElec - 1] + sortedEnergies[nElec - 2]) / 2.0));

		sq_free(idx);
		sq_free(sortedEnergies);
	}

	void Density::calcRawRho(int nPts, int nElec, const double* weights, const std::complex<double>* psi, double* psi2_work, double* rho){
		std::fill_n(rho, nPts, 0);
		vtls::normSqr(nPts * nElec, psi, psi2_work);
		for (int i = 0; i < nElec; i++) {
			double pref = weights[i];
			for (int j = 0; j < nPts; j++)
				rho[j] += psi2_work[i * nPts + j] * pref;
		}
	}

	void DirectDensity::calcRho(int nPts, int nElec, double dx, double* rho) {
		return;
	}

	CylindricalDensity::CylindricalDensity(double center, double radius, double minX) : center(center), radius(radius), minX(minX) {};

	CylindricalDensity::~CylindricalDensity(){
		if (thinning)
			sq_free(thinning);
	}

	void CylindricalDensity::doFirst(int nPts, double dx) {
		first = false;
		if(thinning)
			sq_free(thinning);
		thinning = (double*) sq_malloc(sizeof(double)*nPts);

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

	void CylindricalDensity::calcRho(int nPts, int nElec, double dx, double* rho) {
		if (first)
			doFirst(nPts, dx);

		vtls::seqMulArrays(endIndex-startIndex, &thinning[startIndex], &rho[startIndex]);
	}
	
	GaussianSmoothedDensity::~GaussianSmoothedDensity(){
		if(tempRho)
			sq_free(tempRho);
		if(conv)
			delete conv;
	}

	void GaussianSmoothedDensity::calcRho(int nPts, int nElec, double dx, double* rho) {
		if (first) {
			//Initialize variables
			if(tempRho)
				sq_free(tempRho);
			double* mask = (double*)sq_malloc(sizeof(double)*nPts);
			tempRho = (double*)sq_malloc(sizeof(double)*nPts);
			first = false;

			//Initialize Gaussian mask (in k space)
			for (int i = 0; i < nPts / 2; i++) {
				mask[i] = 1.0 / (sig/dx * std::sqrt(2.0 * PhysCon::pi)) * std::exp(-0.5 / (sig * sig) * (i * i * dx * dx));
				mask[nPts - i - 1] = 1.0 / (sig/dx * std::sqrt(2.0 * PhysCon::pi)) * std::exp(-0.5 / (sig * sig) * ((i+1) * (i+1) * dx * dx));
			}
			if(nPts%2)
				mask[nPts/2] = 1.0 / (sig/dx * std::sqrt(2.0 * PhysCon::pi)) * std::exp(-0.5 / (sig * sig) * (nPts * nPts / 4.0 * dx * dx));
			
			vtls::scaMulArray(nPts, 1.0 / vtlsInt::rSum(nPts, mask, dx), mask); //normalize

			//Initialize FFT for convolution
			if(conv)
				delete conv;
			conv = new vtls::MaskConvolver<double>(nPts, mask);

			sq_free(mask);
		}

		conv->compute(rho);
	}
}