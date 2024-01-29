#include "WfcRhoTools.h"

namespace WfcToRho {
	void calcEnergies(int nelec, int nPts, double dx, std::complex<double>* psi, double* totPot, KineticOperators::KineticOperator* kin, double* energies) {
		if (nelec < 1) {
			std::cout << "calcEnergies: Number of electrons is not finite! Failed to initialize." << std::endl;
			throw -1;
		}
		double* rho = new double[nPts];
		for (int i = 0; i < nelec; i++) {
			vtls::normSqr(nPts, &psi[i * nPts], rho);
			energies[i] = vtlsInt::rSumMul(nPts, rho, totPot, dx) + kin->evaluateKineticEnergy(&psi[i*nPts]);
			//potential energy + kinetic energy
		}
		delete[] rho;

		/*std::complex<double>* sc1 = new std::complex<double>[nPts];
		std::complex<double>* sc2 = new std::complex<double>[nPts];
		for (int i = 0; i < nelec; i++) {
			vtls::secondDerivative(nPts, &psi[i * nPts], sc1, dx);
			vtls::scaMulArray(nPts, -PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me), sc1);
			vtls::seqMulArrays(nPts, totPot, &psi[i * nPts], sc2);
			vtls::addArrays(nPts, sc2, sc1);
			for (int j = 0; j < nPts; j++)
				sc2[j] = std::conj(psi[i * nPts + j]);
			energies[i] = std::real(vtlsInt::simpsMul(nPts, sc1, sc2, dx));
		}
		delete[] sc1, sc2;*/
	}


	void FermiGasDistro::calcWeights(int nelec, double* energies, double* weights) {
		double minE = energies[0], maxE = energies[0];
		//Find min/max energy
		for (int i = 0; i < nelec; i++) {
			if (energies[i] > maxE)
				maxE = energies[i];
			if (energies[i] < minE)
				minE = energies[i];
		}
		//Find center of Fermi slab, and set new 0-energy accoridngly, and convert to Fermi energy difference
		double bottom = (maxE + minE) / 2 - ef / 2;
		for (int i = 0; i < nelec; i++)
			weights[i] = ef - (energies[i] - bottom);
		//Get factored prefactor
		double fact = 0;
		for (int i = 0; i < nelec; i++)
			fact += weights[i];
		fact = 2 / (3 * PhysCon::pi) * PhysCon::me * ef / (PhysCon::hbar * PhysCon::hbar) * nelec / fact;
		//Combine factored prefactor with energy difference
		for (int i = 0; i < nelec; i++)
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

	void FromDOS::calcWeights(int nelec, double* energies, double* weights) {
		//approximate effective width of well... will need to be reconsidered if using non square-ish wells
		//corrects for lost normalized density for larger wells (densities should be O(1))
		//double leff = PhysCon::hbar * 2.0 * PhysCon::pi * nelec / (2.0 * std::sqrt(2.0 * ef * PhysCon::me));

		if (nelec == 1) {
			std::cout << "FromDOS: Only single electron provided, weight set to total DOS from energies[0]-ef to energies[0]" << std::endl;
			weights[0] = leff * (dosISpline(energies[0]) - dosISpline(energies[0]-ef));
			return;
		}

		//sort energy
		long* idx = new long[nelec];
		for (int i = 0; i < nelec; i++)
			idx[i] = i;
		vtls::insertSort_idxs(nelec, energies, idx);

		//represented energies are half-way between adjacent energies
		//eg, if we have energy states E = 0, 1, 3, then the state of energy 1 represents energies 2 <- 0.5
		for (int i = 1; i < nelec - 1; i++) {
			weights[idx[i]] = leff * (dosISpline((energies[i] + energies[i + 1]) / 2.0) - dosISpline((energies[i] + energies[i - 1]) / 2.0));
		}

		//the bottom state represents halfway above and the same amount below
		//eg, if we have energy states E = 0, 1 then state with energy 0 represents energies 0.5 <- -0.5
		weights[idx[0]] = leff * (dosISpline((energies[0] + energies[1]) / 2.0) - dosISpline((3.0 * energies[0] - energies[1]) / 2.0));

		//the top state is similar to the bottom state
		weights[idx[nelec - 1]] = leff * (dosISpline((3.0 * energies[nelec - 1] - energies[nelec - 2]) / 2.0) - dosISpline((energies[nelec - 1] + energies[nelec - 2]) / 2.0));

		delete[] idx;
	}

	void DirectDensity::calcRho(int nPts, int nelec, double dx, double* weights, std::complex<double>* psi, double* rho) {
		if (first) {
			psi2 = new double[nPts * nelec];
			first = 0;
		}

		std::fill_n(rho, nPts, 0);
		vtls::normSqr(nPts * nelec, psi, psi2);
		for (int i = 0; i < nelec; i++) {
			double pref = weights[i];
			for (int j = 0; j < nPts; j++)
				rho[j] += psi2[i * nPts + j] * pref;
		}
	}


	void GaussianSmoothedDensity::calcRho(int nPts, int nelec, double dx, double* weights, std::complex<double>* psi, double* rho) {
		if (first) {
			//Initialize variables
			psi2 = (double*)fftw_malloc(sizeof(double)*nPts*nelec);
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
		}

		std::fill_n(tempRho, nPts, 0);
		vtls::normSqr(nPts * nelec, psi, psi2);
		for (int i = 0; i < nelec; i++) {
			double pref = weights[i];
			for (int j = 0; j < nPts; j++)
				tempRho[j] += psi2[i * nPts + j] * pref;
		}

		conv->compute(tempRho, rho);

		for (int i = 0; i < nPts; i++)
			rho[i] = std::real(tempRho[i]);
	}
}