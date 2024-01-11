#include "AbsorptiveRegions.h"
// OBSOLETE
namespace AbsorptiveRegions {
	AbsorptiveRegionVel::AbsorptiveRegionVel(double dx, double dt, int nPts, int left, int right, double rate) {
		AbsorptiveRegionVel::dx = dx;
		AbsorptiveRegionVel::dt = dt;
		AbsorptiveRegionVel::nPts = nPts;
		AbsorptiveRegionVel::left = left;
		AbsorptiveRegionVel::right = right;
		AbsorptiveRegionVel::rate = rate;
	}

	void AbsorptiveRegionVel::decay(std::complex<double> * psi) {
		double v1, v2;
		if (psi[left] != 0.0)
			v2 = std::abs(PhysCon::hbar*std::imag(std::conj(psi[left])*vtls::firstDerivative(nPts, psi, left, dx)) / (PhysCon::me*psi[left] * psi[left]));
		else
			v2 = 0.0;
		for (int i = left + 1; i < right; i++) {
			v1 = v2;
			if (psi[i] != 0.0)
				v2 = std::abs(PhysCon::hbar*std::imag(std::conj(psi[i])*vtls::firstDerivative(nPts, psi, i, dx)) / (PhysCon::me*psi[i] * psi[i]));
			else
				v2 = 0.0;
			psi[i - 1] /= (1.0 + std::sqrt(v1*rate*dt*dx));
		}
		psi[right - 1] /= (1.0 + std::sqrt(v2*rate*dt*dx));
	}

	AbsorptiveRegionVelSmooth::AbsorptiveRegionVelSmooth(double dx, double dt, int nPts, int inner, int outer, double rate) {
		AbsorptiveRegionVelSmooth::dx = dx;
		AbsorptiveRegionVelSmooth::dt = dt;
		AbsorptiveRegionVelSmooth::nPts = nPts;
		AbsorptiveRegionVelSmooth::inner = inner;
		AbsorptiveRegionVelSmooth::outer = outer;
		AbsorptiveRegionVelSmooth::rate = rate;
		size = std::abs(outer - inner);
		mask = new double[size];
		if (outer > inner) {
			for (int i = 0; i < size; i++) {
				double k = (double)(i + 1) / size;
				mask[i] = (3.0*k *k - 2.0*k*k*k)*rate;
			}
			left = inner;
			right = outer;
		}
		else {
			for (int i = 0; i < size; i++) {
				double k = (double)(size - i) / size;
				mask[i] = (3.0*k*k - 2.0*k*k*k)*rate;
			}
			left = outer;
			right = inner;
		}
	}

	void AbsorptiveRegionVelSmooth::decay(std::complex<double> * psi) {
		double v1, v2;
		if (psi[left] != 0.0)
			v2 = std::abs(PhysCon::hbar*std::imag(std::conj(psi[left])*vtls::firstDerivative(nPts, psi, left, dx)) / (PhysCon::me*psi[left] * psi[left]));
		else
			v2 = 0.0;
		for (int i = left + 1; i < right; i++) {
			v1 = v2;
			if (psi[i] != 0.0)
				v2 = std::abs(PhysCon::hbar*std::imag(std::conj(psi[i])*vtls::firstDerivative(nPts, psi, i, dx)) / (PhysCon::me*psi[i] * psi[i]));
			else
				v2 = 0.0;
			psi[i - 1] /= std::sqrt((1.0 + v1 * mask[i - left - 1] * dt));
		}
		psi[right - 1] /= std::sqrt((1.0 + v2 * mask[size - 1] * dt));
	}

	AbsorptiveRegionManager::AbsorptiveRegionManager() {
		numReg = 0;
	}

	void AbsorptiveRegionManager::addAbsorptiveRegion(AbsorptiveRegion* reg) {
		regs.push_back(reg);
		numReg++;
	}

	void AbsorptiveRegionManager::decay(std::complex<double> * psi) {
		for (int i = 0; i < numReg; i++)
			regs[i]->decay(psi);
	}

	std::complex<double>* getSmoothedTimePhaseDecay(int len, int inner, int outer, double rate) {
		std::complex<double> * mask = new std::complex<double>[len];
		std::fill_n(mask, len, 0.0);
		int size = std::abs(outer - inner);
		double k;
		if (outer > inner) {
			for (int i = 0; i < size; i++) {
				k = (double)(i + 1) / size;
				mask[i + inner] = (924.0*std::pow(k, 13) -
					6006.0*std::pow(k, 12) +
					16380.0*std::pow(k, 11) -
					24024.0*std::pow(k, 10) +
					20020.0*std::pow(k, 9) -
					9009.0*std::pow(k, 8) +
					1716.0*std::pow(k, 7))*(rate*PhysCon::im);
			}
		}
		else {
			for (int i = 0; i < size; i++) {
				k = (double)(size - i) / size;
				mask[i + outer] = (924.0*std::pow(k, 13) -
					6006.0*std::pow(k, 12) +
					16380.0*std::pow(k, 11) -
					24024.0*std::pow(k, 10) +
					20020.0*std::pow(k, 9) -
					9009.0*std::pow(k, 8) +
					1716.0*std::pow(k, 7))*(rate*PhysCon::im);
			}
		}
		return mask;
	}

	double* getSmoothedSpatialDampDecay(int len, int inner, int outer, double rate) {
		double* mask = new double[len];
		std::fill_n(mask, len, 1.0);
		int size = std::abs(outer - inner);
		double k;
		if (outer < inner) {
			for (int i = 0; i < size; i++) {
				k = (double)(i + 1) / size;
				mask[i + outer] = std::pow((924.0 * std::pow(k, 13) -
					6006.0 * std::pow(k, 12) +
					16380.0 * std::pow(k, 11) -
					24024.0 * std::pow(k, 10) +
					20020.0 * std::pow(k, 9) -
					9009.0 * std::pow(k, 8) +
					1716.0 * std::pow(k, 7)) , (rate));
			}
		}
		else {
			for (int i = 0; i < size; i++) {
				k = (double)(size - i) / size;
				mask[i + inner] = std::pow((924.0 * std::pow(k, 13) -
					6006.0 * std::pow(k, 12) +
					16380.0 * std::pow(k, 11) -
					24024.0 * std::pow(k, 10) +
					20020.0 * std::pow(k, 9) -
					9009.0 * std::pow(k, 8) +
					1716.0 * std::pow(k, 7)) , (rate));
			}
		}
		return mask;
	}
}