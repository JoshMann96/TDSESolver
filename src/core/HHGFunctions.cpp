#include "HHGFunctions.h"
#include "PhysCon.h"

namespace HHGFunctions {
	double getZim(double ef, double w) {
		double nEf = ef / PhysCon::auE_ha;
		double nBulk = 1.0 / (3.0*std::pow(PhysCon::pi, 2))*std::pow(2.0*nEf, 2);
		double rs = std::cbrt(3.0 / (4.0*PhysCon::pi*nBulk));
		return (-0.2*rs + 1.25)*PhysCon::a0;
	}

	double getIdealDX(double maxE, double error) {
		return PhysCon::hbar / std::sqrt(2.0*PhysCon::me*maxE)*std::acos(1.0 - error);
	}

	double getIdealDT(double dx) {
		return dx * dx*PhysCon::me / PhysCon::hbar*0.5;
	}

	double getPonderomotiveEnergy(double eMax, double lam) {
		double nc = eMax * PhysCon::qe / (PhysCon::pi * 2 * PhysCon::c / lam);
		return 0.25*nc*nc / PhysCon::me;
	}

	double getKeldyshParameter(double w, double eMax, double lam) {
		return std::sqrt(w / (2.0*getPonderomotiveEnergy(eMax, lam)));
	}
}