#include "stdafx.h"

int main(int argc, char **argv)
{
	double dx = 1e-11;
	double dt = 1e-17;
	double xmax = 1e-8;
	int nPts = (int)(xmax / dx / 2.0) * 2;
	std::cout << nPts << std::endl;
	dx = xmax / (nPts - 1);

	double e0 = 20 * PhysCon::qe;
	double sig = 1e-9;

	double* v0 = new double[nPts];

	for (int i = 0; i < nPts; i++)
		v0[i] = -e0 * std::exp(-std::pow(i * dx - xmax / 2.0, 2) / (2 * sig * sig));

	std::complex<double>** psiPtr = new std::complex<double>*;
	int nelec = -1;

	KineticOperators::GenDisp_PSM* kin = new KineticOperators::GenDisp_PSM_FreeElec(nPts, dx, dt, 1.0);

	kin->findEigenStates(v0, -e0, -e0 / 2.0, psiPtr, &nelec);

	std::complex<double>* psi = *psiPtr;
	std::complex<double>* npsi = new std::complex<double>[nPts * nelec];

	double* spd = new double[nPts];
	std::fill_n(spd, nPts, 1.0);
	double* rho_ptr = new double[nPts];
	std::vector<double> rho(nPts);

	//move center of potential
	for (int i = 0; i < nPts; i++)
		v0[i] = -e0 * std::exp(-std::pow(i * dx - xmax / 3.0, 2) / (2 * sig * sig));

	namespace plt = matplotlibcpp;

	for (int i = 0; i < 1000; i++) {
		kin->stepOS_U2TU(psi, v0, spd, npsi, nelec);

		vtls::copyArray(nPts*nelec, npsi, psi);

		vtls::normSqr(nPts, &psi[nPts], rho_ptr);

		std::copy(rho_ptr, rho_ptr+nPts, rho.begin());

		plt::clf();
		plt::plot(rho);
		plt::title(std::to_string(i));
		plt::draw();
		plt::pause(0.05);
	}


	/*if (argc != 2) {
		std::cout << "Needs one input argument: config file address." << std::endl;
		if (argc > 2)
			std::cout << "If address includes spaces, it needs to be surrounded by quotes." << std::endl;
	}
	else {
		cfgParse::readCFG(argv[1]);
	}

	return 0;*/
}