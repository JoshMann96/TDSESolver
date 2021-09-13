#include "stdafx.h"

int main(int argc, char **argv)
{
	/*double dx = 1e-12;
	double dt = 1e-18;
	double xmax = 1e-8;
	int nPts = (int)(xmax / dx / 2.0) * 2;
	dx = xmax / (nPts - 1);

	double e0 = 10 * PhysCon::qe;
	double k0 = std::sqrt(2.0 * PhysCon::me * e0 / (PhysCon::hbar * PhysCon::hbar));
	double sig = 1e-10;*/


	if (argc != 2) {
		std::cout << "Needs one input argument: config file address." << std::endl;
		if (argc > 2)
			std::cout << "If address includes spaces, it needs to be surrounded by quotes." << std::endl;
	}
	else {
		cfgParse::readCFG(argv[1]);
	}

	return 0;
}