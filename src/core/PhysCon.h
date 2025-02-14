#pragma once

#include <complex>
// Various physical constants.
namespace PhysCon {
	static const bool USING_AU = false;
	using namespace std;
	static const double hbar = USING_AU ? 1.0 : 1.054571800e-34;	//Planck's Reduced Constant (SI)
	static const double me = USING_AU ? 1.0 : 9.10938356e-31;		//Mass of Electron (SI)
	static const double qe = USING_AU ? 1.0 : 1.60217662e-19;		//Charge of Electron (SI)
	static const double c = USING_AU ? 137.035999 : 299792458.0;		//Speed of Light (SI)
	static const double auE_ha = USING_AU ? 1.0 : 4.35974417e-18;		//Hartree Atomic Unit of Energy (SI)
	static const double auE_ry = auE_ha / 2.0;		//Rydberg Atomic Unit of Energy (SI)
	static const std::complex<double> im = 0.0 + 1.0i;			//Imaginary Number (NU)
	static const double pi = M_PI;		//Pi
	static const double e0 = USING_AU ? 0.25/pi : 8.854187817e-12;	//Vacuum Permittivity (SI)
	static const double mu0 = 1.0/(c*c*e0); //Vacuum Magnetic Permeability (SI)
	static const double a0 = USING_AU ? 1.0 : 5.2917721092e-11;	//Bohr Radius (SI)
	static const double k = 1.0/(4.0*pi*e0); //Coulomb electric constant (SI)
};