#pragma once
// Various physical constants.
namespace PhysCon {
	using namespace std;
	static const double hbar = 1.054571800e-34;	//Planck's Reduced Constant (SI)
	static const double me = 9.10938356e-31;		//Mass of Electron (SI)
	static const double qe = 1.60217662e-19;		//Charge of Electron (SI)
	static const double c = 299792458.0;		//Speed of Light (SI)
	static const double e0 = 8.854187817e-12;	//Vacuum Permittivity (SI)
	static const double mu0 = 4.0*3.141592653589e-7; //Vacuum Magnetic Permeability (SI)
	static const double auE_ha = 4.35974417e-18;		//Hartree Atomic Unit of Energy (SI)
	static const double auE_ry = 2.17987236e-18;		//Rydberg Atomic Unit of Energy (SI)
	static const std::complex<double> im = 0.0 + 1.0i;			//Imaginary Number (NU)
	static const double pi = M_PI;		//Pi
	static const double a0 = 5.2917721092e-11;	//Bohr Radius (SI)
	static const double k = 8987551787.3681764; //Coulomb electric constant (SI)
};