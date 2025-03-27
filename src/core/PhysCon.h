/**
 * @file PhysCon.h
 * @brief Contains various physical constants.
 */
#pragma once

#include <complex>
/**
 * @namespace PhysCon
 * @brief Contains various physical constants.
 * One may edit the source with USING_AU=true to use atomic units.
 */
namespace PhysCon {
	static const bool USING_AU = false;
	using namespace std;
	static const double hbar = USING_AU ? 1.0 : 1.054571800e-34;	/// Planck's Reduced Constant
	static const double me = USING_AU ? 1.0 : 9.10938356e-31;		/// Mass of Electron
	static const double qe = USING_AU ? 1.0 : 1.60217662e-19;		/// Charge of Electron
	static const double eV = USING_AU ? 1.0/27.211386245 : 1.60217662e-19;	/// Electron Volt
	static const double c = USING_AU ? 137.035999 : 299792458.0;		/// Speed of Light
	static const double auE_ha = USING_AU ? 1.0 : 4.35974417e-18;		/// Hartree Atomic Unit of Energy
	static const double auE_ry = auE_ha / 2.0;		/// Rydberg Atomic Unit of Energy
	static const std::complex<double> im = 0.0 + 1.0i;			/// Imaginary Number
	static const double pi = M_PI;		/// Pi
	static const double e0 = USING_AU ? 0.25/pi : 8.854187817e-12;	/// Vacuum Permittivity
	static const double mu0 = 1.0/(c*c*e0); /// Vacuum Magnetic Permeability
	static const double a0 = USING_AU ? 1.0 : 5.2917721092e-11;	/// Bohr Radius
	static const double k = 1.0/(4.0*pi*e0); /// Coulomb electric constant
};