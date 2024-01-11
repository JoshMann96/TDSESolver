#pragma once
#include "CommonHeader.h"

// Functions that are useful for HHG simulations. Everything is in SI units.
namespace HHGFunctions {
	// Calculates image plane position relative to center of Jellium potential, given Fermi energy (ef) and work function (w).
	double getZim(double ef, double w);

	// Gets the ideal dx spacing for the simulation given the maximum energy and the expected error at that energy (0-2].
	double getIdealDX(double maxE, double error);

	// Gets ideal dt spacing given dx.
	double getIdealDT(double dx);

	// Calculates ponderomotive energy of light field given maximum electric field strength (eMax) and wavelength (lam).
	double getPonderomotiveEnergy(double eMax, double lam);

	// Calculates the keldysh parameter given work function (w), maximum electric field strength (eMax), and wavelength (lam).
	double getKeldyshParameter(double w, double eMax, double lam);
}