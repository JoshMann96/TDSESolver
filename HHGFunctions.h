#pragma once
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

	// Sets up HHG simulation with the basics.
	// PARAMETERS:
	// emax = Max electric field | lam = Wavelength | tau = Gaussian envelope time constant | peakEnvT = Time in simulation where peak of envelope is
	// pulseSmoothTime = Time at start of simulation where field slowly transitions from nothing to gaussian pulse
	// ef = Fermi energy | w = Work function | jellCenter = Center for jellium potential (not imaginary plane)
	// atomCenter = Center for atomic potential | atomicSpacing = Distance between atoms in the lattice
	// effZ = Effective number of protons in each atom | shieldingLam = Decay factor for shielded atomic potential
	// tipRadius = Radius of curvature for metal tip/ridge | maxEnhancement = Peak enhancement factor at surface of metal
	// duration = Duration of simulation | minX = Minimum x-position | maxX = Maximum x-position | err = Expected error at 100 Up (yes, 100)
	// absEdgeSize = Width of absorptive edges | absEdgeRate = Absorptive edge rate (higher is better at absorbing but may reflect more)
	// fol = Folder where all measurements go into (include final /)
	//
	// See CPP file for measurers used.
	//
	// Typical use:
	// SimulationManager * sim = setupHHGSimulation(inputs);
	// sim->addPotential(any extra potentials desired);
	// sim->addMeasurer(any extra measurers desired);
	// sim->runPrintProgress() or sim->runParallel();
	// delete sim;
	//
	SimulationManager* setupHHGSimulation(double emax, double lam, double tau, double phase, double peakEnvT,
		double pulseSmoothTime, double ef, double w, double jellCenter, double atomCenter, double atomicSpacing,
		double effZ, double shieldingLam, double tipRadius, double maxEnhancement, double duration, double minX,
		double maxX, double err, double absEdgeSize, double absEdgeRate, const char* fol);
	//
	// Typical values:
	//
	// double emax = 20e9;					//Maximum electric field (including enhancement)
	// double lam = 800e-9;					//Wavelength
	// double tau = 8e-15;					//Light pulse duration (e^-(t/tau)^2/2)
	// double phase = 0.0;					//Light phase (0.0 -> sinusoidal, pull out then push back in)
	// double duration = 300e-15;			//Duration of simulation
	// double peakEnvT = 50e-15;			//Time in simulation at which the envelope is at its peak
	// double pulseSmoothTime = 10e-15;		//Time used to smoothly transition from no field to the gaussian
	// double ef = 9.2*PhysCon::qe;			//Fermi energy of metal
	// double w = 6.2*PhysCon::qe;			//Work function of metal
	// double jellCenter = 0.0;				//Center of jellium potential (NOT imaginary plane, but where the imaginary plane is calculated around)
	// double atomCenter = 0.0;				//Center of atomic potential
	// double atomicSpacing = 2.5e-10;		//Atomic spacing of lattice
	// double effZ = 1.74;					//Effective number of protons in atom
	// double shieldingLam = 1e-10;			//Shielding decay parameter
	// double tipRadius = 20e-9;			//Radius of curvature for tip/ridge
	// double maxEnhancement = 12.0;		//Maximum enhancement factor
	// double minX = -150e-10;				//Minimum x-position in system
	// double maxX = 700e-10;				//Maximum x-position in system
	// double err = 0.5;					//Maximum velocity error at 100 Up [from 0 (least, but infinite grid) to 2 (most) error]
	// double absEdgeSize = 100e-10;		//Width of absorptive edges
	// double absEdgeRate = 0.75;			//Rate of decay for absorptive edges

	SimulationManager* setupGasHHGSimulation(double emax, double lam, double tau, double phase, double peakEnvT,
		double pulseSmoothTime, double atomCenter, double atomicSpacing, double effZ, double shieldingLam,
		double duration, double minX, double maxX, double err, double absEdgeSize, double absEdgeRate, const char* fol);

	void addPenetrationPotential(double emax, double lam, double tau, double phase, double peakEnvT,
		double pulseSmoothTime, double ef, double w, std::complex<double> er, double cond, double jellCenter,
		double minX, double maxX, SimulationManager * sim);
}