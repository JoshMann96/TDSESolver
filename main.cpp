#include "stdafx.h"

//const char* fol = "E:/Users/Joshua/Desktop/Unsynced Stuff/PBPL Stuff/cpptesting/"; //LAPTOP
//const char* fol = "E:/Joshua/Desktop/School Stuff (Unsynced)/PBPL HHG/cpptesting/"; //DESKTOP

/*void playWithGaussianWave() {

	const char* fol = "E:/Users/Joshua/Desktop/Unsynced Stuff/PBPL Stuff/cpptesting/"; //LAPTOP
																					   //const char* fol = "E:/Joshua/Desktop/School Stuff (Unsynced)/PBPL HHG/cpptesting/"; //DESKTOP


	double dx = 5e-11;
	double dt = 5e-18;
	std::cout << "dx: " << dx << " m" << std::endl;
	std::cout << "dt: " << dt << " s" << std::endl;
	int n = (int)(700e-10 / dx);

	double maxT = 200e-15;

	SimulationManager * sim = new SingleSimulationManager(n, dx, dt, maxT);

	std::complex<double> * psi = new std::complex<double>[n];
	for (int i = 0; i < n; i++) psi[i] = std::exp(-200.0*std::pow(((double)i - n * 0.5) / n, 2) + PhysCon::im*(double)i*PhysCon::pi*0.1);
	vtls::normalizeSqrNorm(n, psi, dx);
	sim->setPsi(psi);

	double * x = vtls::linspace(n, -100e-10, 600e-10);
	sim->addPotential(new Potentials::JelliumPotential(n, x, x[n / 4], 9.2*PhysCon::qe, 6.2*PhysCon::qe, n - 1));
	sim->addPotential(new Potentials::ShieldedAtomicPotential(n, x, x[n / 4], 2.5e-10, 1.74, 1e-10));

	sim->addMeasurer(new Measurers::BasicMeasurers("test", n, dx, dt, x, fol));
	sim->addMeasurer(new Measurers::Psi2t(n, 2000, 2000, maxT, x, fol));
	sim->addMeasurer(new Measurers::OrigPot(n, fol));
	sim->addMeasurer(new Measurers::Vfunct(n, 1000, 1000, maxT, x, fol));
	sim->addMeasurer(new Measurers::ExpectX(n, x, dx, fol));
	sim->addMeasurer(new Measurers::ExpectP(n, dx, fol));
	sim->addMeasurer(new Measurers::ExpectA(n, dx, fol));

	sim->addMeasurer(new Measurers::VDProbCurrent(n, dx, 9 * n / 10, 0, "emit", fol));
	sim->addMeasurer(new Measurers::VDPsi(9 * n / 10, 0, "emit", fol));


	sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, n * 7 / 8, n, -0.5));
	sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, n / 8, 0, -0.5));


	sim->finishInitialization();

	auto t0 = std::chrono::system_clock::now();

	sim->findGroundState(100, 1e-15*PhysCon::qe);

	sim->runParallel();

	auto tf = std::chrono::system_clock::now();
	std::chrono::duration<double> dur = tf - t0;
	std::cout << "Total Time: " << dur.count() << " s" << std::endl;

	delete sim;
}

void testTiming() {
	int n = 100000;

	std::complex<double> * psi = new std::complex<double>[n];
	for (int i = 0; i < n; i++) psi[i] = std::exp(-200.0*std::pow(((double)i - n / 2.0) / n, 2) + PhysCon::im*(double)i / 5.0);

	double * v = new double[n];
	for (int i = 0; i < n; i++) {
		if (i < 0.3*n || i > 0.7*n)
			v[i] = PhysCon::qe*2000.0;
		else
			v[i] = 0.0;
	}

	double dt = 1e-24;
	double dx = 1e-12;

	TDSEIterator1D sim = TDSEIterator1D(dx, dt, n);

	std::complex<double> * psi2 = new std::complex<double>[n];

	auto t0 = std::chrono::system_clock::now();

	sim.initializeCN();

	for (int i = 0; i < 1000; i++) {
		//sim.stepEuler(psi, v, psi2);
		//sim.stepEuler(psi2, v, psi);
		sim.stepCN(psi, v, v, psi2);
		sim.stepCN(psi2, v, v, psi);
	}

	auto tend = std::chrono::system_clock::now();

	std::chrono::duration<double> dur = tend - t0;

	std::cout << dur.count() << std::endl << std::endl;
	std::cin.ignore();
}

void runHHGSimulations() {
	//All values are in SI units, radians, or unitless
	double emax = 40e9;					//Maximum electric field (including enhancement)
	double lam = 800e-9;				//Wavelength
	double tau = 8e-15;					//Light pulse duration (e^-(t/tau)^2/2)
	double phase = 0.0;					//Light phase (0.0 -> sinusoidal, pull out then push back in)
	double duration = 350e-15;			//Duration of simulation
	double peakEnvT = 50e-15;			//Time in simulation at which the envelope is at its peak
	double pulseSmoothTime = 10e-15;	//Time used to smoothly transition from no field to the gaussian
	double ef = 9.2*PhysCon::qe;		//Fermi energy of metal
	double w = 6.2*PhysCon::qe;			//Work function of metal
	double jellCenter = 0.0;			//Center of jellium potential (NOT imaginary plane, but where the imaginary plane is calculated around)
	double atomCenter = 0.0;			//Center of atomic potential
	double atomicSpacing = 2.5e-10;		//Atomic spacing of lattice
	double effZ = 1.74;					//Effective number of protons in atom
	double shieldingLam = 1e-10;		//Shielding decay parameter
	double tipRadius = 20e-9;			//Radius of curvature for tip/ridge
	double maxEnhancement = 12.0;		//Maximum enhancement factor
	double minX = -300e-10;				//Minimum x-position in system
	double maxX = 700e-10;				//Maximum x-position in system
	double err = 0.5;					//Maximum velocity error at 100 Up [from 0 (least, but infinite grid) to 2 (most) error]
	double absEdgeSize = 100e-10;		//Width of absorptive edges
	double absEdgeRate = 0.75;			//Rate of decay for absorptive edges
	const char* fol = "data/varsimsmetallicstatemixphase/";
	//------------------
	std::complex<double> ers[13] = { -23.36 + 0.77i, -31.66 + 1.0i, -40.76 + 1.26i, -50.84 + 1.53i, -61.69 + 1.98i, -73.54 + 2.27i, -86.43 + 2.86i, -100.05 + 3.55i, -114.55 + 4.03i, -130.34 + 4.53i, -147.2 + 5.42i, -164.82 + 6.41i, -183.23 + 7.52i };
	double cond = 4.1e7;
	//UCLA 1.8 Micron Experiment
	/*
	tau = 40e-15;
	peakEnvT = 150e-15;
	duration = 450e-15;
	tipRadius = 40e-9;
	ef = 5.53*PhysCon::qe;
	w = 5.1*PhysCon::qe;
	duration = 550e-15;
	peakEnvT = 200e-15;
	*//*
	//------------------
	OSSpecificFuncs::createFolder(fol);

	double lams[13] = { 800e-9, 900e-9, 1000e-9, 1100e-9, 1200e-9, 1300e-9, 1400e-9, 1500e-9, 1600e-9, 1700e-9, 1800e-9, 1900e-9, 2000e-9 };
	double emaxes[2] = { 10e9, 20e9 };

	ProgressTracker * prg = new ProgressTracker();
/*
#pragma omp parallel for collapse(2)
	for (int i = 0; i < 13; i++) {
		for (int j = 0; j < 2; j++) {
			std::stringstream ss;
			ss << fol << "lam" << (int)(lams[i] * 1e9) << "/";
			std::string str = ss.str();
			const char* nfol = str.c_str();

			OSSpecificFuncs::createFolder(nfol);

			std::stringstream ss2;
			ss2 << nfol << "emax" << (int)(emaxes[j] * 1e-9) << "/";
			std::string str2 = ss2.str();
			const char* nfol2 = str2.c_str();

			OSSpecificFuncs::createFolder(nfol2);

			SimulationManager* sim = HHGFunctions::setupHHGSimulation(emaxes[j], lams[i], tau, phase, peakEnvT, pulseSmoothTime, ef, w, jellCenter, atomCenter, atomicSpacing,
				effZ, shieldingLam, tipRadius, maxEnhancement, duration, minX, maxX, err, absEdgeSize, absEdgeRate, nfol2);

			//---------------
			//HHGFunctions::addPenetrationPotential(emaxes[j], lams[i], tau, phase, peakEnvT, pulseSmoothTime, ef, w, ers[i], cond, jellCenter, minX, maxX, sim);
			//Potentials::WaveFunctionSelfPotentialJellPotMask * npot = new Potentials::WaveFunctionSelfPotentialJellPotMask(sim->getNumPoints(), vtls::linspace(sim->getNumPoints(), minX, maxX), sim->getDX(), 1e11, atomicSpacing / 5.0, jellCenter, ef, w, sim->getNumPoints() - 1);
			//sim->addPotential(npot);
			//---------------
			sim->finishInitialization();
			//sim->findGroundState(100, 1e-15*PhysCon::qe);
			sim->findMetallicInitialState(ef, w, 50e-15, 10.0);
			//---------------
			//npot->negateGroundEffects(sim->getPsi());
			//---------------
			sim->runParallelLogProgress(prg, i * 2 + j);
			delete sim;
		}
	}
	*//*
	//---------------
	std::complex<double> ers2[3] = { -23.36 + 0.77i, -147.2 + 5.42i, -161.22 + 6.21i };
	//--------------
	double lams2[3] = { 800e-9, 1800e-9, 1880e-9 };
	double emaxes2[15] = { 10e9, 12e9, 14e9, 16e9, 18e9, 20e9, 22e9, 24e9, 26e9, 28e9, 30e9, 32e9, 34e9, 38e9, 40e9 };

#pragma omp parallel for collapse(2)
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 15; j++) {
			std::stringstream ss;
			ss << fol << "lam" << (int)(lams2[i] * 1e9) << "/";
			std::string str = ss.str();
			const char* nfol = str.c_str();

			OSSpecificFuncs::createFolder(nfol);

			std::stringstream ss2;
			ss2 << nfol << "emax" << (int)(emaxes2[j] * 1e-9) << "/";
			std::string str2 = ss2.str();
			const char* nfol2 = str2.c_str();

			OSSpecificFuncs::createFolder(nfol2);

			SimulationManager* sim = HHGFunctions::setupHHGSimulation(emaxes2[j], lams2[i], tau, phase, peakEnvT, pulseSmoothTime, ef, w, jellCenter, atomCenter, atomicSpacing,
				effZ, shieldingLam, tipRadius, maxEnhancement, duration, minX, maxX, err, absEdgeSize, absEdgeRate, nfol2);

			//---------------
			//HHGFunctions::addPenetrationPotential(emaxes2[j], lams2[i], tau, phase, peakEnvT, pulseSmoothTime, ef, w, ers2[i], cond, jellCenter, minX, maxX, sim);
			//Potentials::WaveFunctionSelfPotentialJellPotMask * npot = new Potentials::WaveFunctionSelfPotentialJellPotMask(sim->getNumPoints(), vtls::linspace(sim->getNumPoints(), minX, maxX), sim->getDX(), 1e11, atomicSpacing / 5.0, jellCenter, ef, w, sim->getNumPoints() - 1);
			//sim->addPotential(npot);
			//---------------

			sim->finishInitialization();
			//sim->findGroundState(100, 1e-15*PhysCon::qe);
			sim->findMetallicInitialState(ef, w, 50e-15, 10.0);
			//---------------
			//npot->negateGroundEffects(sim->getPsi());
			//---------------
			sim->runParallelLogProgress(prg, i * 15 + j + 26);
			delete sim;
		}
	}
	delete prg;
}

void runGasHHGSimulation() {
	//All values are in SI units, radians, or unitless
	double emax = 15e9;					//Maximum electric field (including enhancement)
	double lam = 800e-9;				//Wavelength
	double tau = 8e-15;					//Light pulse duration (e^-(t/tau)^2/2)
	double phase = 0.0;					//Light phase (0.0 -> sinusoidal, pull out then push back in)
	double duration = 300e-15;			//Duration of simulation
	double peakEnvT = 50e-15;			//Time in simulation at which the envelope is at its peak
	double pulseSmoothTime = 10e-15;	//Time used to smoothly transition from no field to the gaussian
	double atomCenter = 0.0;			//Center of atomic potential
	double atomicSpacing = 2.5e-10;		//Atomic spacing of lattice (to get atomic potential)
	double effZ = 2.5;					//Effective number of protons in atom
	double shieldingLam = 1e-10;		//Shielding decay parameter
	double minX = -400e-10;				//Minimum x-position in system
	double maxX = 400e-10;				//Maximum x-position in system
	double err = 0.4;					//Maximum velocity error at 100 Up [from 0 (least, but infinite grid) to 2 (most) error]
	double absEdgeSize = 100e-10;		//Width of absorptive edges
	double absEdgeRate = 0.75;			//Rate of decay for absorptive edges
	const char* fol = "data/gas/";

	OSSpecificFuncs::createFolder(fol);

	SimulationManager* sim = HHGFunctions::setupGasHHGSimulation(emax, lam, tau, phase, peakEnvT, pulseSmoothTime, atomCenter, atomicSpacing,
		effZ, shieldingLam, duration, minX, maxX, err, absEdgeSize, absEdgeRate, fol);

	sim->finishInitialization();
	sim->findGroundState(100, 1e-15*PhysCon::qe);
	sim->runParallel();
	delete sim;
}
*/
int main(int argc, char **argv)
{
	//runHHGSimulations();
	/*int n = 1000;
	DFTI_DESCRIPTOR_HANDLE dftiHandle = 0;
	DftiCreateDescriptor(&dftiHandle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n);
	DftiCommitDescriptor(dftiHandle);

	double* p = new double[n];
	std::complex<double>* psi = new std::complex<double>[n];
	for (int i = 0; i < n; i++) {
		psi[i] = std::exp(-(i - n / 2.0+0.5) * (i - n / 2.0+0.5) / 10.0);
	}

	vtls::normSqr(n, psi, p);
	vtlsPrnt::printGraph(n, p);

	DftiComputeForward(dftiHandle, psi);
	std::cout << "___________________________________________" << std::endl;
	vtls::normSqr(n, psi, p);
	vtlsPrnt::printGraph(n, p);
	std::cout << p[0] << " " << p[1] << " " << p[n/2-2] << " " << p[n / 2-1] << " " << p[n / 2] << " " << p[n/2+1] << " " << p[n - 1] << std::endl;
	*///vtlsPrnt::printArray(n, p);

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