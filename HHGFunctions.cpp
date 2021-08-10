#include "stdafx.h"

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

	SimulationManager* setupHHGSimulation(double emax, double lam, double tau, double phase, double peakEnvT,
		double pulseSmoothTime, double ef, double w, double jellCenter, double atomCenter, double atomicSpacing,
		double effZ, double shieldingLam, double tipRadius, double maxEnhancement, double duration, double minX,
		double maxX, double err, double absEdgeSize, double absEdgeRate, const char* fol) {
		double pond = HHGFunctions::getPonderomotiveEnergy(emax, lam);

		double dx = HHGFunctions::getIdealDX(pond*100.0, err);
		double dt = HHGFunctions::getIdealDT(dx);

		int n = (int)((maxX - minX) / dx);

		std::cout << "SingleSimulationManager is no longer supported -- Use MultiSimulationManager with 1 simulation." << std::endl;
		throw "FEATURE REMOVED";
		//SimulationManager * sim = new SingleSimulationManager(n, dx, dt, duration);
		/*
		std::complex<double> * psi = new std::complex<double>[n];
		for (int i = 0; i < n; i++) psi[i] = 1.0;
		vtls::normalizeSqrNorm(n, psi, dx);
		sim->setPsi(psi);

		double * x = vtls::linspace(n, minX, maxX);

		sim->addPotential(new Potentials::JelliumPotentialBacked(n, x, jellCenter, ef, w, minX+absEdgeSize, absEdgeSize, n - 1));
		sim->addPotential(new Potentials::ShieldedAtomicPotential(n, x, atomCenter, atomicSpacing, effZ, shieldingLam));
		sim->addPotential(new Potentials::ElectricFieldProfileToPotential(n,
			new ElectricFieldProfiles::CylindricalToLinearProfile(n, x, getZim(ef, w) + jellCenter, maxX - absEdgeSize, tipRadius, emax, maxEnhancement),
			dx, phase, peakEnvT, lam,
			new Envelopes::SmoothedInitialGaussianEnvelope(tau, peakEnvT, pulseSmoothTime),
			n - 1));

		sim->addMeasurer(new Measurers::BasicMeasurers("test", n, dx, dt, x, fol));

		sim->addMeasurer(new Measurers::Psi2t(n, 2000, 2000, duration, dt, x, fol));
		sim->addMeasurer(new Measurers::OrigPot(n, fol));
		sim->addMeasurer(new Measurers::Vfunct(n, 1000, 1000, duration, x, fol));
		sim->addMeasurer(new Measurers::ExpectA(n, dx, fol));

		sim->addMeasurer(new Measurers::DoubleConst(emax, "emax", fol));
		sim->addMeasurer(new Measurers::DoubleConst(lam, "lam", fol));
		sim->addMeasurer(new Measurers::DoubleConst(tau, "tau", fol));
		sim->addMeasurer(new Measurers::DoubleConst(phase, "phase", fol));
		sim->addMeasurer(new Measurers::DoubleConst(peakEnvT, "peakT", fol));
		sim->addMeasurer(new Measurers::DoubleConst(ef, "fermie", fol));
		sim->addMeasurer(new Measurers::DoubleConst(w, "workf", fol));

		sim->addMeasurer(new Measurers::VDProbCurrent(n, dx, (int)((maxX - absEdgeSize - minX) / (maxX - minX)*n), 0, "emit", fol));
		sim->addMeasurer(new Measurers::VDPsi((int)((maxX - absEdgeSize - minX) / (maxX - minX)*n), 0, "emit", fol));

		sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, (int)((maxX - absEdgeSize - minX) / (maxX - minX) * n), n, -absEdgeRate));
		sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, (int)((absEdgeSize) / (maxX - minX) * n), 0, -absEdgeRate));

		return sim;
		*/
		return NULL;
	}

	SimulationManager* setupGasHHGSimulation(double emax, double lam, double tau, double phase, double peakEnvT,
		double pulseSmoothTime, double atomCenter, double atomicSpacing, double effZ, double shieldingLam, 
		double duration, double minX, double maxX, double err, double absEdgeSize, double absEdgeRate, const char* fol) {
		double pond = HHGFunctions::getPonderomotiveEnergy(emax, lam);

		double dx = HHGFunctions::getIdealDX(pond*100.0, err);
		double dt = HHGFunctions::getIdealDT(dx);

		int n = (int)((maxX - minX) / dx);

		std::cout << "SingleSimulationManager is no longer supported -- Use MultiSimulationManager with 1 simulation." << std::endl;
		throw "FEATURE REMOVED";
		//SimulationManager * sim = new SingleSimulationManager(n, dx, dt, duration);
		/*
		std::complex<double> * psi = new std::complex<double>[n];
		for (int i = 0; i < n; i++) psi[i] = 1.0;
		vtls::normalizeSqrNorm(n, psi, dx);
		sim->setPsi(psi);

		double * x = vtls::linspace(n, minX, maxX);

		sim->addPotential(new Potentials::ShieldedAtomicPotential(n, x, atomCenter, atomicSpacing, effZ, shieldingLam));
		sim->addPotential(new Potentials::ElectricFieldProfileToPotential(n,
			new ElectricFieldProfiles::ConstantFieldProfile(n, x, emax, minX+absEdgeSize, maxX-absEdgeSize),
			dx, phase, peakEnvT, lam,
			new Envelopes::SmoothedInitialGaussianEnvelope(tau, peakEnvT, pulseSmoothTime),
			n - 1));

		sim->addMeasurer(new Measurers::BasicMeasurers("test", n, dx, dt, x, fol));

		sim->addMeasurer(new Measurers::Psi2t(n, 2000, 2000, duration, dt, x, fol));
		sim->addMeasurer(new Measurers::OrigPot(n, fol));
		sim->addMeasurer(new Measurers::Vfunct(n, 1000, 1000, duration, x, fol));
		sim->addMeasurer(new Measurers::ExpectA(n, dx, fol));

		sim->addMeasurer(new Measurers::DoubleConst(emax, "emax", fol));
		sim->addMeasurer(new Measurers::DoubleConst(lam, "lam", fol));
		sim->addMeasurer(new Measurers::DoubleConst(tau, "tau", fol));
		sim->addMeasurer(new Measurers::DoubleConst(phase, "phase", fol));
		sim->addMeasurer(new Measurers::DoubleConst(peakEnvT, "peakT", fol));

		sim->addMeasurer(new Measurers::VDProbCurrent(n, dx, (int)((maxX - absEdgeSize - minX) / (maxX - minX)*n), 0, "emit", fol));
		sim->addMeasurer(new Measurers::VDPsi((int)((maxX - absEdgeSize - minX) / (maxX - minX)*n), 0, "emit", fol));

		sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, (int)((maxX - absEdgeSize - minX) / (maxX - minX) * n), n, -absEdgeRate));
		sim->addTimeRotation(AbsorptiveRegions::getSmoothedTimePhaseDecay(n, (int)((absEdgeSize) / (maxX - minX) * n), 0, -absEdgeRate));

		return sim;
		*/
		return NULL;
	}

	void addPenetrationPotential(double emax, double lam, double tau, double phase, double peakEnvT,
		double pulseSmoothTime, double ef, double w, std::complex<double> er, double cond, double jellCenter,
		double minX, double maxX, SimulationManager * sim) {
		double * x = vtls::linspace(sim->getNumPoints(), minX, maxX);
		sim->addPotential(new Potentials::ElectricFieldProfileToPotential(sim->getNumPoints(),
			new ElectricFieldProfiles::InMetalFieldProfile(sim->getNumPoints(), x, minX, getZim(ef, w) + jellCenter, emax, lam, er, cond),
			sim->getDX(), phase, peakEnvT, lam,
			new Envelopes::SmoothedInitialGaussianEnvelope(tau, peakEnvT, pulseSmoothTime),
			sim->getNumPoints() - 1));
	}
}