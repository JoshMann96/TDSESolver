#pragma once
#include "CORECommonHeader.h" 

namespace Potentials {
	// Electric field profiles for use of creating potentials.
	namespace ElectricFieldProfiles {
		// Template class (no field).
		class ElectricFieldProfile {
		public:
			virtual std::complex<double> * getProfile() = 0;
		};

		// Constant field between minPos and maxPos of strength e.
		class ConstantFieldProfile :
			public ElectricFieldProfile {
		private:
			std::complex<double> * fs;
		public:
			ConstantFieldProfile(int nPts, double * x, double eMax, double minX, double maxX);
			std::complex<double> * getProfile();
		};

		// Uses a cylindrical profile (enhFact*r/x+1) from minX until a linear approximation leads to the field dying out at maxX, at which point the linear approximation is used.
		// It is important for there to be no field at the edge of the simulation (and ideal that it smoothyl dies out, like this linear approximation).
		// Otherwise we may see extraneous plateaus, for instance, in HHG.
		class CylindricalToLinearProfile :
			public ElectricFieldProfile {
		private:
			std::complex<double> * fs;
		public:
			CylindricalToLinearProfile(int nPts, double * x, double minX, double maxX, double r, double eMax, double enhFact);
			~CylindricalToLinearProfile();
			std::complex<double> * getProfile();
		};

		class CylindricalToCutoffProfile :
			public ElectricFieldProfile {
		private:
			std::complex<double> * fs;
		public:
			CylindricalToCutoffProfile(int nPts, double * x, double minX, double maxX, double r, double eMax, double enhFact, double decayLength);
			~CylindricalToCutoffProfile();
			std::complex<double> * getProfile();
		};

		class InMetalFieldProfile :
			public ElectricFieldProfile {
		private:
			std::complex<double> * fs;
		public:
			InMetalFieldProfile(int nPts, double * x, double minX, double maxX, double eMax, double lam, std::complex<double> er, double cond);
			~InMetalFieldProfile();
			std::complex<double> * getProfile();
		};

		class ExponentialToLinearProfile :
			public ElectricFieldProfile {
		private:
			std::complex<double>* fs;
		public:
			ExponentialToLinearProfile(int nPts, double* x, double minX, double maxX, double r, double eMax);
			~ExponentialToLinearProfile();
			std::complex<double>* getProfile();
		};

		class FileFieldProfile :
			public ElectricFieldProfile {
			std::complex<double> * fs;
		public:
			FileFieldProfile(int nPts, double * x, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, const char * fil);
			~FileFieldProfile();
			std::complex<double> * getProfile();
		};
	}

	// Envelopes for use with laser pulses.
	namespace Envelopes {
		// Template class (no pulse).
		class Envelope {
		public:
			virtual double getValue(double t) = 0;
		};

		class GaussianEnvelope :
			public Envelope {
		private:
			double tau;
			double tmax;
		public:
			GaussianEnvelope(double tau, double tmax);
			double getValue(double t);
		};

		//Smooths the beginning of the envelope from zero to gaussian using a 13th order smooth polynomial. Improves stability of the beginning of the simulation.
		class SmoothedInitialGaussianEnvelope :
			public Envelope {
		private:
			double tau;
			double tmax;
			double buf;
		public:
			SmoothedInitialGaussianEnvelope(double tau, double tmax, double bufferTime);
			double getValue(double t);
		};

		class CosSquaredEnvelope :
			public Envelope {
		private:
			double tau;
			double tmax;
			double a_t = 2.0 * std::acos(std::pow(0.5, 0.25)); // cos^4(a_t*(tau/2)/tau)=1/2 s.t. tau = FWHM-power
		public:
			CosSquaredEnvelope(double tau, double tmax);
			double getValue(double t);
		};
	}

	enum PotentialComplexity{
		STATIC, // potential is a constant, only needs to be evaluated once
		DYNAMIC, // potential changes with time but is independent of the wavefunction
		WAVEFUNCTION_DEPENDENT // the potential depends on the wavefunction (including density-functional potentials)
	};

	// Template function for potential (zero potential).
	class Potential {
	public:
		virtual ~Potential() = default;
		virtual void getV(double t, double * targ) = 0;
		virtual void getV(double * rho, std::complex<double> * psi, double t, double * targ) = 0;
		virtual PotentialComplexity getComplexity() = 0;
	};

	//void calcFermiBoxDimensionalityConversion(int nElec, int nPts, double dx, double ef, std::complex<double>* psi, Potential* totPot, double* prefactor);

	// Read potential from binary data file.
	// The first element of the file should be an integer representing the number of points provided in the potential (n).
	// The next n doubles are the positions that these are sampled at, in m.
	// The next n doubles are the potentials, in J.
	class FilePotential :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		FilePotential(int nPts, double * x, double offset, const char * fil, int refPoint);
		~FilePotential();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	class BiasFieldPotential :
		public Potential
	{
	private:
		double * v;
		double tstart, tbuf;
		int nPts;
	public:
		BiasFieldPotential(int nPts, double * x, double tstart, double tbuf, double xmin, double xmax, double xmin_buf, double xmax_buf, double fieldStrength, int refPoint);
		~BiasFieldPotential();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> *  psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	class CoulombPotential :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		CoulombPotential(int nPts, double * x, double ne, double chargePos, double minX, double maxX, int refPoint);
		~CoulombPotential();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> *  psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	class FiniteBox :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		FiniteBox(int nPts, double * x, double left, double right, double vin, int refPoint);
		~FiniteBox();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	// Wacther's Jellium potential.
	class JelliumPotential :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		JelliumPotential(int nPts, double * x, double center, double ef, double w, int refPoint);
		~JelliumPotential();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	class JelliumPotentialBacked :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		JelliumPotentialBacked(int nPts, double * x, double center, double ef, double w, double backStart, double backWidth, int refPoint);
		~JelliumPotentialBacked();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	// Shielded atomic potential, averaged across an infinite plane parallel to metal's surface.
	class ShieldedAtomicPotential :
		public Potential {
	private:
		double * v;
		int nPts;
	public:
		ShieldedAtomicPotential(int nPts, double * x, double center, double latticeSpacing, double zProtons, double decayConst);
		~ShieldedAtomicPotential();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	// Converts an electric field profile and envelope to a light pulse potential.
	class ElectricFieldProfileToPotential :
		public Potential
	{
	private:
		int nPts;
		std::complex<double> * potMask;
		double phase;
		double tmax;
		double w;
		Envelopes::Envelope * env;
	public:
		ElectricFieldProfileToPotential(int nPts, ElectricFieldProfiles::ElectricFieldProfile * fieldProfile, double dx, double phase, double tmax, double lam, Envelopes::Envelope * env, int refPoint);
		~ElectricFieldProfileToPotential();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::DYNAMIC;};
	};

	// Class of potentials which represent relative changes with respect to the initital state.
	// In other words, it is intended to be zero at t=0 and represent the change in potential from the initial state.
	class NonlinearDynamicalPotential :
		public Potential
	{
	private:
	public:
		virtual void negateGroundEffects(double* rho, std::complex<double> * psi) = 0;
		virtual void getV(double t, double * targ) = 0;
		virtual void getV(double* rho, std::complex<double> * psi, double t, double * targ) = 0;
		virtual PotentialComplexity getComplexity() = 0;
	};

	class CurrentIntegrator
	{
	private:
		double integratedFlux, last_t;
		double ** weights, dx;
		int evalPoint, * nElec, nPts;
	public:
		CurrentIntegrator(int nPts, double dx, int evalPoint,  int* nElec, double** weights);
		void integrate(std::complex<double>* psi, double t);
		double getIntegratedFlux(){return integratedFlux;};
	};

	class CylindricalImageCharge :
		public NonlinearDynamicalPotential
	{
	private:
		int nPts, refPoint, posMin, posMax, surfPos;
		double dx, ef, w, rad, * origPot, * potTemp, * genTemp, * lrxr, * myRho, *nsMask, *dethin;
		void calcPot(double* rho, std::complex<double>* psi, double cur_t, double* targ);
		CurrentIntegrator * curInt;
	public:
		CylindricalImageCharge(int nPts, double* x, double dx, double ef, double w, double rad, int* nElec, double** weights, int posMin, int posMax, int surfPos, int refPoint);
		~CylindricalImageCharge();
		void negateGroundEffects(double* rho, std::complex<double>* psi);
		void getV(double t, double* targ);
		void getV(double* rho, std::complex<double>* psi, double t, double* targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::WAVEFUNCTION_DEPENDENT;};
	};

	class PlanarToCylindricalHartree :
		public NonlinearDynamicalPotential
	{
	private:
		int nPts, refPoint, * nElec, posMin, posMax, surfPos;
		double dx, rad, * origPot, * potTemp, * fieldScaler, *myRho, *dethin;
		double originalCharge;
		void calcPot(double* rho, std::complex<double>* psi, double t, double* targ);
		CurrentIntegrator * curInt;
	public:
		PlanarToCylindricalHartree(int nPts, double* x, double dx, double rad, int* nElec, double** weights, int posMin, int posMax, int surfPos, int refPoint);
		~PlanarToCylindricalHartree();
		void negateGroundEffects(double* rho, std::complex<double>* psi);
		void getV(double t, double* targ);
		void getV(double* rho, std::complex<double>* psi, double t, double* targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::WAVEFUNCTION_DEPENDENT;};
	};

	enum class LDAFunctionalType {
		X_SLATER,
		C_PW
	};
	
	class LDAFunctional :
		public NonlinearDynamicalPotential
	{
	private:
		int nPts, refPoint, * nElec;
		double * origPot, * rho, dx;
		void calcPot(double* rho, double* targ);

		LDAFunctionalType typ;
	public:
		LDAFunctional(LDAFunctionalType typ, int nPts, double dx, int refPoint);
		~LDAFunctional();
		void negateGroundEffects(double* rho, std::complex<double>* psi);
		void getV(double t, double* targ);
		void getV(double* rho, std::complex<double>* psi, double t, double* targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::WAVEFUNCTION_DEPENDENT;};
	};

	// Combines multiple potentials of all kinds.
	class CompositePotential :
		public Potential {
	private:
		int nPts;
		int numSPots;
		int numDPots;
		int numWPots;
		Potential ** staticPots;
		Potential ** dynamicPots;
		Potential ** waveFuncDependentPots;
		double * v0;
		double * nv;
	public:
		CompositePotential(int nPts, int numSPots, int numDPots, int numWPots, Potential ** staticPots, Potential ** dynamicPots, Potential ** waveFuncDependentPots);
		~CompositePotential();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity();
	};

	// Manages potentials.
	class PotentialManager :
		public Potential {
	private:
		int nPts;
		std::vector<Potential*> staticPots, dynamicPots, waveFuncDependentPots;
		CompositePotential * pot=nullptr;
		PotentialComplexity myComplex = PotentialComplexity::STATIC;
		Potential ** spots = nullptr, ** dpots = nullptr, ** wpots = nullptr;
	public:
		PotentialManager(int nPts);
		~PotentialManager(){if(pot) delete pot; if(spots) delete[] spots; if(dpots) delete[] dpots; if(wpots) delete[] wpots;};
		void addPotential(Potential * pot);
		void finishAddingPotentials();
		void getV(double t, double * targ);
		void getV(double* rho, std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return myComplex;};
	};
}