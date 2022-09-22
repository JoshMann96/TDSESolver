#pragma once
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
		std::complex<double> * getProfile();
	};

	class CylindricalToCutoffProfile :
		public ElectricFieldProfile {
	private:
		std::complex<double> * fs;
	public:
		CylindricalToCutoffProfile(int nPts, double * x, double minX, double maxX, double r, double eMax, double enhFact, double decayLength);
		std::complex<double> * getProfile();
	};

	class InMetalFieldProfile :
		public ElectricFieldProfile {
	private:
		std::complex<double> * fs;
	public:
		InMetalFieldProfile(int nPts, double * x, double minX, double maxX, double eMax, double lam, std::complex<double> er, double cond);
		std::complex<double> * getProfile();
	};

	class ExponentialToLinearProfile :
		public ElectricFieldProfile {
	private:
		std::complex<double>* fs;
	public:
		ExponentialToLinearProfile(int nPts, double* x, double minX, double maxX, double r, double eMax);
		std::complex<double>* getProfile();
	};

	class FileFieldProfile :
		public ElectricFieldProfile {
		std::complex<double> * fs;
	public:
		FileFieldProfile(int nPts, double * x, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, const char * fil);
		virtual std::complex<double> * getProfile();
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
}

namespace Potentials {
	// Template function for potential (zero potential).
	class Potential {
	public:
		virtual void getV(double t, double * targ, KineticOperators::KineticOperator* kin) = 0;
		virtual void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin) = 0;
		virtual int isDynamic() = 0;
	};

	//void calcFermiBoxDimensionalityConversion(int nelec, int nPts, double dx, double ef, std::complex<double>* psi, Potential* totPot, double* prefactor);

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
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
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
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> *  psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class CoulombPotential :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		CoulombPotential(int nPts, double * x, double ne, double chargePos, double minX, double maxX, int refPoint);
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> *  psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class FiniteBox :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		FiniteBox(int nPts, double * x, double left, double right, double vin, int refPoint);
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
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
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class JelliumPotentialBacked :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		JelliumPotentialBacked(int nPts, double * x, double center, double ef, double w, double backStart, double backWidth, int refPoint);
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	// Shielded atomic potential, averaged across an infinite plane parallel to metal's surface.
	class ShieldedAtomicPotential :
		public Potential {
	private:
		double * v;
		int nPts;
	public:
		ShieldedAtomicPotential(int nPts, double * x, double center, double latticeSpacing, double zProtons, double decayConst);
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
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
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	// Creates a pseudo multi-electron potential based on coulomb's law. Can negate ground state effects and so only changes in probability distribution alters the potential.
	class WaveFunctionSelfPotential :
		public Potential
	{
	private:
		/*int nPts;
		double strength;
		double dx;
		double odimDist;
		int refPoint;
		double * psi2;
		double * mask;
		double * add;
		VSLConvTaskPtr task, *task_ptr;
		int status;
		void initializeConvolution();*/
	public:
		//WaveFunctionSelfPotential(int nPts, double dx, double strength, double otherDimensionDistance, int refPoint);
		virtual void negateGroundEffects(std::complex<double> * psi, KineticOperators::KineticOperator* kin) = 0;
		virtual void getV(double t, double * targ, KineticOperators::KineticOperator* kin) = 0;
		virtual void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin) = 0;
		virtual int isDynamic() = 0;
	};

	//Does as described above, however only outside metal. Does not do real charge distribution...
	class WaveFunctionSelfPotentialJellPotMask :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts;
		double strength;
		double dx;
		double odimDist;
		int refPoint;
		double* psi2;
		double* mask;
		double* stepFunc;
		double* add;
		VSLConvTaskPtr task, * task_ptr;
		int status;
		void initializeConvolution();
	public:
		WaveFunctionSelfPotentialJellPotMask(int nPts, double* x, double dx, double strength, double otherDimensionDistance, double center, double ef, double w, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class SurfaceSpaceCharge :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts;
		double dx, ef;
		int refPoint;
		int* nelecPtr;
		int nelec;
		int first = 1;
		double* prefactor, *rho;
		double* origPot;
		int posMin, posMax;
		void calcPot(std::complex<double>* psi, double* targ, KineticOperators::KineticOperator* kin);
		void doFirst(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		Potential* totPot;
		WfcToRho::Weight* wght;
		WfcToRho::Density* dens;
		double* fldTot, *psi2;
		std::complex<double>* psi0;
	public:
		SurfaceSpaceCharge(int nPts, double dx, double ef, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class FullCylindricalSpaceCharge :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts, refPoint, *nelecPtr, nelec, first=1, posMin, posMax, surfPos;
		double dx, ef, r, *prefactor, *origPot, *potTemp, *psi2, *lrxr, *rho;
		void calcPot(std::complex<double>* psi, double* targ, KineticOperators::KineticOperator* kin);
		void doFirst(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		Potential* totPot;
		WfcToRho::Weight* wght;
		WfcToRho::Density* dens;
	public:
		FullCylindricalSpaceCharge(int nPts, double * x, double dx, double ef, double r, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int surfPos, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class LinearBulkCylindricalFieldSpaceCharge :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts, refPoint, * nelecPtr, nelec, first = 1, posMin, posMax, surfPos, trackInnerLoss;
		double dx, dt, ef, rad, * prefactor, * origPot, * potTemp, * genTemp, * psi2, * lrxr, * rho;
		//variables for charge loss prevention
		double lostCharge = 0.0, chargeCenterD, chargeWidthD;
		int chargeCenter, chargeWidth;
		void calcPot(std::complex<double>* psi, double* targ, KineticOperators::KineticOperator* kin);
		void doFirst(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		Potential* totPot;
		WfcToRho::Weight* wght;
		WfcToRho::Density* dens;
	public:
		LinearBulkCylindricalFieldSpaceCharge(int nPts, double* x, double dx, double dt, double ef, double rad, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int surfPos, int trackInnerLoss, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class CylindricalImageCharge :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts, refPoint, * nelecPtr, nelec, first = 1, posMin, posMax, surfPos;
		double dx, dt, ef, w, rad, * prefactor, * origPot, * potTemp, * genTemp, * psi2, * lrxr, * rho, *nsMask;
		double emittedCharge = 0.0;
		void calcPot(std::complex<double>* psi, double* targ, KineticOperators::KineticOperator* kin);
		void doFirst(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		Potential* totPot;
		WfcToRho::Weight* wght;
		WfcToRho::Density* dens;
	public:
		CylindricalImageCharge(int nPts, double* x, double dx, double dt, double ef, double w, double rad, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int surfPos, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class DielectricBulkCylindricalFieldSpaceCharge :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts, refPoint, * nelecPtr, nelec, first = 1, posMin, posMax, surfPos;
		double dx, dt, ef, rad, * prefactor, * origPot, * potTemp, * genTemp, * psi2, * lrxr, *xs, * rho, wellWidth, dampRate, *xfsf, *xfsf2;
		//double pol, lpol;
		void calcPot(std::complex<double>* psi, double* targ, KineticOperators::KineticOperator* kin);
		void doFirst(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		Potential* totPot;
		WfcToRho::Weight* wght;
		WfcToRho::Density* dens;
	public:
		DielectricBulkCylindricalFieldSpaceCharge(int nPts, double* x, double dx, double dt, double ef, double rad, double wellWidth, double dampRate, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int posMin, int posMax, int surfPos, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class LinearBulkCylSectionFieldSpaceCharge :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts, refPoint, * nelecPtr, nelec, first = 1, surfPos;
		double dx, ef, rad, * prefactor, * origPot, * potTemp, * genTemp, * psi2, * rho, *mMat;
		void calcPot(std::complex<double>* psi, double* targ, KineticOperators::KineticOperator* kin);
		void doFirst(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		Potential* totPot;
		WfcToRho::Weight* wght;
		WfcToRho::Density* dens;
	public:
		LinearBulkCylSectionFieldSpaceCharge(int nPts, double* x, double dx, double ef, double rad, double theta0, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int surfPos, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	class OhmicRetardingPotential :
		public WaveFunctionSelfPotential
	{
	private:
		int nPts, refPoint, * nelecPtr, nelec, first=1;
		double dx, * prefactor, * probCur, *mask, resistivity, *origPot;
		std::complex<double>* temp;
		void calcProbCur(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void doFirst(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void calcPot(std::complex<double>* psi, double* targ, KineticOperators::KineticOperator* kin);
		Potential* totPot;
		WfcToRho::Weight* wght;
		WfcToRho::Density* dens;
	public:
		OhmicRetardingPotential(int nPts, double dx, double transLen, double resistivity, int* nelec, Potential* totPot, WfcToRho::Weight* wght, WfcToRho::Density* dens, int surfPos, int refPoint);
		void negateGroundEffects(std::complex<double>* psi, KineticOperators::KineticOperator* kin);
		void getV(double t, double* targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double>* psi, double t, double* targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
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
		CompositePotential(int nPts, int numSPots, int numDPots, int numWPots, KineticOperators::KineticOperator* kin, Potential ** staticPots, Potential ** dynamicPots, Potential ** waveFuncDependentPots);
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};

	// Manages potentials.
	class PotentialManager :
		public Potential {
	private:
		int nPts;
		std::vector<Potential*> staticPots;
		std::vector<Potential*> dynamicPots;
		std::vector<Potential*> waveFuncDependentPots;
		CompositePotential * pot;
	public:
		PotentialManager(int nPts);
		void addPotential(Potential * pot);
		void finishAddingPotentials(KineticOperators::KineticOperator* kin);
		void getV(double t, double * targ, KineticOperators::KineticOperator* kin);
		void getV(std::complex<double> * psi, double t, double * targ, KineticOperators::KineticOperator* kin);
		int isDynamic();
	};
}