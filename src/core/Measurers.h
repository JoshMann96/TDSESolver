#pragma once
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"
#include "MathTools.h"


// Measurers used for use with a simulation manager (or not).
namespace Measurers {

	std::fstream openFile(const char* fol);

	// Template class.
	class Measurer
	{
	private:
		int isTerminated=0;
		virtual void terminate() = 0;
	public:;
		virtual ~Measurer() = default;
		// Required function that takes a measurement whenever called.
		virtual int measure(int step, std::complex<double> * psi, double * v, double t) = 0;
		// Required function which terminates the measurer (closes file).
		void kill(){
			if(!isTerminated){
				isTerminated=1;
				terminate();
			}
		}
		virtual int getIndex() = 0;
		virtual int isHeavy() = 0;
	private:
		int index;
		const char* fname = "";
	};

	// Writes a double constant to a file for future reference.
	class DoubleConst :
		public Measurer
	{
	private:
		double c;
		std::fstream fil;
		int index = -2;
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex(){ return index; };

		DoubleConst(double c, const char* filName, const char* fol);

		~DoubleConst();

		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Writes text (8 chars required) to file.
	class Header :
		public Measurer
	{
	private:
		std::fstream fil;
		int index = -1;
		const char* fname = "head.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		Header(const char* title, const char* fol);

		~Header();

		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the number of points in simulation.
	class NPts :
		public Measurer {
	private:
		std::fstream fil;
		int index = 0;
		const char* fname = "nPts.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };
		NPts(int nPts, const char* fol);
		~NPts();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the number of time steps in simulation.
	class NSteps :
		public Measurer {
	private:
		std::fstream fil;
		int index = 1;
		const char* fname = "nSteps.dat";
		int steps = 0;
		double tmea = -1;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		NSteps(int numSteps, const char* fol);
		~NSteps();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records dx spacing.
	class DX :
		public Measurer {
	private:
		std::fstream fil;
		int index = 2;
		const char* fname = "dx.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		DX(double dx, const char* fol);
		~DX();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records dt spacing.
	class DT :
		public Measurer {
	private:
		std::fstream fil;
		int index = 3;
		const char* fname = "dt.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		DT(double dt, const char* fol);
		~DT();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records x-array.
	class XS :
		public Measurer {
	private:
		std::fstream fil;
		int index = 4;
		const char* fname = "xs.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		XS(int len, double* xs, const char* fol);
		~XS();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records t-array.
	class TS :
		public Measurer {
	private:
		std::fstream fil;
		int index = 5;
		const char* fname = "ts.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		TS(const char* fol);
		~TS();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the original potential at beginning of simulation.
	class OrigPot :
		public Measurer
	{
	private:
		std::fstream fil;
		int index = 6;
		int n;
		const char* fname = "v0.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		OrigPot(int n, const char* fol);

		~OrigPot();

		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the wave function's probability distribution some number of times during the simulation, and also downsamples it.
	class Psi2t :
		public Measurer
	{
	private:
		std::fstream fil;
		int index = 9;
		int nPts;
		int *nElec;
		int nx, nt;
		int numSteps;
		int curIdx;
		double * psi2b;
		double * psi2s;
		double * xs;
		double * ts;
		const char* fname = "psi2t.dat";

		int *measSteps;

		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		Psi2t(int nPts, int nx, int nt, int numSteps, double maxT, double * x, int* nElec, const char* fol);

		~Psi2t();

		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records expectation value of energy.
	class ExpectE :
		public Measurer {
	private:
		std::fstream fil;
		int index = 10;
		const char* fname = "expectE.dat";
		int nPts, *nElec;
		double* rho;
		double dx;
		void terminate();
		KineticOperators::KineticOperator** kin;
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectE(int len, double dx, int* nElec, const char* fol, KineticOperators::KineticOperator** kin);
		~ExpectE();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records expectation value of position.
	class ExpectX :
		public Measurer {
	private:
		std::fstream fil;
		int index = 11;
		const char* fname = "expectX.dat";
		double* x;
		int nPts, *nElec;
		double* scratch;
		double dx;
		void terminate();
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectX(int len, double* xs, double dx, int* nElec, const char* fol);
		~ExpectX();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records expectation value of momentum (fairly computationally expensive).
	class ExpectP :
		public Measurer {
	private:
		std::fstream fil;
		int index = 12;
		const char* fname = "expectP.dat";
		int nPts, *nElec;
		std::complex<double> *scratch1, *scratch2;
		double dx;
		void terminate();
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectP(int len, double dx, int* nElec, const char* fol);
		~ExpectP();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records expectation value of acceleration.
	class ExpectA :
		public Measurer {
	private:
		std::fstream fil;
		int index = 13;
		const char* fname = "expectA.dat";
		int nPts, *nElec;
		double *scratch1, *scratch2;
		double dx;
		void terminate();
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectA(int nPts, double dx, int* nElec, const char* fol);
		~ExpectA();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the total probability remaining in simulation.
	class TotProb :
		public Measurer
	{
	private:
		double dx;
		std::fstream fil;
		double * psi2;
		int nPts, *nElec;
		int index = 16;
		const char* fname = "totProb.dat";
		void terminate();
	public:
		int isHeavy() { return 1; };

		int getIndex() { return index; };

		TotProb(int n, double dx, int* nElec, const char* fol);

		~TotProb();

		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the probability current at the virtual detector position (index).
	class VDProbCurrent :
		public Measurer {
	private:
		double dx;
		std::fstream fil;
		int nPts, *nElec;
		int vdPos;
		int index = 14;
		int vdNum;
		const char* fname = "jrd.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDProbCurrent(int n, double dx, int *nElec, int vdPos, int vdNum, const char* name, const char* fol);
		~VDProbCurrent();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the wave function at the virtual detector position (index).
	class VDPsi :
		public Measurer {
	private:
		std::fstream fil;
		int nPts, *nElec;
		int vdPos;
		int index = 15;
		int vdNum;
		const char* fname = "psird.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDPsi(int* nElec, int vdPos, int vdNum, const char* name, const char* fol);
		~VDPsi();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	class VDPot :
		public Measurer {
	private:
		std::fstream fil;
		int n;
		int vdPos;
		int index = 21;
		int vdNum;
		int curStep = -1;
		const char* fname = "vrd.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDPot(int vdPos, int vdNum, const char* name, const char* fol);
		~VDPot();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	//Records the Fourier Transform in time at the VD position and one grid point to the right
	//More efficient way of storing information for creating flux spectra
	class VDFluxSpec :
		public Measurer {
	private:
		std::fstream fil;
		int vdPos, *nElec, nSamp, first=1;
		int index = 24;
		int vdNum, nPts;
		double ct;
		double dw, tmax, tukeyAl=0.05;
		const char* fname = "fluxspecvd.dat";
		std::complex<double>* wfcs0 = nullptr, * wfcs1 = nullptr, *phss, cumPotPhs, *phaseCalcExpMul, *temp;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDFluxSpec(int nPts, int vdPos, int vdNum, int* nElec, int nsamp, double emax, double tmax, const char* name, const char* fol);
		~VDFluxSpec();
		int measure(int step, std::complex<double>* psi, double* v, double t);
	};

	// Records the entire wave function at sample time
	class PsiT :
		public Measurer {
	private:
		std::fstream fil;
		int nPts;
		double meaT;
		int index = 19;
		int vdNum, *nElec;
		const char* fname = "psit.dat";
		int done = 0;
		double curTime=-1;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		PsiT(int n, double meaT, int* nElec, int vdNum, const char* name, const char* fol);
		~PsiT();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the entire potential at sample time
	class PotT :
		public Measurer {
	private:
		std::fstream fil;
		int n;
		double meaT;
		int index = 20;
		int vdNum;
		const char* fname = "pott.dat";
		int done = 0;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		PotT(int n, double meaT, int vdNum, const char* name, const char* fol);
		~PotT();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Records the potential function some number of times during the simulation, and also downsamples it.
	class Vfunct :
		public Measurer
	{
	private:
		std::fstream fil;
		int index = 17;
		int nPts;
		int nx;
		double maxT;
		int nt;
		int curIdx;
		int* measSteps;
		double * vs;
		double * xs;
		double * ts;
		const char* fname = "Vfunct.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		Vfunct(int nPts, int nx, int nt, int numSteps, double maxT, double * x, const char* fol);

		~Vfunct();

		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	class NElec :
		public Measurer {
	private:
		std::fstream fil;
		int index = 22;
		int first = 1;
		int* nElec;
		const char* fname = "nElec.dat";
		char* nfil = nullptr;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		NElec(int* nElec, const char* fol);
		~NElec();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	class ExpectE0 :
		public Measurer {
	private:
		std::fstream fil;
		int index = 23;
		const char* fname = "expectE0.dat";
		int nPts, *nElec;
		double dx;
		double tmea;
		double* rho;
		int first = 1;
		void terminate();
		KineticOperators::KineticOperator** kin;
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		ExpectE0(int nPts, double dx, int* nElec, const char* fol, KineticOperators::KineticOperator** kin);
		~ExpectE0();
		int measure(int step, std::complex<double>* psi, double* v, double t);
	};

	class WfcRhoWeights :
		public Measurer {
	private:
		std::fstream fil;
		int index = 25;
		const char* fname = "wghts.dat";
		int *nElec, nPts;
		double dx;
		int first = 1;
		WfcToRho::Weight* wght;
		KineticOperators::KineticOperator ** kin;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		WfcRhoWeights(int* nelecPtr, int nPts, double dx, WfcToRho::Weight* wght, KineticOperators::KineticOperator** kin, const char* fol);
		~WfcRhoWeights();
		int measure(int step, std::complex<double>* psi, double* v, double t);
	};

	// Includes several basic measurements. See cpp file for specifics.
	class BasicMeasurers :
		public Measurer {
	private:
		std::vector<Measurer*> meas;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return INT_MIN + 1; };
		BasicMeasurers(int nPts, int numSteps, double dx, double dt, double * xs, const char* fol);
		~BasicMeasurers();
		int measure(int step, std::complex<double> * psi, double * v, double t);
	};

	// Manages multiple measurements. Ideal if using more than one.
	class MeasurementManager :
		public Measurer {
	private:
		std::fstream fil;
		int index = INT_MAX;
		std::vector<Measurer*> meas;
		const char* fname;
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		MeasurementManager(const char* fname);

		~MeasurementManager();

		void addMeasurer(Measurer * m);

		int measure(int step, std::complex<double> * psi, double * v, double t);
	};
}