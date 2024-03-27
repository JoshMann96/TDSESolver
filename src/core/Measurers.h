#pragma once
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"

// Measurers used for use with a simulation manager (or not).
namespace Measurers {

	std::fstream openFile(const char* fol);

	// Template class.
	class Measurer
	{
	private:
		int isTerminated=0;
		virtual void terminate() = 0;
	public:
		//Measurer();
		//~Measurer();
		// Required function that takes a measurement whenever called.
		virtual int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) = 0;
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
		int index = INT_MIN;
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

		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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

		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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
		NSteps(const char* fol);
		~NSteps();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records t-array.
	class TS :
		public Measurer {
	private:
		std::fstream fil;
		int index = 5;
		const char* fname = "ts.dat";
		double tmea = -1;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		TS(const char* fol);
		~TS();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records the original potential at beginning of simulation.
	class OrigPot :
		public Measurer
	{
	private:
		std::fstream fil;
		int index = 6;
		int n;
		int measured = 0;
		const char* fname = "v0.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		OrigPot(int n, const char* fol);

		~OrigPot();

		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records the wave function's probability distribution some number of times during the simulation, and also downsamples it.
	class Psi2t :
		public Measurer
	{
	private:
		std::fstream fil;
		int index = 9;
		int n;
		int nx;
		double maxT;
		int nt;
		int curIts;
		double * psi2b;
		double * psi2s;
		double * xs;
		double * ts;
		double interval;
		double mulT;
		double dt;
		const char* fname = "psi2t.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		Psi2t(int n, int nx, int nt, double maxT, double dt, double * x, const char* fol);

		~Psi2t();

		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records expectation value of energy.
	class ExpectE :
		public Measurer {
	private:
		std::fstream fil;
		int index = 10;
		const char* fname = "expectE.dat";
		int nPts;
		double* rho;
		double dx;
		void terminate();
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectE(int len, double dx, const char* fol);
		~ExpectE();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records expectation value of position.
	class ExpectX :
		public Measurer {
	private:
		std::fstream fil;
		int index = 11;
		const char* fname = "expectX.dat";
		double* x;
		int nPts;
		double* scratch;
		double dx;
		void terminate();
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectX(int len, double* xs, double dx, const char* fol);
		~ExpectX();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records expectation value of momentum (fairly computationally expensive).
	class ExpectP :
		public Measurer {
	private:
		std::fstream fil;
		int index = 12;
		const char* fname = "expectP.dat";
		int nPts;
		std::complex<double> *scratch1, *scratch2;
		double dx;
		void terminate();
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectP(int len, double dx, const char* fol);
		~ExpectP();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records expectation value of acceleration.
	class ExpectA :
		public Measurer {
	private:
		std::fstream fil;
		int index = 13;
		const char* fname = "expectA.dat";
		int nPts;
		double *scratch1, *scratch2;
		double dx;
		void terminate();
	public:
		int isHeavy() { return 1; };
		int getIndex() { return index; };
		ExpectA(int len, double dx, const char* fol);
		~ExpectA();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records the total probability remaining in simulation.
	class TotProb :
		public Measurer
	{
	private:
		double dx;
		std::fstream fil;
		double * psi2;
		int n;
		int index = 16;
		const char* fname = "tprob.dat";
		void terminate();
	public:
		int isHeavy() { return 1; };

		int getIndex() { return index; };

		TotProb(int n, double dx, const char* fol);

		~TotProb();

		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records the probability current at the virtual detector position (index).
	class VDProbCurrent :
		public Measurer {
	private:
		double dx;
		std::fstream fil;
		int n;
		int pos;
		int index = 14;
		int vdNum;
		const char* fname = "jrd.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDProbCurrent(int n, double dx, int vdPos, int vdNum, const char* name, const char* fol);
		~VDProbCurrent();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records the wave function at the virtual detector position (index).
	class VDPsi :
		public Measurer {
	private:
		std::fstream fil;
		int n;
		int pos;
		int index = 15;
		int vdNum;
		const char* fname = "psird.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDPsi(int vdPos, int vdNum, const char* name, const char* fol);
		~VDPsi();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	class VDPot :
		public Measurer {
	private:
		std::fstream fil;
		int n;
		int pos;
		int index = 21;
		int vdNum;
		const char* fname = "vrd.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDPot(int vdPos, int vdNum, const char* name, const char* fol);
		~VDPot();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	//Records the Fourier Transform in time at the VD position and one grid point to the right
	//More efficient way of storing information for creating flux spectra
	class VDFluxSpec :
		public Measurer {
	private:
		std::fstream fil;
		int pos, nElec, *nelecPtr, nsamp, first=1;
		int index = 24;
		int vdNum;
		int celec;
		double ct;
		double dw, tmax, tukeyAl=0.05;
		const char* fname = "fluxspecvd.dat";
		std::complex<double>* wfcs0, * wfcs1, *phss, cumPotPhs, *phaseCalcExpMul, *temp;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		VDFluxSpec(int vdPos, int vdNum, int* nElec, int nsamp, double emax, double tmax, const char* name, const char* fol);
		~VDFluxSpec();
		int measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records the entire wave function at sample time
	class PsiT :
		public Measurer {
	private:
		std::fstream fil;
		int n;
		double meaT;
		int index = 19;
		int vdNum;
		const char* fname = "psit.dat";
		int done = 0;
		double curTime=-1;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		PsiT(int n, double meaT, int vdNum, const char* name, const char* fol);
		~PsiT();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	// Records the potential function some number of times during the simulation, and also downsamples it.
	class Vfunct :
		public Measurer
	{
	private:
		std::fstream fil;
		int index = 17;
		int n;
		int nx;
		double maxT;
		int nt;
		int curIts;
		double * vs;
		double * xs;
		double * ts;
		double interval;
		const char* fname = "Vfunct.dat";
		void terminate();
	public:
		int isHeavy() { return 0; };

		int getIndex() { return index; };

		Vfunct(int n, int nx, int nt, double maxT, double * x, const char* fol);

		~Vfunct();

		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	class ElectronNumber :
		public Measurer {
	private:
		std::fstream fil;
		int index = 22;
		int first = 1;
		int* nElec;
		const char* fname = "nElec.dat";
		char* nfil;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		ElectronNumber(int* nElec, const char* fol);
		~ElectronNumber();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
	};

	class ExpectE0 :
		public Measurer {
	private:
		std::fstream fil;
		int index = 23;
		const char* fname = "expectE0.dat";
		int nPts;
		double dx;
		double tmea;
		double* rho;
		int first = 1;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		ExpectE0(int len, double dx, const char* fol);
		~ExpectE0();
		int measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin);
	};

	class WfcRhoWeights :
		public Measurer {
	private:
		std::fstream fil;
		int index = 25;
		const char* fname = "wghts.dat";
		int *nelecPtr, nPts;
		double dx;
		int first = 1;
		WfcToRho::Weight* wght;
		void terminate();
	public:
		int isHeavy() { return 0; };
		int getIndex() { return index; };
		WfcRhoWeights(int* nelecPtr, int nPts, double dx, WfcToRho::Weight* wght, const char* fol);
		~WfcRhoWeights();
		int measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin);
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
		BasicMeasurers(int nPts, double dx, double dt, double * xs, const char* fol);
		~BasicMeasurers();
		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);
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

		int measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin);

		int measureMany(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin, int nElec, int nPts);
	};
}