/**
 * @file Measurers.h
 * @brief Contains classes for measuring various properties of the wave function and potential.
 */
#pragma once
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"
#include "MathTools.h"


/**
 * @namespace Measurers
 * @brief Contains classes for measuring various properties of the wave function and potential.
 */
namespace Measurers {

	/// Status of a Measurer::measure operation.
	enum MeasurerStatus {
		/// The measurement was successful.
		SUCCESS,
		/// The measurement failed (error).
		FAIL,
		/// The Measurer::measure call was the last and the Measurer is ready to be destructed.
		ALL_DONE
	};

	/**
	 * Opens a file fstream for writing binary data.
	 * @param fil The path to the file to open.
	 * @return The file stream.
	 */
	std::fstream openFile(const char* fil);

	/**
	 * Opens a file fstream for writing binary data.
	 * @param args The path to the file to open. Must be a list of strings to be concatenated.
	 * @return The file stream.
	 */
	std::fstream openFile(std::initializer_list<const char*> args);

	/// @copydoc openFile(std::initializer_list<const char*> args)
	std::fstream openFile(std::list<const char*> args);

	/// Template for a measurer class.
	class Measurer
	{
	protected:
		bool needsDens = false;
		std::fstream fil;
		const char* ext = ".dat";
		int index;
	public:;
		/// Default constructor.
		Measurer() = default;

		/**
		 * Constructor.
		 * @param index The index of the measurer.
		 */
		Measurer(int index) : index(index) {};

		/**
		 * Constructor.
		 * @param index The index of the measurer.
		 * @param fname The name of the file to write to. The extension '.dat' will be appended and the fstream \a fil will be opened.
		 */
		Measurer(int index, const char* fname) : index(index) {
			if (fname){
				fil = openFile({fname, ext});
				fil.write(reinterpret_cast<char*>(&index), sizeof(int));
			}
		}

		/**
		 * Constructor.
		 * @param index The index of the measurer.
		 * @param fnameArgs The path to the file to open. Must be a list of strings to be concatenated.
		 */
		Measurer(int index, std::initializer_list<const char*> fnameArgs) : index(index) {
			open(fnameArgs);
		}

		/// Destructor, closes the fstream if it is open.
		~Measurer(){
			if(fil.is_open())
				fil.close();
		}

		/**
		 * Opens a file fstream for writing binary data.
		 * @param fnameArgs The path to the file to open. Must be a list of strings to be concatenated.
		 */
		void open(std::initializer_list<const char*> fnameArgs){
			if(fil.is_open()){
				std::cerr << "Warning: File stream already open! Closing it before opening a new one." << std::endl;
				fil.close();
			}

			std::list<const char*> args(fnameArgs);
			args.push_back(ext);
			fil = openFile(args);
			fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		}

		/**
		 * Writes data to the file.
		 * @param data The data to write.
		 * @param size The size of the data to write.
		 */
		void write(void* data, size_t size){
			fil.write(reinterpret_cast<char*>(data), size);
		}

		/// @copydoc write(void* data, size_t size)
		void write(const void* data, size_t size){
			fil.write(reinterpret_cast<const char*>(data), size);
		}

		/**
		 * Closes the file stream.
		 * This is called in the destructor, but can be called manually if needed.
		 */
		void close(){
			if(fil.is_open())
				fil.close();
		}
		
		/**
		 * Performs a measurement.
		 * @param step The current time step.
		 * @param psi (in) The wave function.
		 * @param v (in) The potential.
		 * @param t The current time.
		 * @return The status of the measurement. SUCCESS if successful, FAIL if not, ALL_DONE if all measurements are complete and this measurer is ready to be destructed.
		 */
		virtual MeasurerStatus measure(int step, const std::complex<double> * psi, const double * v, double t) = 0;
		
		/**
		 * Gets the index of the measurer.
		 * @return The index.
		 */
		int getIndex(){return index;};

		/**
		 * Checks if the measurer requires that the density be calculated.
		 * @return True if the density is needed, false if not.
		 */
		bool needsDensity(){return needsDens;};
	};

	/// Records a constant value to file.
	class DoubleConst :
		public Measurer
	{
	private:
		double c;
		std::fstream fil;
	public:
		/**
		 * Constructor.
		 * @param c The constant value to record.
		 * @param name The name of the file to write to. The extension '.dat' will be appended.
		 * @param fol The folder to write to.
		 */
		DoubleConst(double c, const char* name, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t){return MeasurerStatus::ALL_DONE;};
	};

	/// Writes text (8 chars required) to file. Output file is head.dat.
	class Header :
		public Measurer
	{
	private:
		std::fstream fil;
		const char* fname = "head";
	public:
		/**
		 * Constructor.
		 * @param title The text to write. 8 characters required.
		 * @param fol The folder to write to.
		 */
		Header(const char* title, const char* fol);
		
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t) { return MeasurerStatus::ALL_DONE; };
	};

	/// Records the number of grid points in a simulation. Output file is nPts.dat.
	class NPts :
		public Measurer {
	private:
		std::fstream fil;
		const char* fname = "nPts";
	public:
		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param fol The folder to write to.
		 */
		NPts(int nPts, const char* fol);
		
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t) { return MeasurerStatus::ALL_DONE; };
	};

	/// Records the number of time steps in simulation. Output file is nSteps.dat.
	class NSteps :
		public Measurer {
	private:
		std::fstream fil;
		
		const char* fname = "nSteps";
		int steps = 0;
		double tmea = -1;
	public:
		/**
		 * Constructor.
		 * @param fol The folder to write to.
		 */
		NSteps(const char* fol);

		~NSteps();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records dx (spatial) spacing. Output file is dx.dat.
	class DX :
		public Measurer {
	private:
		std::fstream fil;
		const char* fname = "dx";
	public:
		/**
		 * Constructor.
		 * @param dx The spacing.
		 * @param fol The folder to write to.
		 */
		DX(double dx, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t) { return MeasurerStatus::ALL_DONE; };
	};

	/// Records dt (temporal) spacing. Output file is dt.dat.
	class DT :
		public Measurer {
	private:
		std::fstream fil;
		const char* fname = "dt";
	public:
		/**
		 * Constructor.
		 * @param dt The spacing.
		 * @param fol The folder to write to.
		 */
		DT(double dt, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t) { return MeasurerStatus::ALL_DONE; };
	};

	/// Records array of x positions. Output file is xs.dat.
	class XS :
		public Measurer {
	private:
		std::fstream fil;
		const char* fname = "xs";
	public:
		/**
		 * Constructor.
		 * @param len The number of positions.
		 * @param xs The array of positions.
		 * @param fol The folder to write to.
		 */
		XS(int len, const double* xs, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t) { return MeasurerStatus::ALL_DONE; };
	};

	/// Records array of time step times. Output file is ts.dat.
	class TS :
		public Measurer {
	private:
		std::fstream fil;
		const char* fname = "ts";
	public:
		/**
		 * Constructor.
		 * @param fol The folder to write to.
		 */
		TS(const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the original potential at beginning of simulation. Output file is v0.dat.
	class OrigPot :
		public Measurer
	{
	private:
		std::fstream fil;
		int n;
		const char* fname = "v0";
	public:
		/**
		 * Constructor.
		 * @param n The number of grid points.
		 * @param fol The folder to write to.
		 */
		OrigPot(int n, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the absolute value squared of the wave function, downsampling in both space and time. Output file is psi2t.dat.
	class Psi2t :
		public Measurer
	{
	private:
		std::fstream fil;
		
		int nPts;
		const int *nElec;
		int nx, nt;
		int numSteps;
		int curIdx;
		double * psi2b;
		double * psi2s;
		double * xs;
		double * ts;
		const char* fname = "psi2t";

		int *measSteps;

		
	public:
		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param nx The number of spatial points to downsample to.
		 * @param nt The number of time points to downsample to.
		 * @param numSteps The number of time steps.
		 * @param x (in) The array of spatial positions.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 */
		Psi2t(int nPts, int nx, int nt, int numSteps, const double * x, const int* nElec, const char* fol);

		~Psi2t();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records expectation value of energy for each wavefunction at each step. Output file is expectE.dat.
	class ExpectE :
		public Measurer {
	private:
		std::fstream fil;
		
		const char* fname = "expectE";
		int nPts;
		const int* nElec;
		double* rho;
		double dx;
		KineticOperators::KineticOperator * const* kin;
	public:
		/**
		 * Constructor.
		 * @param len The number of grid points.
		 * @param dx The spatial spacing.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 * @param kin (in) The kinetic operator.
		 */
		ExpectE(int len, double dx, const int* nElec, const char* fol, KineticOperators::KineticOperator * const* kin);
		~ExpectE();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records expectation value of position for each wavefunction at each step. Output file is expectX.dat.
	class ExpectX :
		public Measurer {
	private:
		std::fstream fil;
		
		const char* fname = "expectX";
		const double* x;
		int nPts;
		const int* nElec;
		double* scratch;
		double dx;
	public:
		/**
		 * Constructor.
		 * @param len The number of grid points.
		 * @param xs (in) The array of spatial positions.
		 * @param dx The spatial spacing.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 */
		ExpectX(int len, const double* xs, double dx, const int* nElec, const char* fol);
		~ExpectX();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records expectation value of momentum (fairly computationally expensive) for each wavefunction at each step. Output file is expectP.dat.
	class ExpectP :
		public Measurer {
	private:
		std::fstream fil;
		
		const char* fname = "expectP";
		int nPts;
		const int* nElec;
		std::complex<double> *scratch1, *scratch2;
		double dx;
	public:
		/**
		 * Constructor.
		 * @param len The number of grid points.
		 * @param dx The spatial spacing.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 */
		ExpectP(int len, double dx, const int* nElec, const char* fol);

		~ExpectP();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records expectation value of acceleration for each wavefunction at each step. Output file is expectA.dat.
	class ExpectA :
		public Measurer {
	private:
		std::fstream fil;
		
		const char* fname = "expectA";
		int nPts;
		const int* nElec;
		double *scratch1, *scratch2;
		double dx;
	public:
		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param dx The spatial spacing.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 */
		ExpectA(int nPts, double dx, const int* nElec, const char* fol);
		~ExpectA();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the total probability remaining in simulation for each wavefunction at each step. Output file is totProb.dat.
	class TotProb :
		public Measurer
	{
	private:
		double dx;
		std::fstream fil;
		double * psi2;
		int nPts;
		const int* nElec;
		
		const char* fname = "totProb";
	public:
		/**
		 * Constructor.
		 * @param n The number of grid points.
		 * @param dx The spatial spacing.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 */
		TotProb(int n, double dx, const int* nElec, const char* fol);
		~TotProb();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the probability current at the virtual detector position (index) for each wavefunction at each step. Output file is {vdNum}jrd.dat.
	class VDProbCurrent :
		public Measurer {
	private:
		double dx;
		std::fstream fil;
		int nPts;
		const int* nElec;
		int vdPos;
		
		const char* fname = "jrd";
	public:
		/**
		 * Constructor.
		 * @param n The number of grid points.
		 * @param dx The spatial spacing.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param vdPos The position (index) of the virtual detector.
		 * @param vdNum The number (index) of the virtual detector. Used for distinguishing between multiple detectors.
		 * @param name A 4-character descriptor for the detector to be saved within the data file.
		 * @param fol The folder to write to.
		 */
		VDProbCurrent(int n, double dx, const int *nElec, int vdPos, int vdNum, const char* name, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the wave function's value at the virtual detector position (index) for each wavefunction at each step. Output file is {vdNum}psird.dat.
	class VDPsi :
		public Measurer {
	private:
		std::fstream fil;
		int nPts;
		const int* nElec;
		int vdPos;
		
		const char* fname = "psird";
	public:
		/**
		 * Constructor.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param vdPos The position (index) of the virtual detector.
		 * @param vdNum The number (index) of the virtual detector. Used for distinguishing between multiple detectors.
		 * @param name A 4-character descriptor for the detector to be saved within the data file.
		 * @param fol The folder to write to.
		 */
		VDPsi(const int* nElec, int vdPos, int vdNum, const char* name, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the potential at the virtual detector position (index) at each step. Output file is {vdNum}vrd.dat.
	class VDPot :
		public Measurer {
	private:
		std::fstream fil;
		int n;
		int vdPos;
		
		int vdNum;
		int curStep = -1;
		const char* fname = "vrd";
	public:
		/**
		 * Constructor.
		 * @param vdPos The position (index) of the virtual detector.
		 * @param vdNum The number (index) of the virtual detector. Used for distinguishing between multiple detectors.
		 * @param name A 4-character descriptor for the detector to be saved within the data file.
		 * @param fol The folder to write to.
		 */
		VDPot(int vdPos, int vdNum, const char* name, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/**
	 * Records the discrete Fourier transform of the wave function at the virtual detector position (index) for each wavefunction at each step.
	 * This is a memory and disk friendly method of recording the flux spectrum. The discrete Fourier transform is updated at each step, and only saved at the end.
	 * Applies a temporal Tukey window (alpha = 0.05) according to the maximum time provided.
	 * Output file is {vdNum}fluxspecvd.dat.
	 */
	class VDFluxSpec :
		public Measurer {
	private:
		std::fstream fil;
		int vdPos, nSamp;
		bool first = true;
		const int* nElec;
		
		int nPts;
		double ct;
		double dw, tmax, tukeyAl=0.05;
		const char* fname = "fluxspecvd";
		std::complex<double>* wfcs0 = nullptr, * wfcs1 = nullptr, *phss, cumPotPhs, *phaseCalcExpMul, *temp;
	public:
		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param vdPos The position (index) of the virtual detector.
		 * @param vdNum The number (index) of the virtual detector. Used for distinguishing between multiple detectors.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param nSamp The number of samples to take (in frequency space) for the Fourier transform.
		 * @param emax The maximum energy to consider in the Fourier transform.
		 * @param tmax The maximum time.
		 * @param name A 4-character descriptor for the detector to be saved within the data file.
		 * @param fol The folder to write to.
		 */
		VDFluxSpec(int nPts, int vdPos, int vdNum, const int* nElec, int nSamp, double emax, double tmax, const char* name, const char* fol);

		~VDFluxSpec();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the entire wave function at a sample time. Output file is {vdNum}psit.dat.
	class PsiT :
		public Measurer {
	private:
		std::fstream fil;
		int nPts;
		double meaT;
		
		const int* nElec;
		const char* fname = "psit";
		bool done = false;
		double curTime=-1;
	public:
		/**
		 * Constructor.
		 * @param n The number of grid points.
		 * @param meaT The time to record the wave function.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param vdNum The number (index) of the virtual detector. Used for distinguishing between multiple detectors.
		 * @param name A 4-character descriptor for the detector to be saved within the data file.
		 * @param fol The folder to write to.
		 */
		PsiT(int n, double meaT, const int* nElec, int vdNum, const char* name, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the entire potential at a sample time. Output file is {vdNum}pott.dat.
	class PotT :
		public Measurer {
	private:
		std::fstream fil;
		int n;
		double meaT;
		
		int vdNum;
		const char* fname = "pott";
		bool done = false;
	public:
		/**
		 * Constructor.
		 * @param n The number of grid points.
		 * @param meaT The time to record the potential.
		 * @param vdNum The number (index) of the virtual detector. Used for distinguishing between multiple detectors.
		 * @param name A 4-character descriptor for the detector to be saved within the data file.
		 * @param fol The folder to write to.
		 */
		PotT(int n, double meaT, int vdNum, const char* name, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the potential, downsampling in both space and time. Output file is vfunct.dat.
	class Vfunct :
		public Measurer
	{
	private:
		std::fstream fil;
		
		int nPts;
		int nx;
		double maxT;
		int nt;
		int curIdx;
		int* measSteps;
		double * vs;
		double * xs;
		double * ts;
		const char* fname = "Vfunct";
	public:
		/**
		 * Constructor.
		 * @param potNum The number of grid points.
		 * @param nPts The number of spatial points.
		 * @param nx The number of spatial points to downsample to.
		 * @param nt The number of time points to downsample to.
		 * @param numSteps The number of time steps.
		 * @param maxT The maximum time.
		 * @param x (in) The array of spatial positions.
		 * @param fol The folder to write to.
		 */
		Vfunct(int potNum, int nPts, int nx, int nt, int numSteps, double maxT, const double * x, const char* fol);

		~Vfunct();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the total number of wavefunctions (electrons) in the simulation. Output file is nElec.dat.
	class NElec :
		public Measurer {
	private:
		std::fstream fil;
		
		bool first = true;
		const int* nElec;
		const char* fname = "nElec";
		const char* fol;
	public:
		/**
		 * Constructor.
		 * @param nElec (in) The pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 */
		NElec(int* nElec, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Records the expectation value of the energy at the beginning of the simulation. Output file is expectE0.dat.
	class ExpectE0 :
		public Measurer {
	private:
		std::fstream fil;
		
		const char* fname = "expectE0";
		int nPts;
		const int* nElec;
		double dx;
		double tmea;
		double* rho;
		bool first = true;
		KineticOperators::KineticOperator * const* kin;
	public:
		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param dx The spatial spacing.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param fol The folder to write to.
		 * @param kin (in) The kinetic operator.
		 */
		ExpectE0(int nPts, double dx, const int* nElec, const char* fol, KineticOperators::KineticOperator * const* kin);
		~ExpectE0();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};


	/// Records the weights of the wave functions used to calculate their 3D density from their probability density. Output file is wghts.dat.
	class WfcRhoWeights :
		public Measurer {
	private:
		std::fstream fil;
		
		const char* fname = "wghts";
		const int *nElec;
		bool first = true;
		double * const * weights;
	public:
		/**
		 * Constructor.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param weights (in) Pointer to the weights of the wave functions.
		 * @param fol The folder to write to.
		 */
		WfcRhoWeights(const int* nElec, double * const * weights, const char* fol);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Uses the GNUPlotter to plot the density during the simulation. Output is a plot window.
	class DensityPlotter :
		public Measurer {
	private:
		plotting::GNUPlotter* plotter=nullptr;
		int nPts, stepsPerPlot;
		const int* nElec;
		WfcToRho::Density *const dens;
		double *const*wght;
		const double *xs;
		double dx, *rho=nullptr;
		bool pause;
	public:
		

		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param nElec (in) Pointer to the number of electrons. Must be determined by the time 'measure' is called.
		 * @param dx The spatial spacing.
		 * @param xs (in) The array of spatial positions.
		 * @param dens (in) The density calculator object.
		 * @param wght (in) Pointer to the weights of the wave functions.
		 * @param stepsPerPlot The number of time steps to wait before updating the plot. Default is 1.
		 * @param pause Whether to pause and wait for user input after each plot. Default is true.
		 */
		DensityPlotter(int nPts, const int* nElec, double dx, const double* xs, WfcToRho::Density *const dens, double * const * wght, int stepsPerPlot=1, bool pause=true);
		
		~DensityPlotter();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};
	
	/// Uses the GNUPlotter to plot the potential during the simulation. Output is a plot window.
	class PotentialPlotter :
		public Measurer {
	private:
		plotting::GNUPlotter* plotter=nullptr;
		int nPts, stepsPerPlot;
		const double *xs;
		bool pause;
	public:
		

		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param xs (in) The array of spatial positions.
		 * @param stepsPerPlot The number of time steps to wait before updating the plot. Default is 1.
		 * @param pause Whether to pause and wait for user input after each plot. Default is true.
		 */
		PotentialPlotter(int nPts, const double* xs, int stepsPerPlot=1, bool pause=true);

		~PotentialPlotter();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Includes a few basic measurements: nPts, nSteps, dx, dt
	class BasicMeasurers :
		public Measurer {
	private:
		std::vector<Measurer*> meas;
	public:
		/**
		 * Constructor.
		 * @param nPts The number of grid points.
		 * @param dx The spatial spacing.
		 * @param dt The temporal spacing.
		 * @param fol The folder to write to.
		 */
		BasicMeasurers(int nPts, double dx, double dt, const char* fol);

		~BasicMeasurers();
		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};

	/// Manages multiple measurers. Intended for use by the SimulationManager.
	class MeasurementManager :
		public Measurer {
	private:
		std::fstream fil;
		int index = INT_MAX;
		std::vector<Measurer*> meas;
		const char* fname;
	public:
		

		/**
		 * Constructor.
		 * @param fname The folder to write to.
		 */
		MeasurementManager(const char* fname);

		~MeasurementManager();

		/**
		 * Adds a measurer to the manager.
		 * @param m The measurer to add.
		 */
		void addMeasurer(Measurer * m);

		MeasurerStatus measure(int step, const std::complex<double> * psi, const double* v, double t);
	};
}