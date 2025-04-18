/**
 * @file SimulationManager.h
 * @brief Header file for the SimulationManager class.
 */
#pragma once
#include "CORECommonHeader.h"
#include "Measurers.h"
#include "Potentials.h"
#include "WfcRhoTools.h"
#include "KineticOperator.h"

/// A cyclic integer class that wraps around a maximum value. Useful for managing the local wavefunction and potential history.
template<typename T>
class cyclic_int
{
	static_assert(std::is_integral<T>::value, "cyclic_int only works with integral types");
protected:
	T val, max;
public:
    /// Default constructor initializes the value to 0 and the maximum to 0.
	cyclic_int() : val(0), max(0) {};

	/**
	 * Constructor initializes the value to 0 and the maximum to the given value.
	 * @param max The maximum value (exclusive) for the cyclic integer.
	 */
	cyclic_int(T max) : val(0), max(max) {};

	/**
	 * Constructor initializes the value to the given value and the maximum to the given value.
	 * @param val The initial value of the cyclic integer.
	 * @param max The maximum value (exclusive) for the cyclic integer.
	 */
	cyclic_int(T val, T max) : val(val), max(max) {};

	inline void increment() { val = (val + 1) % max; };
	inline cyclic_int& operator++() { increment(); return *this; }; //prefix
	inline cyclic_int operator++(int) { cyclic_int c = *this; increment(); return c; }; //postfix
	inline cyclic_int operator+(T n) { cyclic_int c(max); c.val = (val + n) % max; return c; };
	inline cyclic_int operator+(int n) { cyclic_int c(max); c.val = (val + n) % max; return c; };
	inline cyclic_int& operator+=(T n) { val = (val + n) % max; return *this; };
	inline cyclic_int& operator=(T n) { val = n % max; return *this; };
	inline operator int() const { return val; };
	inline operator long() const { return val; };
};

/// Tool for tracking progress, calling a callback function when a full percentage of the task is done.
class ProgressTracker
{
	private:
		std::function<void(int)> progCallback;
		int percDone = 0;
		size_t nSteps;
	public:
	/**
	 * Constructor initializes the progress tracker with a callback function. The callback function must take in an integer (0-100) as an argument.
	 * #reset must be called to set the total number of steps before use.
	 * @param callback The callback function to be called when progress is made.
	 */
	ProgressTracker(std::function<void(int)> callback) : progCallback(callback) {};

	/**
	 * Constructor initializes the progress tracker with a callback function and the total number of steps.
	 * @param callback The callback function to be called when progress is made.
	 * @param nSteps The total number of steps for the task.
	 */
	ProgressTracker(std::function<void(int)> callback, size_t nSteps) : progCallback(callback), nSteps(nSteps), percDone(0) {};

	/**
	 * Updates the progress tracker with the current step. The callback function is called when a full percentage of the task is completed.
	 * If multiple percentages are completed in one step, the callback function is called multiple times with each missed percentage.
	 * @param step The current step of the task.
	 */
	void update(size_t step) {
		if(nSteps < 0)
			throw std::runtime_error("ProgressTracker::update: Number of steps not set!");
		while (step*(long)100 / nSteps > percDone) {
			if (progCallback != nullptr)
				progCallback(percDone);
			percDone++;
		}
	};

	/**
	 * Resets the progress tracker with the total number of steps.
	 * @param nSteps The total number of steps for the task.
	 */
	void reset(size_t nSteps) {
		this->nSteps = nSteps;
		percDone = 0;
	};
};

/** Manages simulation by calling Potentials, Measurers, and KineticOperators (with corresponding numerical methods for time integration) 
 * 		using task parallelism as well as the parallelism used within each object.
 * Stores a brief history of the potential and wavefunctions so that measurements can be done in parallel.
 * For linear systems the potential calculation may be done in parallel as well.
 */
class SimulationManager
{
private:
	Potentials::PotentialManager * pot;
	Measurers::MeasurementManager * meas;
	WfcToRho::Weight* wght = nullptr;
	WfcToRho::Density* dens = nullptr;
	KineticOperators::KineticOperator* kin = nullptr;

	double *ts, *x, dt, dx;
	double **vs, **rhos, *spatialDamp;
	size_t nPts, nElec;
	bool calcDensity = false;
	cyclic_int<size_t> index;
	size_t* step;
	std::complex<double> *scratch1, *scratch2;

	bool wavefunctionInitialized = false;

	WfcToRho::NormalizationScheme normScheme = WfcToRho::UNNORMALIZED;

	/**
	 * Calculate the potential with provided allocated memory.
	 * @param rho The density array to be used for potential calculation. It should be of size nPts. It will contain the density after the operation of calcDensity is true.
	 * @param psi The wavefunction array to be used for potential calculation. It should be of size nPts * nElec. It will contain the wavefunction at the current time step.
	 * @param t The current time in the simulation.
	 * @param v The output array to store the calculated potential.
	 * @return Time in microseconds taken to calculate the potential.
	 */
	size_t calculatePotential(double* rho, const std::complex<double>* psi, double t, double* v);

	/**
	 * Calculate the potential from the raw density array.
	 * @param rho The raw density array to be used for potential calculation. It should be of size nPts.
	 * @param psi The wavefunction array to be used for potential calculation. It should be of size nPts * nElec.
	 * @param t The current time in the simulation.
	 * @param v The output array to store the calculated potential.
	 * @return Time in microseconds taken to calculate the potential.
	 */
	size_t calculatePotentialFromRawRho(double* rho, const std::complex<double>* psi, double t, double* v);

	/**
	 * Updates the potential for the given index.
	 * @param idx The index of the potential to be updated. It should be in the range [0, HISTORY_LENGTH).
	 * @return Time in microseconds taken to update the potential.
	 */
	size_t updatePotential(int idx);

	/**
	 * Measures the wavefunction and potential at the given index.
	 * @param idx The index of the measurement to be made. It should be in the range [0, HISTORY_LENGTH).
	 * @return Time in microseconds taken to perform the measurement.
	 */
	size_t measure(int idx);

	std::complex<double> **psis;

	/// Frees the memory contained in the psis multidimensional array.
	void freePsis(){
		for(int i = 0; i < 4; i++){
			if(psis[i]){
				sq_free(psis[i]);
				psis[i] = nullptr;
			}
		}
	}

	double* weights = nullptr;

	/// Calculate the weights. 
	void calcWeights();

	const int HISTORY_LENGTH = 4; // TODO: See how much history is necessary for each scheme (psi/pot/mea asynchronous vs pot/mea asynchronous vs only mea asynchronous)

	ProgressTracker progTracker;

public:

	/**
	 * Constructor initializes the simulation manager with the number of points, spacing, and a callback function for progress tracking.
	 * @param nPts The number of points in the simulation.
	 * @param xMin The minimum x-coordinate of the simulation.
	 * @param dx The spacing between points in the simulation.
	 * @param dt The time step for the simulation.
	 * @param callback The callback function to be called when progress is made. It must take in an integer (0-100) as an argument. nullptr for no callback.
	 */
	SimulationManager(size_t nPts, double xMin, double dx, double dt, std::function<void(int)> callback = nullptr);

	/**
	 * Constructor initializes the simulation manager with the x-coordinate range, spacing, and a callback function for progress tracking.
	 * @param xMin The minimum x-coordinate of the simulation.
	 * @param xMax The maximum x-coordinate of the simulation.
	 * @param dx The spacing between points in the simulation.
	 * @param dt The time step for the simulation.
	 * @param callback The callback function to be called when progress is made. It must take in an integer (0-100) as an argument. nullptr for no callback.
	 */
	SimulationManager(double xMin, double xMax, double dx, double dt, std::function<void(int)> callback = nullptr) :
		SimulationManager((size_t) ((xMax - xMin) / dx), xMin, dx, dt, callback) {};
		
	/**
	 * Constructor initializes the simulation manager with the x-coordinate range, number of points, and a callback function for progress tracking.
	 * @param xMin The minimum x-coordinate of the simulation.
	 * @param xMax The maximum x-coordinate of the simulation.
	 * @param nPts The number of points in the simulation.
	 * @param dt The time step for the simulation.
	 * @param callback The callback function to be called when progress is made. It must take in an integer (0-100) as an argument. nullptr for no callback.
	 */
	SimulationManager(double xMin, double xMax, size_t nPts, double dt, std::function<void(int)> callback = nullptr) :
		SimulationManager(nPts, xMin, (xMax - xMin) / nPts, dt, callback) {};

	~SimulationManager();

	/**
	 * Adds a Measurer to the simulation.
	 * @param nMeas The measurer to be added to the simulation.
	 */
	void addMeasurer(Measurers::Measurer * nMeas);
	
	/** Adds a Potential to simulation.
	 * @param nPot The potential to be added to the simulation.
	 */
	void addPotential(Potentials::Potential * nPot);


	/** 
	 * Adds (rather, multiplies) spatial absorptive region, to be applied at each time step.
	 * Intended to be used for absorbing boundary conditions.
	 * @param arr The array of spatial damping values to be added to the simulation.
	 */
	void addSpatialDamp(const double* arr);

	/**
	 * Sets the weight calculator to be used in the simulation.
	 * @param nwght The weight calculator to be used in the simulation.
	 */
	void setWeight(WfcToRho::Weight* nwght) { wght = nwght; }

	/**
	 * Sets the Density calculator to be used in the simulation.
	 * @param ndens The Density calculator to be used in the simulation.
	 */
	void setDensity(WfcToRho::Density* ndens) { dens = ndens; }

	/**
	 * Gets the Weight calculator used in the simulation.
	 * @return The Weight calculator used in the simulation.
	 */
	WfcToRho::Weight* getWeight () const { return wght; }

	/**
	 * Gets the Density calculator used in the simulation.
	 * @return The Density calculator used in the simulation.
	 */
	WfcToRho::Density* getDensity () const { return dens; }

	/**
	 * Calculates the energies of the wavefunctions for the requested step.
	 * If the step requested is not in the history, it will throw an error.
	 * @param step The step to be considered.
	 * @param energies (out) The array to store the calculated energies.
	 * @throw std::runtime_error if the step is not found in the history.
	 */
	void calcEnergies(size_t step, double* energies) const;

	/**
	 * Returns a pointer to the weights presently being used.
	 * @return A pointer to the weights.
	 */
	double* getWeightValues() const { return weights; }

	/**
	 * Sets the kinetic operator to be used in the simulation.
	 * @param nkin The kinetic operator to be used in the simulation.
	 */
	void setKineticOperator(KineticOperators::KineticOperator* nkin) { kin = nkin; }

	/**
	 * Runs the simulation for \a nSteps iterations.
	 * This function evaluates the type of KineticOperator and Potential and calls the appropriate run function.
	 * @param nSteps The number of steps to run the simulation for.
	 * @throw std::runtime_error if the kinetic operator is not a valid type.
	 * @details If the kinetic operator is the KineticOperators::CrankNicolson method, it will call #runCN_L for linear potentials or #runCN_NL for nonlinear potentials.
	 * If the kinetic operator is a pseudospectral method (KineticOperators::KineticOperator_PSM), it will call #runEPS_U2TU for linear potentials or #runEPS_UW2TUW for nonlinear potentials.
	 */
	void run(size_t nSteps){
		// is the kinetic operator Crank-Nicolson?
		KineticOperators::CrankNicolson* kin_fdm = dynamic_cast<KineticOperators::CrankNicolson*>(kin);
		if(kin_fdm != nullptr){
			if (canAsyncCalcPot())
				runCN_L(nSteps);
			else
				runCN_NL(nSteps);
			return;
		}

		// is the kinetic operator [explicit] pseudospectral?
		KineticOperators::KineticOperator_PSM* kin_ps = dynamic_cast<KineticOperators::KineticOperator_PSM*>(kin);
		if(kin_ps != nullptr){
			if (canAsyncCalcPot())
				runEPS_U2TU(nSteps);
			else
				runEPS_UW2TUW(nSteps);
			return;
		}
		
		std::cerr << "SimulationManager::run: Kinetic operator is not a valid type!" << std::endl;
		throw std::runtime_error("SimulationManager::run: Kinetic operator is not a valid type!");
	}

	/**
	 * Runs \a nSteps iterations using an explicit pseudospectral method with a \a linear potential.
	 * Half of the potential phase is applied first, then the kinetic phase is applied, and finally the other half of the potential phase is applied.
	 * Both measurements and potential calculations are done with task parallelism if possible.
	 * @param nSteps The number of steps to run the simulation for.
	 */
	void runEPS_U2TU(size_t nSteps);

	/**
	 * Runs \a nSteps iterations using an explicit pseudospectral method with a \a nonlinear potential.
	 * Half of the potential phase is applied first, then the kinetic phase is applied, and finally the other half of the potential phase is applied.
	 * The potential is recalculated after the kinetic phase. Therefore, only measurements are done with task parallelism.
	 * @param nSteps The number of steps to run the simulation for.
	 */
	void runEPS_UW2TUW(size_t nSteps);

	/**
	 * Runs \a nSteps iterations using the Crank-Nicolson method with a \a linear potential.
	 * The potential provided to the solver is the average of the potential at the current and next time step.
	 * Both measurements and potential calculations are done with task parallelism if possible.
	 * @param nSteps The number of steps to run the simulation for.
	 */
	void runCN_L(size_t nSteps); // Assumes linear potential

	/**
	 * Runs \a nSteps iterations using the Crank-Nicolson method with a \a nonlinear potential.
	 * The potential provided to the solver is the average of the potential at the current and next time step.
	 * To estimate the potential at the next time step, the wavefunction is propagated to the next time step using the current potential.
	 * The potential is then recalculated using the new wavefunction, and the average is then taken.
	 * Only measurements are done with task parallelism.
	 * @param nSteps The number of steps to run the simulation for.
	 */
	void runCN_NL(size_t nSteps);

	/**
	 * Finds the eigenstates of the system using the given energy range.
	 * Nonlinear potentials assume a neutral charge distribution -- this function does not find a self-consistent solution.
	 * @param emin The minimum energy of the eigenstates to be found.
	 * @param emax The maximum energy of the eigenstates to be found.
	 * @warning If a pseudospectral method is used, the corresponding Hamiltonian is dense and the resulting calculation takes a lot of memory and time.
	 */
	void findEigenStates(double emin, double emax);

	/**
	 * Finds the eigenstates of the system assuming inhomogeneous boundary conditions are in place.
	 * This is only intended to work for finite difference schemes with supported boundary conditions.
	 * @param nElec The number of electrons in the system.
	 * @param energies (in) The eigenstate energies. The boundary conditions must be consistent with these energies.
	 * @throw std::runtime_error if the kinetic operator is not a finite difference method.
	 * @see KineticOperators::KineticOperator_FDM::findInhomogeneousEigenStates
	 */
	void findInhomogeneousEigenStates(size_t nElec, const double* energies);

	/**
	 * Sets the wavefunction to be used in the simulation. If nElec is not set, it will assume there is only 1 electron.
	 * @param npsi The wavefunction to be used in the simulation.
	 * @param norm The normalization scheme to be used for the wavefunction. Default is WfcToRho::UNNORMALIZED.
	 */
	void setPsi(std::complex<double>* npsi, WfcToRho::NormalizationScheme norm = WfcToRho::UNNORMALIZED);

	/// Iterates the simulation index.
	void iterateIndex();

	/**
	 * Returns the number of grid points in the simulation.
	 * @return The number of grid points in the simulation.
	 */
	size_t getNumPoints() const {return nPts;};

	/**
	 * Returns the index which floors the given x-coordinate.
	 * @param xp The x-coordinate to be found.
	 * @return The index of the x-coordinate in the simulation.
	 */
	size_t findXIdx(double xp){return vtls::findValue(nPts, x, xp);}

	/**
	 * Returns the grid spacing in the simulation.
	 * @return The grid spacing in the simulation.
	 */
	double getDX() const {return dx;};

	/**
	 * Returns the time step size in the simulation.
	 * @return The time step size in the simulation.
	 */
	double getDT() const {return dt;};

	/**
	 * Returns the grid array in the simulation.
	 * @return The grid array in the simulation.
	 */
	double* getX() const {return x;};

	/**
	 * Returns a pointer to the wavefunction at the present index.
	 * @return The wavefunction.
	 */
	std::complex<double>* getPsi() const {
		if(!wavefunctionInitialized)
			return nullptr;
		return psis[index];
	};

	/**
	 * Returns a pointer to the density at the present index.
	 * @return The density.
	 */
	double* getRho() {
		if(!wavefunctionInitialized)
			return nullptr;
		if(!calcDensity)
			dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
		return rhos[index];
	};

	/**
	 * Gets the number of electrons currently in the simulation.
	 * @return The number of electrons in the simulation.
	 */
	size_t getNElec() const {return nElec;};

	/**
	 * Returns a pointer to the number of electrons in the simulation.
	 * @return The number of electrons in the simulation.
	 */
	size_t* getNElecPtr() {return &nElec;};

	/**
	 * Determines if the potential can be calculated asynchronously (if it is linear).
	 * @return True if the potential can be calculated asynchronously, false otherwise.
	 */
	bool canAsyncCalcPot() const { return pot->getComplexity() != Potentials::PotentialComplexity::WAVEFUNCTION_DEPENDENT; }

	/**
	 * Returns a pointer to the weights in the simulation.
	 * @return The weights in the simulation.
	 */
	double*const* getWeightsPtr() { return &weights; }

	/**
	 * Returns a pointer to the KineticOperator being used in the simulation.
	 * @return The KineticOperator.
	 */
	KineticOperators::KineticOperator** getKin() { return &kin; }

	/**
	 * Returns a pointer to the PotentialManager being used in the simulation.
	 * @return The PotentialManager.
	 */
	Potentials::Potential* getPotPointer() const { return pot; }

	/**
	 * Finds the electrical surface of an initial state using first-order perturbation theory.
	 * This function calculates the electrical centroid of the electron density.
	 * @param minPos The minimum position of the surface.
	 * @param maxPos The maximum position of the surface.
	 * @return The index of the electrical surface.
	 * @warning This function is not fully tested yet and is likely not working.
	 */
	size_t findElectricalSurfaceCentroidRule(size_t minPos, size_t maxPos);
};