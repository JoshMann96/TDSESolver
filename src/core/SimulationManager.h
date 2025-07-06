/**
 * @file SimulationManager.h
 * @brief Header file for the SimulationManager class.
 */
#pragma once
#include "CORECommonHeader.h"
#include "Measurers.h"
#include "Potentials.h"
#include "Densities.h"
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
	inline operator size_t() const { return val; };
};

/// Tool for tracking progress, calling a callback function wit the current progress (as a double) of the current run call.
class ProgressTracker
{
	private:
		std::function<void(double)> progCallback;
		size_t nSteps = -1, totalCalls = -1, numCalls = 0;
	public:

	/**
	 * Constructor initializes the progress tracker with a callback function and the total number of steps.
	 * #reset must be called to set the total number of steps before use.
	 * @param callback The callback function to be called when progress is made.
	 * @param nCalls The number of calls to the callback function. First call (0) and last call (1.0) are included. Default is 101.
	 */
	ProgressTracker(std::function<void(double)> callback, size_t nCalls = 101) : progCallback(callback) { setNumCallbacks(nCalls); };

	/**
	 * Updates the progress tracker with the current step. The callback function is called when a full percentage of the task is completed.
	 * If multiple percentages are completed in one step, the callback function is called multiple times with each missed percentage.
	 * @param step The current step of the task.
	 */
	void update(size_t step) {
		if(nSteps < 0)
			throw std::runtime_error("ProgressTracker::update: Number of steps not set!");
		// while(step / nSteps > step / totalCalls),  multiplly denominators, use double to avoid overflow
		while ((double)step * totalCalls >= (double)numCalls * nSteps) {
			numCalls++;
			if (progCallback != nullptr)
				progCallback((double) (step) / (double) (nSteps));
		}
	};

	/**
	 * Sets the number of times throughout a whole calculation that the callback function will be called.
	 * @param callback The callback function to be called when progress is made. First call (0) and last call (1.0) are included.
	 */
	void setNumCallbacks(size_t nCalls) { totalCalls = nCalls; };

	/**
	 * Resets the progress tracker with the total number of steps.
	 * @param nSteps The maximum step value for the current run. The last call to update should be equal to this value.
	 */
	void reset(size_t nSteps) {
		this->nSteps = nSteps;
		numCalls = 0;
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
	Densities::Weight* wght = nullptr;
	Densities::Density* dens = nullptr;
	KineticOperators::KineticOperator* kin = nullptr;

	double *ts, *x, dt, dx;
	double **vs, **rhos, **curs, *spatialDamp;
	size_t nPts, nElec;
	bool calcDensityForPot = false;
	bool calcDensityForMeas = false;
	bool calcCurrentForPot = false;
	bool calcCurrentForMeas = false;
	cyclic_int<size_t> index;
	size_t* step;
	std::complex<double> *scratch1, *scratch2;

	bool wavefunctionInitialized = false;
	bool potentialAvailable = false;
	bool weightsCalculated = false;

	Densities::NormalizationScheme normScheme = Densities::UNNORMALIZED;

	/**
	 * Calculate the potential with provided allocated memory.
	 * @param rho The density array to be used for potential calculation. It should be of size nPts. It will contain the density after the operation of calcDensity is true.
	 * @param cur The current density array to be used for potential calculation. It should be of size nPts. It will contain the current density at the current time step.
	 * @param psi The wavefunction array to be used for potential calculation. It should be of size nPts * nElec. It will contain the wavefunction at the current time step.
	 * @param t The current time in the simulation.
	 * @param v The output array to store the calculated potential.
	 * @param virt If true, the potential is calculated for a virtual step (i.e., it does not affect future calls to getV*).
	 * @return Time in microseconds taken to calculate the potential.
	 */
	size_t calculatePotential(double* rho, double* cur, const std::complex<double>* psi, double t, double* v, bool virt);

	/**
	 * Calculate the potential from the raw density array.
	 * @param rho The raw density array to be used for potential calculation. It should be of size nPts.
	 * @param cur The current density array to be used for potential calculation. It should be of size nPts.
	 * @param psi The wavefunction array to be used for potential calculation. It should be of size nPts * nElec.
	 * @param t The current time in the simulation.
	 * @param v The output array to store the calculated potential.
	 * @param virt If true, the potential is calculated for a virtual step (i.e., it does not affect future calls to getV*).
	 * @return Time in microseconds taken to calculate the potential.
	 */
	size_t calculatePotentialFromRawRhoCur(double* rho, double* cur, const std::complex<double>* psi, double t, double* v, bool virt);

	/**
	 * Updates the potential for the given index.
	 * @param idx The index of the potential to be updated. It should be in the range [0, HISTORY_LENGTH).
	 * @param virt If true, the potential is updated for a virtual step.
	 * @return Time in microseconds taken to update the potential.
	 */
	size_t updatePotential(int idx, bool virt);

	/**
	 * Measures the wavefunction and potential at the given index.
	 * @param idx The index of the measurement to be made. It should be in the range [0, HISTORY_LENGTH).
	 * @return Time in microseconds taken to perform the measurement.
	 */
	size_t measure(int idx);

	/**
	 * Helper function for nonlinear Crank-Nicolson SCF iterations for updating the mean potential.
	 * Accounts for the possibile usage of GPU acceleration.
	 * @param kin_cn The Crank-Nicolson kinetic operator to be used for the SCF iterations.
	 * @param meanPot The array to store the mean potential. It should be of size nPts.
	 */
    void updateMeanPotCNNL(KineticOperators::CrankNicolson *kin_cn, double *meanPot);

	std::complex<double> **psis;

	/// Frees the memory contained in the psis multidimensional array.
	void freePsis(){
		for(int i = 0; i < HISTORY_LENGTH; i++){
			if(psis[i]){
				sq_free(psis[i]);
				psis[i] = nullptr;
			}
		}
		wavefunctionInitialized = false;
	}

	double* weights = nullptr;

	/// Calculate the weights. 
	void calcWeights();

	static constexpr int HISTORY_LENGTH = 4;

	ProgressTracker progTracker;

public:

	/**
	 * Constructor initializes the simulation manager with the number of points, spacing, and a callback function for progress tracking.
	 * This implementation of the constructor uses the left boundary position, the grid spacing, and the number of grid points to generate the grid.
	 * @param nPts The number of points in the simulation.
	 * @param xMin The minimum x-coordinate of the simulation.
	 * @param dx The spacing between points in the simulation.
	 * @param dt The time step for the simulation.
	 * @param callback The callback function to be called when progress is made. Default is nullptr.
	 * 					It will be called \a numCallbackCalls times with doubles ranging from 0 to 100.0, inclusive, nullptr for no callback.
	 * @param numCallbackCalls The number of times the callback function will be called. Default is 101.
	 */
	SimulationManager(size_t nPts, double xMin, double dx, double dt, std::function<void(double)> callback = nullptr, size_t numCallbackCalls = 101);

	/**
	 * Constructor initializes the simulation manager with the x-coordinate range, spacing, and a callback function for progress tracking.
	 * This implementation of the constructor uses the left and right boundary positions and the grid spacing to generate the grid. The right boundary position may not actually be included exactly.
	 * @param xMin The minimum x-coordinate of the simulation.
	 * @param xMax The maximum x-coordinate of the simulation.
	 * @param dx The spacing between points in the simulation.
	 * @param dt The time step for the simulation.
	 * @param callback The callback function to be called when progress is made. Default is nullptr.
	 * 					It will be called \a numCallbackCalls times with doubles ranging from 0 to 1.0, inclusive, nullptr for no callback.
	 * @param numCallbackCalls The number of times the callback function will be called. Default is 101.
	 */
	SimulationManager(double xMin, double xMax, double dx, double dt, std::function<void(double)> callback = nullptr, size_t numCallbackCalls = 101) :
		SimulationManager( ( (size_t) ((xMax - xMin) / dx) ) + 1, xMin, dx, dt, callback, numCallbackCalls) {};
		
	/**
	 * Constructor initializes the simulation manager with the x-coordinate range, number of points, and a callback function for progress tracking.
	 * This implementation of the constructor uses the left and right boundary positions and the number of points to generate the grid. The right boundary position may not actually be included exactly.
	 * @param xMin The minimum x-coordinate of the simulation.
	 * @param xMax The maximum x-coordinate of the simulation.
	 * @param nPts The number of points in the simulation.
	 * @param dt The time step for the simulation.
	 * @param callback The callback function to be called when progress is made. Default is nullptr.
	 * 					It will be called \a numCallbackCalls times with doubles ranging from 0 to 1.0, inclusive, nullptr for no callback.
	 * @param numCallbackCalls The number of times the callback function will be called. Default is 101.
	 */
	SimulationManager(double xMin, double xMax, size_t nPts, double dt, std::function<void(double)> callback = nullptr, size_t numCallbackCalls = 101) :
		SimulationManager(nPts, xMin, (double)((xMax - xMin) / (nPts-1)), dt, callback, numCallbackCalls) {};

	virtual ~SimulationManager();

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
	void setWeight(Densities::Weight* nwght) { wght = nwght; }

	/**
	 * Sets the Density calculator to be used in the simulation.
	 * @param ndens The Density calculator to be used in the simulation.
	 */
	void setDensity(Densities::Density* ndens) { dens = ndens; }

	/**
	 * Gets the Weight calculator used in the simulation.
	 * @return The Weight calculator used in the simulation.
	 * @throw std::runtime_error if the weight calculator is not set.
	 */
	Densities::Weight* getWeight () const { 
		if(!wght) 
			throw std::runtime_error("SimulationManager::getWeight: Weight calculator not set!");
		return wght; 
	}

	/**
	 * Gets the Density calculator used in the simulation.
	 * @return The Density calculator used in the simulation.
	 * @throw std::runtime_error if the density calculator is not set.
	 */
	Densities::Density* getDensity () const { 
		if(!dens) 
			throw std::runtime_error("SimulationManager::getDensity: Density calculator not set!");
		return dens; 
	}

	/** 
	 * Returns whether the wavefunction has been initialized.
	 * @return True if the wavefunction has been initialized, false otherwise.
	 */
	bool wavefunctionIsInitialized() const { return wavefunctionInitialized; }

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
	 * If the weights have not been calculated yet, it will try to calculate them first.
	 * @return A pointer to the weights.
	 */
	double* getWeightValues() { 
		if(!weightsCalculated)
			calcWeights();
		return weights;
	}

	/**
	 * Sets the kinetic operator to be used in the simulation.
	 * @param nkin The kinetic operator to be used in the simulation.
	 */
	void setKineticOperator(KineticOperators::KineticOperator* nkin) { kin = nkin; }

	/**
	 * Runs the simulation for \a nSteps iterations.
	 * This function evaluates the type of KineticOperator and Potential and calls the appropriate run function.
	 * @param nSteps The number of steps to run the simulation for.
	 * @param scfIts The number of self-consistent field iterations to perform. Presently only applies to nonlinear Crank-Nicolson calculations. Default is 8.
	 * @param scfTol The tolerance for the self-consistent field iterations. Default is 1e-6.
	 * @throw std::runtime_error if the kinetic operator is not a valid type.
	 * @details If the kinetic operator is the KineticOperators::CrankNicolson method, it will call #runCN_L for linear potentials or #runCN_NL for nonlinear potentials.
	 * If the kinetic operator is a pseudospectral method (KineticOperators::KineticOperator_PSM), it will call #runEPS_U2TU for linear potentials or #runEPS_UW2TUW for nonlinear potentials.
	 * 
	 * For nonlinear potentials with the Crank-Nicolson method, SCF logic is as follows:
	 * If \a scfTol > 0.0 (default behavior) then SCF iterations continue until $\frac{\Delta t}{\hbar}\max_j{|V_j'-V_j|} < scfTol$ or if \a scfIts is reached.
	 * If \a scfTol = 0.0 and scfIts = 0 then no SCF iterations are performed.
	 * If \a scfTol = 0.0 and scfIts != 0 then SCF is performed for \a scfIts iterations.
	 */
	void run(size_t nSteps, size_t scfIts = 8, double scfTol = 1e-6) {
		// is the kinetic operator Crank-Nicolson?
		KineticOperators::CrankNicolson* kin_fdm = dynamic_cast<KineticOperators::CrankNicolson*>(kin);
		if(kin_fdm != nullptr){
			if (canAsyncCalcPot()){
				std::cout << "SimulationManager::run: Using linear Crank-Nicolson." << std::endl;
				runCN_L(nSteps);
			}
			else{
				std::cout << "SimulationManager::run: Using nonlinear Crank-Nicolson." << std::endl;
				runCN_NL(nSteps, scfIts, scfTol);
			}
			return;
		}

		// is the kinetic operator [explicit] pseudospectral?
		KineticOperators::KineticOperator_PSM* kin_ps = dynamic_cast<KineticOperators::KineticOperator_PSM*>(kin);
		if(kin_ps != nullptr){
			if (canAsyncCalcPot()){
				std::cout << "SimulationManager::run: Using linear explicit pseudospectral method." << std::endl;
				runEPS_U2TU(nSteps);
			}
			else{
				std::cout << "SimulationManager::run: Using nonlinear explicit pseudospectral method." << std::endl;
				runEPS_UW2TUW(nSteps);
			}
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
	 * @param scfIts The number of self-consistent field iterations to perform. Default is 8.
	 * @param scfTol The tolerance for the self-consistent field iterations. Default is 1e-6.
	 * 
	 * @details
	 * If \a scfTol > 0.0 (default behavior) then SCF iterations continue until $\frac{\Delta t}{\hbar}\max_j{|V_j'-V_j|} < scfTol$ or if \a scfIts is reached.
	 * If \a scfTol = 0.0 and scfIts = 0 then no SCF iterations are performed.
	 * If \a scfTol = 0.0 and scfIts != 0 then SCF is performed for \a scfIts iterations.
	 */
    void runCN_NL(size_t nSteps, size_t scfIts = 8, double scfTol = 1e-6);

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
	 * Sets the wavefunction to be used in the simulation. If nElec is not set, it will assume there is only 1 state.
	 * @param npsi The wavefunction to be used in the simulation.
	 * @param norm The normalization scheme to be used for the wavefunction. Default is Densities::UNNORMALIZED.
	 */
	void setPsi(const std::complex<double>* npsi, Densities::NormalizationScheme norm = Densities::UNNORMALIZED);

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
			throw std::runtime_error("SimulationManager::getPsi: Wavefunction not initialized!");
		return psis[index];
	};

	/**
	 * Returns a pointer to the density at the present index.
	 * The wavefunction must be initialized and the density calculator must be set.
	 * @return The density, \a npts elements.
	 */
	double* getRho() {
		assert(wavefunctionInitialized);
		if(!calcDensityForPot)
			if(!dens)
				throw std::runtime_error("SimulationManager::getRho: Density calculator not set!");
			else
				dens->calcRho(nPts, nElec, dx, weights, psis[index], rhos[index]);
		return rhos[index];
	};

	/**
	 * Returns a pointer to the current density at the present index.
	 * The wavefunction must be initialized and the kinetic operator must be set.
	 * If the current is not yet calculated, it will be calculated using the wavefunction and weights.
	 * @return The current density, \a npts elements.
	 */
	double* getCur(){
		assert(wavefunctionInitialized);
		if(!calcCurrentForPot)
			if(!dens)
				throw std::runtime_error("SimulationManager::getCur: Density calculator not set!");
			else{
				kin->calcRawCurrent(psis[index], weights, curs[index], nElec);
				dens->applyProfile(nPts, nElec, dx, curs[index]);
			}
		return curs[index];
	};

	/**
	 * Returns a pointer to the potential at the present index.
	 * If the potential is not yet calculated, it will be calculated using the current density and wavefunction. 
	 * The call will be virtual and will not alter the potentials' internal states.
	 * If the wavefunction is not initialized, it will calculate the bare potential.
	 * @return The potential, \a npts elements.
	 */
	double* getV() {
		if(!potentialAvailable){ // potential not yet calculated -- calculate it!
			if(wavefunctionInitialized) // states set but potential not yet calculated
				calculatePotential(rhos[index], curs[index], psis[index], ts[index], vs[index], true);
			else // states not set, calculate initial potential
				pot->getVBare(ts[index], vs[index]);
		}
		return vs[index];
	}

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
	bool canAsyncCalcPot() const { return (!calcDensityForPot && !calcCurrentForPot); }

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
	 * Sets the number of callback calls for the progress tracker.
	 * @param nCalls The number of callback calls to be set.
	 */
	void setNumCallbackCalls(size_t nCalls) { progTracker.setNumCallbacks(nCalls); }

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