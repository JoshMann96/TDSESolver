/**
 * @file KineticOperator.h
 * @brief Kinetic operators which manifest different representations and iteration schemes for the time-dependent Schrodinger equation (TDSE).
 */
#pragma once
#include "CORECommonHeader.h"
#include "MathTools.h"
#include "FDBCs.h"
#include "CuTridiagSolver.h"

/**
 * @namespace KineticOperators
 * @brief Contains classes for different kinetic operators used in the time-dependent Schrodinger equation (TDSE) solution.
 */
namespace KineticOperators {

	/**
	 * @brief Base class for kinetic operators.
	 * @details This class defines the interface for kinetic operators used in solution of the TDSE.
	 * They have the onus of performing time stepping (though individual signatures/implementations vary), evaluating the kinetic energy, and finding eigenstates.
	 */
	class KineticOperator
	{
	public:
		/**
		 * Calculate the kinetic energy for a wavefunction
		 * @param psi The wavefunction for which to calculate the kinetic energy
		 * @return The kinetic energy of the wavefunction
		 */
		virtual double evaluateKineticEnergy(const std::complex<double>* psi) = 0;
		
		/**
		 * Find the eigenstates of the system using this kinetic operator's basis.
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param emin The minimum energy to search for eigenstates
		 * @param emax The maximum energy to search for eigenstates
		 * @param states (out) The eigenstates found. Memory is allocated by the function and must be freed by the user. See implementation notes for element count.
		 * @param allocator (in) A pointer to the allocator to be used.
		 * @param nEigs (out) The number of eigenstates found, allocated by the caller
		 */
		virtual void findEigenStates(const double* v, double emin, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs) = 0;
	};

	/// Base class for kinetic operators that use a PSM (pseudo-spectral method) for time-stepping. The PSM of choice is the split-step Fourier method.
	class KineticOperator_PSM :
		public KineticOperator
	{
	protected:
		uint fftwPlanPolicy;

		/**
		 * Constructor for KineticOperator_PSM -- all it does is set the FFTW plan policy.
		 * @param fftwPlanPolicy The FFTW plan policy to use for FFTW plans (e.g. FFTW_ESTIMATE, FFTW_MEASURE, etc.)
		 */
		KineticOperator_PSM(uint fftwPlanPolicy) : fftwPlanPolicy(fftwPlanPolicy) {};
	public:

		/**
		 * Performs a half-step of the potential followed by a full step of the kinetic operator (half-potential, full-kinetic), returning the result in real space. This is typically run \a before stepOS_UW.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec (in) The number of electrons in the system
		 * @note This function is useful when updating the potential immediately after the kinetic phase for nonlinear systems, increasing the accuracy of the time-stepping.
		 */
		virtual void stepOS_UW2T(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) = 0;

		/**
		 * Performs a half-step of the potential (half-potential), returning the result in real space. This is typically run \a after stepOS_UW2T.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec (in) The number of electrons in the system
		 * @note This function is useful when updating the potential immediately after the kinetic phase for nonlinear systems, increasing the accuracy of the time-stepping.
		 */
		virtual void stepOS_UW(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) = 0;

		/**
		 * Performs a full step of the kinetic operator (half-potential, full-kinetic, half-potential), returning the result in real space. This is fine for linear systems.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec (in) The number of electrons in the system
		 * @note This function is typically used for linear systems, as it performs a full step of the kinetic operator.
		 * @note It is not recommended for nonlinear systems, as the potential cannot be updated after the kinetic propagation step.
		 */
		virtual void stepOS_U2TU(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) = 0;
	};

	/// Base class for kinetic operators that use the PSM (phase-space method) for time-stepping, and allowing for a general dispersion relation.
	class GenDisp_PSM :
		public KineticOperator_PSM
	{
	protected:
		/**
		 * Constructor for GenDisp_PSM.
		 * @param nPts The number of points in the spatial grid
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param fftwPlanPolicy The FFTW plan policy to use for FFTW plans (e.g. FFTW_ESTIMATE, FFTW_MEASURE, etc.)
		 */
		GenDisp_PSM(size_t nPts, double dx, double dt, uint fftwPlanPolicy=FFTW_PATIENT) : nPts(nPts), dx(dx), dt(dt), osKineticEnergy((std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts)), KineticOperator_PSM(fftwPlanPolicy) {};
	public:
		~GenDisp_PSM();
		
		/// @copydoc KineticOperator_PSM::stepOS_UW2T
		void stepOS_UW2T(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec);
		
		/// @copydoc KineticOperator_PSM::stepOS_UW
		void stepOS_UW(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec);

		/// @copydoc KineticOperator_PSM::stepOS_U2TU
		void stepOS_U2TU(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec);

		/// Frees the internally allocated operator matrix (if any) and resets the needMat flag.
		void freeOpMat() {
			if (opMat)
				sq_free(opMat); opMat = nullptr;
			needMat = true;
		}

		/// @copydoc KineticOperator::findEigenStates
		/// \a states will have \a nPts*nPts elements.
		void findEigenStates(const double* v, double emin, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs);

		/// @copydoc KineticOperator::evaluateKineticEnergy
		double evaluateKineticEnergy(const std::complex<double>* psi);

		/**
		 * Defines the dispersion relation on the reciprocal (momentum-space) grid.
		 * @param kinIn (in) The kinetic energy values to set, nPts elements
		 */
		void set_osKineticEnergy(std::complex<double>* kinIn) {
			vtls::copyArray(nPts, kinIn, osKineticEnergy); needMat = true;
		}
	private:
		bool firstStepAll = true, firstStepOne = true, needMat = true;
		//DFTI_DESCRIPTOR_HANDLE dftiHandle = 0, dftiHandleMat = 0, dftiHandleKin = 0;
		fftw_plan fftwAllForward=NULL, fftwAllBackward=NULL, fftwOneForward=NULL, fftwOneBackward=NULL;

		size_t nPts, nElec;
		std::complex<double> *osKineticPhase = nullptr, * osPotentialPhase = nullptr, *opMat = nullptr;
		std::complex<double>* osKineticEnergy = nullptr;
		double dx, dt;

		/// Allocates and calculates the full kinetic part of the Hamiltonian (a dense matrix).
		void calcOpMat();

		/**
		 * Initializes the FFTW plans for the full kinetic operator (all electrons).
		 * @param nElec The number of electrons in the system
		 * @note This function must be called before any time-stepping FFTW plans are executed.
		 */
		void initializeAllFFT(size_t nElec);

		/// Initializes the FFTW plans for a single array (for transforming the potential to reciprocal space).
		void initializeOneFFT();

		/**
		 * Executes a forward Fourier transform on the input.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, \a nPts*nElec elements
		 * @note This function assumes that the input data is in real space and will be transformed to reciprocal space.
		 */
		void executeAllFFTForward(std::complex<double>* targ);

		/**
		 * Executes a backward Fourier transform on the input.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, \a nPts*nElec elements
		 * @note This function assumes that the input data is in reciprocal space and will be transformed to real space.
		 */
		void executeAllFFTBackward(std::complex<double>* targ);

		/**
		 * Executes a forward Fourier transform on the input for a 1D array.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, nPts elements
		 * @note This function assumes that the input data is in real space and will be transformed to reciprocal space.
		 */
		void executeOneFFTForward(std::complex<double>* targ);

		/**
		 * Executes a backward Fourier transform on the input for a 1D array.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, nPts elements
		 * @note This function assumes that the input data is in reciprocal space and will be transformed to real space.
		 */
		void executeOneFFTBackward(std::complex<double>* targ);
		
	};

	/// A free electron (with an effective mass) dispersion relation using the PSM.
	class GenDisp_PSM_FreeElec :
		public GenDisp_PSM
	{
	public:
		/**
		 * Constructor for GenDisp_PSM_FreeElec.
		 * @param nPts The number of points in the spatial grid
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param m_eff The effective mass of the electron (in atomic units, so 1 is the free electron mass)
		 * @param fftwPlanPolicy The FFTW plan policy to use for FFTW plans (e.g. FFTW_ESTIMATE, FFTW_MEASURE, etc.)
		 */
		GenDisp_PSM_FreeElec(size_t nPts, double dx, double dt, double m_eff, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	/** A general polynomial dispersion relation using the PSM.
	 *  Permits general dispersion relations of the form \f$ E(k) = \sum_{i=0}^{nPoly} a_i k^i \f$.
	 */
	class GenDisp_PSM_Series :
		public GenDisp_PSM
	{
	public:
		/**
		 * Constructor for GenDisp_PSM_Series.
		 * @param nPts The number of points in the spatial grid
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param nPoly The number of polynomial coefficients (should be 2 for linear, 3 for quadratic, etc.)
		 * @param polyCoeffs (in) The polynomial coefficients
		 * @param fftwPlanPolicy The FFTW plan policy to use for FFTW plans (e.g. FFTW_ESTIMATE, FFTW_MEASURE, etc.)
		 */
		GenDisp_PSM_Series(size_t nPts, double dx, double dt, size_t nPoly, const double* polyCoeffs, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	/** A general mathematical expression for the dispersion relation using the PSM.
	 *  Permits any mathematical expression to be used as a dispersion relation, e.g. "k^2 + 0.5*k^3".
	 *  @deprecated This is deprecated and will throw a runtime error if used.
	 */
	class GenDisp_PSM_MathExpr :
		public GenDisp_PSM
	{
	public:
		/// @deprecated
		GenDisp_PSM_MathExpr(size_t nPts, double dx, double dt, std::string expr, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	/**
	 * @brief Base class for kinetic operators that use the PSM (phase-space method) for time-stepping, allowing for spatially-dependent dispersion relations.
	 * @details As the operator is no longer diagonal in reciprocal space, the exponent of the operator is approximated as a series.
	 * This is inefficient and somewhat inaccurate. Furthermore, the operation is no longer unitary, so care must be taken when using this operator.
	 * This mixed description includes several dispersion relations and a mask for each. The sum of all the masks should ideally be 1 across the system.
	 * \f$ K(x,k) = \sum_{d=0}^{n_{disp}-1} M_d(x) \cdot K_d(k) \f$
	 */
	class NonUnifGenDisp_PSM :
		public KineticOperator_PSM
	{
	protected:
		/**
		 * Constructor for NonUnifGenDisp_PSM.
		 * @param nPts The number of points in the spatial grid
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param nDisp The number of dispersion relations (e.g. for different regions of the system)
		 * @param expOrder The order of the expansion for the operator (higher is more accurate but less efficient)
		 * @param forceNormalization Whether to maintain the norm of the wavefunction after each step (1 for yes, 0 for no)
		 * @param fftwPlanPolicy The FFTW plan policy to use for FFTW plans (e.g. FFTW_ESTIMATE, FFTW_MEASURE, etc.)
		 */
		NonUnifGenDisp_PSM(size_t nPts, double dx, double dt, size_t nDisp, size_t expOrder, bool forceNormalization, uint fftwPlanPolicy=FFTW_PATIENT) : 
			nPts(nPts), dx(dx), dt(dt), nDisp(nDisp), expOrder(expOrder), forceNorm(forceNormalization), 
			osKineticEnergy((std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nDisp)),
			osKineticMask((double*)sq_malloc(sizeof(double)*nPts*nDisp)),
			KineticOperator_PSM(fftwPlanPolicy) {};
	public:
		~NonUnifGenDisp_PSM();

		/// @copydoc KineticOperator_PSM::stepOS_UW2T
		void stepOS_UW2T(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec);
		
		/// @copydoc KineticOperator_PSM::stepOS_UW
		void stepOS_UW(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec);

		/// @copydoc KineticOperator_PSM::stepOS_U2TU
		void stepOS_U2TU(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec);

		/// Frees the internally allocated operator matrix (if any) and resets the needMat flag.
		void freeOpMat() {
			if (opMat)
				sq_free(opMat); opMat = NULL;
			needMat = true;
		}

		/// @copydoc KineticOperator::findEigenStates
		/// \a states will have \a nPts*nPts elements.
		void findEigenStates(const double* v, double emin, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs);

		/// @copydoc KineticOperator::evaluateKineticEnergy
		double evaluateKineticEnergy(const std::complex<double>* psi);

		/**
		 * Defines the dispersion relation on the reciprocal (momentum-space) grid.
		 * Takes in multiple dispersion relations for each region, and with a mask to define how much of that dispersion relation applies to a gridpoint:
		 * \f$ K(x,k) = \sum_{d=0}^{n_{disp}-1} M_d(x) \cdot K_d(k) \f$
		 * @param kinIn (in) The kinetic energy values to set, nPts*nDisp elements
		 * @param maskIn (in) The mask for the kinetic energy, nPts*nDisp elements. Ideally, the sum over the dispersion axis should be 1.
		 */
		void set_osKineticEnergy(const std::complex<double>* kinIn, const double* maskIn) {
			vtls::copyArray(nPts * nDisp, kinIn, osKineticEnergy); needMat = true; firstStepOne = true;
			vtls::copyArray(nPts * nDisp, maskIn, osKineticMask);
			//take square root, as is required for this method
			for (size_t i = 0; i < nPts * nDisp; i++)
				osKineticEnergy[i] = std::sqrt(osKineticEnergy[i]);
		}
	private:
		bool firstStepAll = true, firstStepOne = true, needMat = true;
		//DFTI_DESCRIPTOR_HANDLE dftiHandle = 0, dftiHandleMat = 0, dftiHandleKin = 0;
		fftw_plan fftwAllForward=NULL, fftwAllBackward=NULL, fftwOneForward=NULL, fftwOneBackward=NULL;

		size_t nPts, nElec, nDisp, expOrder;
		bool forceNorm;
		std::complex<double>* osPotentialPhase = nullptr, * opMat = nullptr;
		std::complex<double>* osKineticEnergy = nullptr;
		std::complex<double>* tempPsi = nullptr, *tempPsiCum = nullptr;
		double* osKineticMask = nullptr, *norms = nullptr;
		double dx, dt;

		/// Allocates and calculates the full kinetic part of the Hamiltonian (a dense matrix).
		void calcOpMat();

		/**
		 * Initializes the FFTW plans for the full kinetic operator (all electrons).
		 * @param nElec The number of electrons in the system
		 * @note This function must be called before any time-stepping FFTW plans are executed.
		 */
		void initializeAllFFT(size_t nElec);

		/// Initializes the FFTW plans for a single array (for transforming the potential to reciprocal space).
		void initializeOneFFT();

		/**
		 * Executes a forward Fourier transform on the input.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, nPts*nElec elements
		 * @note This function assumes that the input data is in real space and will be transformed to reciprocal space.
		 */
		void executeAllFFTForward(std::complex<double>* targ);

		/**
		 * Executes a backward Fourier transform on the input.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, nPts*nElec elements
		 * @note This function assumes that the input data is in reciprocal space and will be transformed to real space.
		 */
		void executeAllFFTBackward(std::complex<double>* targ);

		/**
		 * Executes a forward Fourier transform on the input for a 1D array.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, nPts elements
		 * @note This function assumes that the input data is in real space and will be transformed to reciprocal space.
		 */
		void executeOneFFTForward(std::complex<double>* targ);

		/**
		 * Executes a backward Fourier transform on the input for a 1D array.
		 * @param targ (in/out) The target array to transform, which will be overwritten with the transformed data, nPts elements
		 * @note This function assumes that the input data is in reciprocal space and will be transformed to real space.
		 */
		void executeOneFFTBackward(std::complex<double>* targ);
	};

	/// A kinetic operator that uses the PSM (phase-space method) for time-stepping, allowing for a spatially-dependent effective mass.
	class NonUnifGenDisp_PSM_EffMassBoundary :
		public NonUnifGenDisp_PSM
	{
	public:
		//meff_r and meff_l are relative effective masses (1 for electron rest mass)
		/**
		 * Constructor for NonUnifGenDisp_PSM_EffMassBoundary.
		 * @param nPts The number of points in the spatial grid
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param expOrder The order of the expansion for the operator (higher is more accurate but less efficient)
		 * @param forceNormalization Whether to maintain the norm of the wavefunction after each step (1 for yes, 0 for no)
		 * @param meff_l The effective mass on the left side of the boundary (in atomic units, so 1 is the free electron mass)
		 * @param meff_r The effective mass on the right side of the boundary (in atomic units, so 1 is the free electron mass)
		 * @param transRate The transition rate (decay constant, 1/(length)) across the effective mass boundary. Larger values mark a sharper transition and less accuracy.
		 * @param transPos The position of the transition (in grid points, 0 to nPts-1), where the dispersion relation is a superposition of the two effective masses.
		 * @param edgeRate The rate at which the effective mass transitions near the boundaries of the calculation.
		 * @param fftwPlanPolicy The FFTW plan policy to use for FFTW plans (e.g. FFTW_ESTIMATE, FFTW_MEASURE, etc.)
		 */
		NonUnifGenDisp_PSM_EffMassBoundary(size_t nPts, double dx, double dt, size_t expOrder, bool forceNormalization, double meff_l, double meff_r, double transRate, size_t transPos, double edgeRate, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	/// A kinetic operator that uses the PSM (phase-space method) for time-stepping, allowing for a mathematical expression for the dispersion relation.
	/// @deprecated This class is deprecated and will throw a runtime error if used.
	class NonUnifGenDisp_PSM_MathExprBoundary :
		public NonUnifGenDisp_PSM
	{
	public:
		/// @deprecated
		NonUnifGenDisp_PSM_MathExprBoundary(size_t nPts, double dx, double dt, size_t expOrder, bool forceNormalization, size_t nDisp, std::vector<std::string> exprs, double* transRates, size_t* transPoss, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	/**
	 * @brief Base class for kinetic operators that use a finite-difference method (FDM) for time-stepping.
	 * @details This class defines the interface for kinetic operators that use a finite-difference method (FDM) for time-stepping.
	 * The FDM of choice is the Crank-Nicolson method, which is unconditionally stable and second-order accurate in both time and space.
	 * The possible boundary conditions are defined in the FDBCs namespace, and can be set for the left and right boundaries of the system.
	 */
	class KineticOperator_FDM :
		public KineticOperator
	{
	protected:
		FDBCs::BoundaryCondition *lbc, *rbc;
		size_t nPts;
	public:
		/**
		 * Constructor for KineticOperator_FDM.
		 * @param nPts The number of points in the spatial grid
		 * @param lbc The left boundary condition to use
		 * @param rbc The right boundary condition to use
		 */
		KineticOperator_FDM(size_t nPts, FDBCs::BoundaryCondition* lbc, FDBCs::BoundaryCondition* rbc) : nPts(nPts), lbc(lbc), rbc(rbc) {};
		
		/**
		 * Performs a time step without finalizing the boundary conditions for this step.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec The number of electrons in the system
		 * @note This function is useful for systems where the potential is nonlinear and an approximation of the wavefunction at the next step is desired.
		 * @note If using the GPU, this function will also neither gather the full wavefunction nor override the present state on the GPU.
		 */
		virtual void stepVirtual(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) = 0; // timestep without iterating BCs
		
		/**
		 * Performs a time step and finalizes the boundary conditions for this step.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec The number of electrons in the system
		 */
		virtual void step(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec) = 0;
		
		/**
		 * Finds the eigenstates of the system using this kinetic operator's basis when one of the boundary conditions are inhomogeneous.
		 * This thereby finds the eigenstates of the open system.
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param es (in) The energy values to use for the calculation, \a nPts elements
		 * @param states (out) The eigenstates found, allocated by the caller, must be of size \a nPts*nElec
		 * @param nElec The number of electrons in the system
		 * @throw std::runtime_error if both boundary conditions are not inhomogeneous.
		 */
		virtual void findInhomogeneousEigenStates(const double* v, const double* es, std::complex<double>* states, size_t nElec) = 0;

		/**
		 * For when the system has not been time-integrated, this function calls the underlying boundary conditions' projectHistory function.
		 * @param psi (in) The wavefunction to project, \a nPts*nElec elements (only the left and right boundaries are referenced)
		 * @param phsL (in) The phase advance for the left boundary, \a nElec elements
		 * @param phsR (in) The phase advance for the right boundary, \a nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param nElec The number of electrons in the system
		 */
		void projectHistory(const std::complex<double>* psi, const std::complex<double>* phsL, const std::complex<double>* phsR, const double* v, size_t nElec);

		/**
		 * Calculates the (raw, unprojected) density of the system using the weights provided.
		 * This uses the state on the GPU.
		 * @param weights (in) The weights to use for the calculation, \a nElec elements
		 * @param rho (out) The density calculated, \a nPts elements
		 * @param virt (in) Whether to calculate the density from the virtual state or the regular state on the GPU.
		 * @return true if the calculation was successful, false otherwise (e.g. if the GPU is not being used, this will return false and one must use the CPU instead).
		 */
		virtual bool calcRawRhoByDevice(const double* weights, double* rho, bool virt) = 0;

		/**
		 * Sets the left and right boundary conditions for the system.
		 * @param bc The boundary condition to set
		 * @param side The side of the system to set the boundary condition for (left or right)
		 * @note This function allows for the boundary conditions to be changed dynamically during the simulation.
		 */
		void setBC(FDBCs::BoundaryCondition* bc, FDBCs::BCSide side) { 
			switch (side) {
				case FDBCs::BCSide::LEFT:
					lbc = bc;
					break;
				case FDBCs::BCSide::RIGHT:
					rbc = bc;
					break;
			}
		};
	};

	/**
	 * @brief A kinetic operator that uses the Crank-Nicolson method for time-stepping.
	 * @details This class implements the Crank-Nicolson method for time-stepping, which is unconditionally stable and second-order accurate in both time and space.
	 * The possible boundary conditions are defined in the FDBCs namespace, and can be set for the left and right boundaries of the system.
	 */
	class CrankNicolson:
		public KineticOperator_FDM
	{
		bool useCuda;
		double dx, dt, m_eff;
		std::complex<double> lhsOffDiag0, lhsDiag0, rhsDiag0, rhsOffDiag, potmul; // elements of LHS tridiagonal matrix
		std::complex<double> *d, *ud, *ld, *r_d=nullptr;

		std::complex<double> *rbct=nullptr, *lbct=nullptr, *bct1=nullptr, *bct2=nullptr;
	
		cudaTridiagonalSolverSystem *cuSolver = nullptr;

		/**
		 * Performs a single time step of the Crank-Nicolson method.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec The number of electrons in the system
		 * @param isVirtual (in) Whether to perform a virtual time step (without finalizing the boundary conditions)
		 */
		void _step(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec, bool isVirtual);
	public:
		// note: if useCuda is true, then the CUDA solver will be used for the tridiagonal system
		//       the present state of the system will be managed internally
		/**
		 * Constructor for CrankNicolson iterator.
		 * @param nPts The number of points in the spatial grid
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param m_eff The effective mass of the electron (in atomic units, so 1 is the free electron mass)
		 * @param leftBC The left boundary condition to use
		 * @param rightBC The right boundary condition to use
		 * @param useCuda Whether to use the CUDA solver for the tridiagonal system
		 * @note If useCuda is true, then the CUDA solver will be used for the tridiagonal system and the present state of the system will be managed internally.
		 * @note If useCuda is false, then the CUDA solver will not be used and the present state of the system will not be managed internally.
		 */
		CrankNicolson(size_t nPts, double dx, double dt, double m_eff, FDBCs::BoundaryCondition* leftBC, FDBCs::BoundaryCondition* rightBC, bool useCuda=true);
		
		~CrankNicolson(){
			sq_free(d);
			sq_free(ud);
			sq_free(ld);
			if(r_d)
				sq_free(r_d);
			
			if(rbct)
				sq_free(rbct);
			if(lbct)
				sq_free(lbct);
			if(bct1)
				sq_free(bct1);
			if(bct2)
				sq_free(bct2);

			if(cuSolver)
				delete cuSolver;
		};

		/**
		 * Performs a time step of the Crank-Nicolson method, finalizing the boundary conditions for this step.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec The number of electrons in the system
		 */
		void step(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec){
			_step(psi0, v, spatialDamp, targ, nElec, false);
		};

		/**
		 * Performs a time step of the Crank-Nicolson method without finalizing the boundary conditions for this step.
		 * @param psi0 (in) The initial wavefunction, \a nPts*nElec elements
		 * @param v (in) The potential to use for the calculation, \a nPts elements
		 * @param spatialDamp (in) The spatial damping to apply (for absorptive BCs), \a nPts elements
		 * @param targ (out) The target wavefunction after the time step, \a nPts*nElec elements
		 * @param nElec The number of electrons in the system
		 */
		void stepVirtual(const std::complex<double>* psi0, const double* v, const double* spatialDamp, std::complex<double>* targ, size_t nElec){
			_step(psi0, v, spatialDamp, targ, nElec, true);
		}

		/**
		 * @copydoc KineticOperator_FDM::findEigenStates
		 * \a states will have \a nPts*(*nElec) elements.
		 * @note This uses the tridiagonal representation of the Hamiltonian associated with the three-point stencil of the laplacian.
		 * Because it only finds eigenstates of the Hamiltonian within the closed system, it is not suitable for open (specifically, inhomogeneous) systems.
		 */
		void findEigenStates(const double* v, double emin, double emax, std::complex<double>** states, void* (*allocator)(size_t), size_t* nEigs);

		/// @copydoc KineticOperator_FDM::findInhomogeneousEigenStates
		void findInhomogeneousEigenStates(const double* v, const double* es, std::complex<double>* states, size_t nElec);

		/// @copydoc KineticOperator::evaluateKineticEnergy
		double evaluateKineticEnergy(const std::complex<double>* psi);

		/// @copydoc KineticOperator_FDM::calcRawRhoByDevice
		bool calcRawRhoByDevice(const double* weights, double* rho, bool virt);

		/**
		 * Calculates the wavenumber associated with a given phase advance per time step according to the Crank-Nicolson dispersion relation.
		 * @param phase The phase advance per time step, as a complex number
		 * @param v The potential to use for the calculation
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param m_eff The effective mass of the electron (in atomic units, so 1 is the free electron mass)
		 * @return The wavenumber associated with the given phase advance per time step
		 */
		static double wavenumberFromPhase(std::complex<double> phase, double v, double dx, double dt, double m_eff){
			dx /= PhysCon::a0;
			dt *= PhysCon::auE_ha/PhysCon::hbar;
			v  /= PhysCon::auE_ha;
			
			double cosine = 1.0 - m_eff*dx*dx*( 2.0/dt*std::tan(std::arg(phase)/2.0) - v ) ;

			if(std::abs(cosine) > 1.0)
				throw std::runtime_error("Crank-Nicolson iteration phase is too large.");
			return 1.0/dx/PhysCon::a0 * std::acos(cosine);
		};

		/**
		 * Calculates the wavenumber associated with a given \a total energy according to the Crank-Nicolson dispersion relation.
		 * @param energy The energy to use for the calculation
		 * @param v The potential to use for the calculation
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param m_eff The effective mass of the electron (in atomic units, so 1 is the free electron mass)
		 * @return The wavenumber associated with the given energy
		 */
		static double wavenumberFromEnergy(double energy, double v, double dx, double dt, double m_eff){
			energy /= PhysCon::auE_ha;
			dx /= PhysCon::a0;
			dt *= PhysCon::auE_ha/PhysCon::hbar;
			v  /= PhysCon::auE_ha;
			
			//std::cout << "energy: " << energy << std::endl;
			//std::cout << "v: " << v << std::endl;

			double cosine = 1.0 - m_eff*dx*dx*( energy - v ) ;

			//std::cout << "cosine: " << cosine << std::endl;

			if(std::abs(cosine) > 1.0)
				throw std::runtime_error("Crank-Nicolson iteration phase is too large.");
			return 1.0/dx/PhysCon::a0 * std::acos(cosine);
		};

		/**
		 * Calculates the wavenumber associated with a given \a total energy according to the Crank-Nicolson dispersion relation.
		 * @param energy The energy to use for the calculation
		 * @param v The potential to use for the calculation
		 * @return The wavenumber associated with the given energy
		 */
		double wavenumberFromEnergy(double energy, double v){
			return wavenumberFromEnergy(energy, v, dx, dt, m_eff);
		};

		/**
		 * Calculates the phase advance associated with a given wavenumber according to the Crank-Nicolson dispersion relation.
		 * @param k0 The wavenumber to use for the calculation
		 * @param v The potential to use for the calculation
		 * @param dx The spatial grid spacing
		 * @param dt The time step size
		 * @param m_eff The effective mass of the electron (in atomic units, so 1 is the free electron mass)
		 * @return The phase advance associated with the given wavenumber
		 */
		static std::complex<double> phaseAdvanceFromWavenumber(double k0, double v, double dx, double dt, double m_eff){
			k0 *= PhysCon::a0;
			dx /= PhysCon::a0;
			dt *= PhysCon::auE_ha/PhysCon::hbar;
			v  /= PhysCon::auE_ha;
			std::complex<double> num = 1.0-0.5*dt*PhysCon::im*( (1.0-std::cos(k0*dx))/dx/dx/m_eff + v );
			return num/std::conj(num);
		};

		/**
		 * Calculates the phase advance associated with a given \a total energy according to the Crank-Nicolson dispersion relation.
		 * @param e The energy to use for the calculation
		 * @param dt The time step size
		 * @return The phase advance associated with the given energy
		 */
		static std::complex<double> phaseAdvanceFromEnergy(double e, double dt){
			std::complex<double> num = 1.0-0.5*dt*PhysCon::im*e/PhysCon::hbar;
			return num/std::conj(num);
		};
	};
}