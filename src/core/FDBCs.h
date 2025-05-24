/**
 * @file FDBCs.h
 * @brief Finite-Difference Boundary Conditions (FDBCs) for the time-dependent Schrödinger equation
 */
#pragma once
#include <complex>
#include "PhysCon.h"
#include "blas.h"

/**
 * @namespace FDBCs
 * @brief Contains classes for finite-difference boundary conditions used in the time-dependent Schrödinger equation (TDSE) solution.
 */
namespace FDBCs
{
	/**
	 * @brief Cyclic array class
	 * @details This class is a cyclic array. It is used to store and continually update the boundary wavefunction history.
	 * It maintains a rotating counter to minimize the number of memory operations.
	 */
    template <typename T>
    class CyclicArray
    {
		template <typename U>
		friend class CyclicArray; // allow access to other CyclicArray instances' private members when template types differ
    private:
        size_t size, idx;
        T* arr;
		/**
		 * Returns the local index of the i-th element in the CyclicArray.
		 * @param i The index to access, can be negative or larger than size.
		 * @return The local index in the cyclic array, wrapping around if necessary.
		 */ 
		template <typename U>
		size_t localIdx(U i) {
			static_assert(std::is_integral<U>::value, "Template parameter U must be an integral type.");
			return (idx + i) % size;
		};
    public:
		/**
		 * Construct a CyclicArray object
		 * @param size The size of the array
		 */
        CyclicArray(size_t size) : size(size), idx(size-1), arr(static_cast<T *>(sq_malloc(sizeof(T) * size))) {};
		
		/**
		 * Construct a CyclicArray object and initialize it with a value
		 * @param size The size of the array
		 * @param val The value to initialize the array with
		 */
		CyclicArray(size_t size, T val) : CyclicArray(size) { std::fill_n(arr, size, val); };

        ~CyclicArray() { sq_free(arr); };

		/**
		 * Push a value into the array at the present index, pushing elements to the left.
		 * @param val The value to push
		 */
        void push(T val) { step(); arr[idx] = val; };

		/**
		 * Push a value into the array at the present index, pushing elements to the right.
		 * @param val The value to push
		 */
		void pushBack(T val) { stepBack(); arr[idx] = val; };

		/**
		 * Set the i-th element in the array
		 * @param i The index of the element to set
		 * @param val The value to set the element to
		 */
		template <typename U>
        void set(U i, T val) { 
			static_assert(std::is_integral<U>::value, "Template parameter U must be an integral type.");
			arr[localIdx(i)] = val; 
		};

		/**
		 * Step the index forward
		 */
        void step() { idx = (idx + 1) % size; };

		/**
		 * Step the index backward
		 */
		void stepBack() { idx = (idx - 1 + size) % size; };

		/**
		 * Get the i-th element in the array
		 * @param i The index of the element to get
		 * @return The value of the i-th element
		 */
		template <typename U>
        T get(U i) { 
			static_assert(std::is_integral<U>::value, "Template parameter U must be an integral type.");
			return arr[localIdx(i)]; 
		};

		/**
		 * Print the array
		 * @details Prints the array in order, starting from the present index.
		 */
        void print() { for (size_t i = 0; i < size; i++) std::cout << arr[i] << " "; std::cout << std::endl; };

		/**
		 * Multiply the array by a scalar
		 * @param val The scalar to multiply by
		 */
		void mul(T val) { for (size_t i = 0; i < size; i++) arr[i] *= val; };

		/**
		 * Perform the inner product with another CyclicArray
		 * @param arr The array to perform the inner product with
		 * @return The inner product
		 * @throws std::invalid_argument If the arrays are not the same size
		 */
		template <typename U>
		decltype(std::declval<T&>()* std::declval<U&>()) inner(CyclicArray<U>* arr) {
			if (this->size != arr->size)
				throw std::invalid_argument("Arrays must be of the same size.");

			if (this->idx == arr->idx) // aligned arrays lead to faster code
				return this->aligned_inner(arr);
			else{
				decltype(std::declval<T&>()* std::declval<U&>()) sum = 0;
				for (size_t i = 0; i < this->size; i++)
					sum += this->get(i) * arr->get(i);
				return sum;
			}
		}

		/**
		 * Perform the inner product with another CyclicArray, assuming the arrays are aligned (at the same index)
		 * This increases performance slightly.
		 * @param arr The array to perform the inner product with
		 * @return The inner product
		 * @throws std::invalid_argument If the arrays are not the same size
		 */
		template <typename U>
		decltype(std::declval<T&>()* std::declval<U&>()) aligned_inner(CyclicArray<U>* arr) {
			if (this->size != arr->size)
				throw std::invalid_argument("Arrays must be of the same size.");

			decltype(std::declval<T&>()* std::declval<U&>()) sum = 0;
			if (this->idx == arr->idx) // aligned arrays lead to faster code
				for (size_t i = 0; i < this->size; i++)
					sum += this->arr[i] * arr->arr[i];
			else
				throw std::invalid_argument("Arrays must be aligned.");
			return sum;
		}

		/**
		 * Perform the inner product with an array.
		 * Uses two call to cblas instead of an explicit loop, making for much faster code.
		 * @param arr The array to perform the inner product with
		 * @return The inner product
		 */
		template <typename U>
		decltype(std::declval<T&>()* std::declval<U&>()) inner(U* arr) {
			decltype(std::declval<T&>()* std::declval<U&>()) sum = 0;
			if constexpr (std::is_same_v<U, T>) //types are the same
				if constexpr (std::is_same_v<T, std::complex<double>>){ // complex double inner product
					std::complex<double> temp;
					cblas_zdotu_sub(this->size - this->idx, this->arr + idx, 1, arr, 1, &temp);
					sum += temp;
					cblas_zdotu_sub(this->idx, this->arr, 1, arr + this->size - this->idx, 1, &temp);
					sum += temp;
				}
			else
				for (size_t i = 0; i < this->size; i++)
					sum += this->get(i) * arr[i];
			return sum;
		}
    };

	/// The side of the physical system that the boundary condition is applied to.
	enum class BCSide {
		LEFT, RIGHT
	};

	/**
	 * @brief Boundary condition class
	 * @details This class is the base class for all boundary conditions. It defines the interface for all boundary conditions.
	 */
	class BoundaryCondition
	{
	public:
		/// Get the diagonal element of the LHS matrix at the boundary position
		virtual std::complex<double> getLHSEle() = 0;
		
		/// Get the element of the LHS matrix adjacent to the boundary position
		virtual std::complex<double> getLHSAdjEle() = 0;

		/** 
		 * Get the RHS value for the boundary condition
		 * @param psibd (in) The boundary wavefunction, \a nElec elements
		 * @param psiad (in) The wavefunction adjacent to the boundary, \a nElec elements
		 * @param vb The potential at the boundary
		 * @param res (out) The RHS value for the boundary condition, \a nElec elements
		 * @param nElec The number of electrons
		 */
		virtual void getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, size_t nElec) = 0; // RHS value for the condition
		
		/**
		 * Finish the present step. This \a must \a be called after the RHS has been calculated and before proceeding with the next time step.
		 * @param psibd (in) The boundary wavefunction, \a nElec elements
		 * @param psiad (in) The wavefunction adjacent to the boundary, \a nElec elements
		 * @param vb The potential at the boundary
		 */
		virtual void finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb) = 0;
		
		/**
		 * Prepare the next step. This \a must \a be called before the RHS is calculated.
		 * @param psibd (in) The boundary wavefunction, \a nElec elements
		 * @param psiad (in) The wavefunction adjacent to the boundary, \a nElec elements
		 * @param vb The potential at the boundary
		 */
		virtual void prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb) = 0;
		
		/** 
		 * Fill the history of the boundary wavefunction assuming the phase change over time provided in \a historialPhaseAdvance
		 * @param psibd (in) The boundary wavefunction, \a nElec elements
		 * @param historialPhaseAdvance (in) The phase advance of the wavefunction at the boundary over each time step (according to the total energy), \a nElec elements
		 * @param vb The potential at the boundary
		 */
		virtual void fillHistory(const std::complex<double>* psibd, const std::complex<double>* historialPhaseAdvance, double vb) = 0;
		
		
		/**
		 * Get the steady-state right-hand-side value for the boundary condition. For finding the open system eigenstate. This is intended for the standard stationary Schrodinger equation.
		 * @param kin The kinetic energy (E-V) at the boundary
		 * @return The steady-state right-hand-side value for the boundary condition
		 */
		virtual std::complex<double> getSteadyRHS(double kin) = 0; // assuming the wavefunction is an eigenstate of the open system, what is the right-hand-side value in the first row?
		
		/**
		 * Get the steady-state left-hand-side diagonal element for the boundary condition. For finding the open system eigenstate. This is intended for the standard stationary Schrodinger equation.
		 * @param kin The kinetic energy (E-V) at the boundary
		 * @return The steady-state left-hand-side diagonal element for the boundary condition
		 */
		virtual std::complex<double> getSteadyLHSEle(double kin) = 0; // '' what is the diagonal first-row LHS component?
		
		/**
		 * Get the steady-state left-hand-side adjacent element for the boundary condition. For finding the open system eigenstate. This is intended for the standard stationary Schrodinger equation.
		 * @param kin The kinetic energy (E-V) at the boundary
		 * @return The steady-state left-hand-side adjacent element for the boundary condition
		 */
		virtual std::complex<double> getSteadyLHSAdjEle(double kin) = 0; // '' what is the first-row LHS component adjacent to the diagonal?
	};

	/// Boundary condition which is the same for all orbitals
	class CommonBC :
		public BoundaryCondition
	{
	public:
		/// @copydoc BoundaryCondition::getLHSEle
		void getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, size_t nElec) {
			for (size_t i = 0; i < nElec; i++)
				res[i] = getRHS(vb);
		}

		/// @copydoc BoundaryCondition::fillHistory
		void fillHistory(const std::complex<double>* psibd, const std::complex<double>* historialPhaseAdvance, double vb) { return; };
		
		/// @copydoc BoundaryCondition::prepareStep
		void prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb) {};
		
		/// @copydoc BoundaryCondition::finishStep
		void finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb) {};

		/// @copydoc BoundaryCondition::getRHS
		virtual std::complex<double> getRHS(double vb) = 0; // RHS value for the condition
	};

	/// Boundary condition which is time-independent (and also the same for all wavefunctions)
	class TimeIndependentBC :
		public CommonBC
	{
	public:
		/// @copydoc BoundaryCondition::getSteadyRHS
		std::complex<double> getSteadyRHS(double kin) {
			return getRHS(kin);
		};

		/// @copydoc BoundaryCondition::getSteadyLHSEle
		std::complex<double> getSteadyLHSEle(double kin) {
			return getLHSEle();
		};

		/// @copydoc BoundaryCondition::getSteadyLHSAdjEle
		std::complex<double> getSteadyLHSAdjEle(double kin) {
			return getLHSAdjEle();
		};
	};

	/// Boundary condition with constant values on the boundaries
	class DirichletBC :
		public TimeIndependentBC
	{
	private:
		std::complex<double> bdVal;
	public:
		/**
		 * Construct a Dirichlet boundary condition with a constant value (for all wavefunctions)
		 * @param bdVal The value to set on the boundary
		 */
		DirichletBC(std::complex<double> bdVal) : bdVal(bdVal) {};

		/// @copydoc TimeIndependentBC::getLHSEle
		std::complex<double> getLHSEle() { return 1.0; };

		/// @copydoc TimeIndependentBC::getLHSAdjEle
		std::complex<double> getLHSAdjEle() { return 0.0; };

		/// @copydoc TimeIndependentBC::getRHS
		std::complex<double> getRHS(double vb) { return bdVal; };
	};

	/// Boundary condition with a constant derivative on the boundary
	class NeumannBC :
		public TimeIndependentBC
	{
	private:
		std::complex<double> bdDer;
		double dx;
		double direction;
	public:
		/**
		 * Construct a Neumann boundary condition with a constant derivative (for all wavefunctions)
		 * @param bdDer The derivative to set on the boundary
		 * @param dx The grid spacing
		 * @param side The side of the boundary (left or right)
		 */
		NeumannBC(std::complex<double> bdDer, double dx, BCSide side) : bdDer(bdDer), dx(dx), direction(side == BCSide::LEFT ? 1 : -1) {};
		
		/// @copydoc TimeIndependentBC::getLHSEle
		std::complex<double> getLHSEle() { return -1.0 / dx * direction; };
		
		/// @copydoc TimeIndependentBC::getLHSAdjEle
		std::complex<double> getLHSAdjEle() { return 1.0 / dx * direction; };
		
		/// @copydoc TimeIndependentBC::getRHS
		std::complex<double> getRHS(double vb) { return bdDer; };
	};

	// TODO: Add effective mass for DTBCs
	/**
	 * @brief Homogeneous Discrete Transparent Boundary Condition (HDTBC)
	 * @details This class implements the homogeneous discrete transparent boundary condition for the Crank-Nicolson method.
	 * It is used to model the open boundary condition for a wavefunction in a finite difference scheme.
	 * The potential at the boundary must be constant.
	 * For more information, see the reference: https://doi.org/10.4310/CMS.2003.v1.n3.a7
	 */
	class UniformHDTransparentBC :
		public BoundaryCondition
	{
	protected:
        size_t order, nElec;
        double dx, dt;
        CyclicArray<std::complex<double>> **psis;
        std::complex<double> kernel0, *kernel;

		bool kernelCalculated = false;
		double kernelVb;

        /**
		 * Calculate the kernel for the boundary condition
		 * @param vb The potential at the boundary
		 * @param dt The time step size. Default is zero, which means the BoundaryCondition's internal value is used (if applicable).
		 */
        void calcKernel(double vb, double dt = 0);
	public:
		/**
		 * Construct a UniformHDTransparentBC object
		 * @param order The order of the boundary condition (length of truncated history to use)
		 * @param nElec The number of electrons (wavefunctions) in the system
		 * @param dx The grid spacing
		 * @param dt The time step
		 */
		UniformHDTransparentBC(size_t order, size_t nElec, double dx, double dt);
        
		~UniformHDTransparentBC();

		/// @copydoc BoundaryCondition::getLHSEle
		std::complex<double> getLHSEle() { 
			if(kernelCalculated) 
				return -kernel0;
			else
				throw std::runtime_error("Kernel not calculated yet. Call prepareStep() first.");
		};

		/// @copydoc BoundaryCondition::getLHSAdjEle
		std::complex<double> getLHSAdjEle() { return 1.0; };
		
		/// @copydoc BoundaryCondition::getRHS
		void getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, size_t nElec);
        
		/// @copydoc BoundaryCondition::finishStep
		void finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb);
		
		/// @copydoc BoundaryCondition::prepareStep
		void prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb);
		
		/// @copydoc BoundaryCondition::getRHS
		void printKernel() { for (size_t i = 0; i < order; i++) std::cout << kernel[i] << " "; std::cout << std::endl; };
		
		/// @copydoc BoundaryCondition::fillHistory
		void fillHistory(const std::complex<double>* psibd, const std::complex<double>* historialPhaseAdvance, double vb);
		
		/// @copydoc BoundaryCondition::getSteadyRHS
		std::complex<double> getSteadyRHS(double kin);

		/// @copydoc BoundaryCondition::getSteadyLHSEle
		std::complex<double> getSteadyLHSEle(double kin);

		/// @copydoc BoundaryCondition::getSteadyLHSAdjEle
		std::complex<double> getSteadyLHSAdjEle(double kin);
	};

	/**
	 * @brief Inhomogeneous Discrete Transparent Boundary Condition (IDTBC)
	 * @details This class implements the inhomogeneous discrete transparent boundary condition for the Crank-Nicolson method.
	 * It is used to model the open boundary condition for a wavefunction in a finite difference scheme.
	 * The inhomogeneity may be used to include an incoming current at the boundary.
	 * For more information, see the reference: https://doi.org/10.4310/CMS.2003.v1.n3.a7
	 */
	class UniformIDTransparentBC :
		public UniformHDTransparentBC
	{
	private:
		std::complex<double> * phaseAdvance, * phs, * adjphs, *ihpsi, *hompsi;
	public:
		/**
		 * Construct a UniformIDTransparentBC object
		 * @param order The order of the boundary condition (length of truncated history to use)
		 * @param nElec The number of electrons (wavefunctions) in the system
		 * @param dx The grid spacing
		 * @param dt The time step
		 * @param k0 (in) The wavevector of the wavefunction(s) at the boundary, \a nElec elements
		 * @param vb The potential at the boundary
		 */
		UniformIDTransparentBC(size_t order, size_t nElec, double dx, double dt, const double* k0, double vb);
		
		~UniformIDTransparentBC() { sq_free(phaseAdvance); sq_free(phs); sq_free(adjphs); sq_free(ihpsi); sq_free(hompsi); };
		
		/// @copydoc BoundaryCondition::prepareStep
		void prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb);
		
		/// @copydoc BoundaryCondition::getRHS
		void getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, size_t nElec);
		
		/// @copydoc BoundaryCondition::finishStep
		void finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb);
		
		/// @copydoc BoundaryCondition::fillHistory
		void fillHistory(const std::complex<double>* psibd, const std::complex<double>* historialPhaseAdvance, double vb);

		/// @copydoc BoundaryCondition::getSteadyRHS
		std::complex<double> getSteadyRHS(double kin);
	};
}