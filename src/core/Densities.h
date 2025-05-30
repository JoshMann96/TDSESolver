/**
 * @file Densities.h
 * @brief Tools for calculating the electron density from the single-particle wavefunctions.
 */
#pragma once
#include "CORECommonHeader.h"
#include "MathTools.h"
#include "KineticOperator.h"

/**
 * @namespace Densities
 * @brief Contains classes and functions for calculating electron density from wavefunctions.
 */
namespace Densities {

	/// Normalization schemes -- whether the wavefunction is normalized or if the wavefunction's absolute magnitude holds physical significance.
	enum NormalizationScheme{
		/// The wavefunction is unnormalized, typically for inhomogeneous/open systems.
		UNNORMALIZED,
		/// The wavefunction is normalized such that the integral of |psi|^2 over the entire simulation space equals 1.
		NORMALIZED
	};

	/**
	 * Base class for weight calculation.
	 * Instances of Weight are responsible for calculating the weights for each electron based on their energies and the corresponding physical model.
	 * The weights calculated may be used to find the density: \f$ \rho = \sum_{i=1}^{n} w_i |\psi_i|^2 \f$ , where \f$ w_i \f$ are the weights and \f$ |\psi_i|^2 \f$ is the squared magnitude of the wavefunction for electron \a i.
	 */
	class Weight {
	public:
		/**
		 * Calculate the weights for the provided set of energies.
		 * @param nElec Number of electrons.
		 * @param energies (in) Array of energies for each electron, size nElec.
		 * @param weights (out) Array to store the calculated weights for each electron, size nElec.
		 * @param norm Normalization scheme to use (UNNORMALIZED or NORMALIZED).
		 * @details The weights calculated may be used to find the density: \f$ \rho = \sum_{i=1}^{n} w_i |\psi_i|^2 \f$, where \f$ w_i \f$ are the weights and \f$ |\psi_i|^2 \f$ is the squared magnitude of the wavefunction for electron \a i.
		 */
		virtual void calcWeights(size_t nElec, const double* energies, double* weights, NormalizationScheme norm) = 0;
	};

	class UniformWeight :
		public Weight
	{
	private:
		double weight;
	public:
		/**
		 * Constructor for UniformWeight.
		 * @param weight The weight to be assigned to each state.
		 */
		UniformWeight(double weight) : weight(weight) {}

		/**
		 * Sets the weights to be a constant value for all states.
		 * @param nElec Number of electrons.
		 * @param energies (in) Array of energies for each electron, size nElec.
		 * @param weights (out) Array to store the calculated weights for each electron, size nElec.
		 * @param norm Normalization scheme to use (UNNORMALIZED or NORMALIZED).
		 */
		void calcWeights(size_t nElec, const double* energies, double* weights, NormalizationScheme norm) {
			for(size_t i = 0; i < nElec; i++)
				weights[i] = weight;
		}
	};

	/**
	 * Produces weights for a bound Fermi gas, where the wavefunctions are \a initially normalized and the system is homogeneous (aside from absorptive BCs).
	 * @details 
	 * This class calculates the weights for each electron based on the Fermi energy level and the energies of the electrons in a bound system.
	 * The bottom of the well is taken to be the midpoint of the maximum and minimum energies provided, minus half of the Fermi energy.
	 * The weight of each state is then calculated as: 
	 * \f$ w_\nu = \frac{2}{3\pi} \frac{m_e E_f}{\hbar^2} \frac{N_e}{\sum_j E_f-E_j} (E_f-E_\nu) \f$ with \f$N_e\f$ the number of states and energies are relative to the bottom of the well.
	 * For a set of wavefunctions corresponding to the eigenstates of a finite well, this results in a nearly flat-top density.
	 */
	class BoundFermiGas :
		public Weight
	{
	private:
		double ef;
	public:
		/**
		 * Constructor for BoundFermiGas.
		 * @param ef The Fermi energy.
		 */
		BoundFermiGas(double ef) : ef(ef) {}
		/**
		 * @copydoc Weight::calcWeights(size_t, const double*, double*, NormalizationScheme)
		 * @throw std::runtime_error if the normalization scheme is UNNORMALIZED
		 */
		void calcWeights(size_t nElec, const double* energies, double* weights, NormalizationScheme norm);
	};

	/**
	 * Produces weights for a semi-infinite Fermi gas, where the wavefunctions correspond to eigenstates of the open system as defined by inhomogeneous boundary conditions.
	 * @details
	 * This class performs the same calculation as BoundFermiGas, but assumes that the wavefunctions have an incoming component of magnitude 1 according to an inhomogeneous boundary condition.
	 * @warning This class expects the boundary conditions to impose incoming planewaves with a form \f$e^{ikx}\f$. It does not check if the boundary condition actually has unity magnitude.
	 */
	class SemiInfiniteFermiGas :
		public Weight
	{
	private:
		double ef;
	public:
		/**
		 * Constructor for SemiInfiniteFermiGas.
		 * @param ef The Fermi energy.
		 */
		SemiInfiniteFermiGas(double ef) : ef(ef) {}

		/**
		 * @copydoc Weight::calcWeights(size_t, const double*, double*, NormalizationScheme)
		 * @throw std::runtime_error if the normalization scheme is NORMALIZED
		 */
		void calcWeights(size_t nElec, const double* energies, double* weights, NormalizationScheme norm);
	};

	/**
	 * Produces weights for a system based on the density of states (DOS) from a file.
	 * @details
	 * This class reads a file containing the density of states (DOS) data and uses it to calculate the weights for each electron based on their energies.
	 * The DOS is interpolated using a cardinal cubic B-spline to provide smooth weights.
	 * The file is a binary format with contents:
	 *  - (int32)\f$\times 1\f$ : number of samples in the DOS data, \a n
	 *  - (double)\f$\times n\f$ : The energy samples IN ELECTRONVOLTS relative to the Fermi level.
	 *  - (double)\f$\times n\f$ : The corresponding DOS values at those energies in \f$\mathrm{\#/m}^3\mathrm{eV}\f$.
	 */
	class FromDOS :
		public Weight
	{
	private:
		double ef, leff;
		boost::math::interpolators::cardinal_cubic_b_spline<double> dosISpline;
	public:
		/**
		 * Constructor for FromDOS.
		 * @param fl The Fermi level (relative to vacuum) for the model system, typically \f$-W\f$.
		 * @param ef The Fermi energy.
		 * @param leff The effective well size (in the same units as the DOS).
		 * @param fil The filename containing the DOS data in binary format.
		 */
		FromDOS(double fl, double ef, double leff, const char* fil);

		/**
		 * @copydoc Weight::calcWeights(size_t, const double*, double*, NormalizationScheme)
		 * @warning This is intended for use with normalized wavefunctions only. Though, it is possible to use it with unnormalized wavefunctions if you know what you're doing.
		 */
		void calcWeights(size_t nElec, const double* energies, double* weights, NormalizationScheme norm);
	};

	/**
	 * Base class for calculating the electron density from wavefunctions and weights calculated from a child of Weight .
	 */
	class Density {
	private:
		double* psi2_work = nullptr;
	public:
		virtual ~Density() { if (psi2_work) sq_free(psi2_work); }

		/**
		 * Calculate the raw electron density, without any post-processing for geometry considerations.
		 * @param nPts Number of grid points in the spatial domain.
		 * @param nElec Number of electrons (wavefunctions).
		 * @param weights (in) Array of weights for each electron, size nElec.
		 * @param psi (in) Array of wavefunctions, size nPts * nElec.
		 * @param psi2_work (in/out) Workspace for squared magnitudes of wavefunctions, size nPts * nElec.
		 * @param rho (out) Array to store the calculated raw density, size nPts.
		 */
		static void calcRawRho(size_t nPts, size_t nElec, const double* weights, const std::complex<double>* psi, double* psi2_work, double* rho);

		/**
		 * Calculate the electron density from wavefunctions and weights, applying any necessary post-processing (e.g., geometry considerations).
		 * @param nPts Number of grid points in the spatial domain.
		 * @param nElec Number of electrons (wavefunctions).
		 * @param dx The grid spacing in the spatial domain.
		 * @param weights (in) Array of weights for each electron, size nElec.
		 * @param psi (in) Array of wavefunctions, size nPts * nElec.
		 * @param rho (out) Array to store the calculated processed density, size nPts.
		 */
		void calcRho(size_t nPts, size_t nElec, double dx, const double* weights, const std::complex<double>* psi, double* rho) {
			if (!psi2_work) 
				psi2_work = (double*)sq_malloc(sizeof(double) * nPts * nElec);
			calcRawRho(nPts, nElec, weights, psi, psi2_work, rho);
			calcRho(nPts, nElec, dx, rho);
		};

		/**
		 * Takes the raw density and processes it to obtain the final electron density with any further geometric considerations.
		 * @param nPts Number of grid points in the spatial domain.
		 * @param nElec Number of electrons (wavefunctions).
		 * @param dx The grid spacing in the spatial domain.
		 * @param rho (in/out) Array of raw density values, size nPts. This will be modified to contain the processed density.
		 */
		virtual void calcRho(size_t nPts, size_t nElec, double dx, double* rho) = 0;
	};

	/// Performs no post-processing on the raw density, simply returning it as is.
	class DirectDensity :
		public Density
	{
	private:
		bool first = true;
	public:
		void calcRho(size_t nPts, size_t nElec, double dx, double* rho);
	};

	/// A density calculator which models a region of cylindrical geometry such that the density decreases further away from the cylinder.
	class CylindricalDensity :
		public Density
	{
	private:
		double center, radius, minX;
		size_t startIndex, endIndex;
		double* thinning=nullptr;
		bool first = true;
		size_t mynPts = 0;
		/**
		 * Initializes the calculation.
		 * @param nPts The number of grid points in the spatial domain.
		 * @param dx The grid spacing in the spatial domain.
		 */
		void doFirst(size_t nPts, double dx);
	public:
		/**
		 * Constructor for CylindricalDensity.
		 * @param center The center of the cylindrical region along the x-axis.
		 * @param radius The radius of the cylindrical region.
		 * @param minX The minimum x-coordinate of the grid, used to determine the start and end indices for thinning.
		 */
		CylindricalDensity(double center, double radius, double minX);

		~CylindricalDensity();
		void calcRho(size_t nPts, size_t nElec, double dx, double* rho);
	};

	// TODO: Make version which is not periodic
	/**
	 * Smooths the density by convolution against a Gaussian.
	 * Uses FFTs to perform the convolution efficiently. This is therefore only suitable for periodic boundary conditions.
	 */
	class GaussianSmoothedDensityPBC :
		public Density
	{
	private:
		bool first = true;
		double *tempRho=nullptr, sig;
		vtls::MaskConvolver<double>* conv = nullptr;
		Density* baseDens = nullptr;
		size_t mynPts = 0;
	public:
		/**
		 * Constructor for GaussianSmoothedDensityPBC.
		 * @param sig The standard deviation of the Gaussian used for smoothing.
		 */
		GaussianSmoothedDensityPBC(double sig) : sig(sig) {}

		/**
		 * Constructor for GaussianSmoothedDensityPBC.
		 * @param sig The standard deviation of the Gaussian used for smoothing.
		 * @param baseDens The base density calculator to use for the raw density calculation.
		 */
		GaussianSmoothedDensityPBC(double sig, Density* baseDens) : sig(sig), baseDens(baseDens) {}

		~GaussianSmoothedDensityPBC();
		void calcRho(size_t nPts, size_t nElec, double dx, double* rho);
	};

	class SmallKernelConvolver :
		public Density
	{
	private:
		double* mask = nullptr, *temp = nullptr;
		size_t maskLen = 0;
		size_t mynPts = 0;
		Density* baseDens = nullptr;
	public:
		/**
		 * Constructor for SmallKernelConvolver.
		 * Uses a cosine-squared kernel which is centered about the middle of the mask.
		 * @param maskLen The length of the kernel mask.
		 */
		SmallKernelConvolver(size_t maskLen);

		/**
		 * Constructor for SmallKernelConvolver.
		 * Uses a cosine-squared kernel which is centered about the middle of the mask.
		 * @param maskLen The length of the kernel mask.
		 * @param baseDens The base density calculator to use for the raw density calculation.
		 */
		SmallKernelConvolver(size_t maskLen, Density* baseDens);

		~SmallKernelConvolver();

		void calcRho(size_t nPts, size_t nElec, double dx, double* rho);
	};
}