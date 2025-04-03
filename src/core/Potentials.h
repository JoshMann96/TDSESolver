/**
 * @file Potentials.h
 * @brief Contains classes for potentials and potential management.
 */
#pragma once
#include "CORECommonHeader.h" 
#include "Measurers.h"
#include <cstdarg>

/**
 * @namespace Potentials
 * @brief Contains classes for potentials and potential management.
 */
namespace Potentials {
	/**
	 * @namespace ElectricFieldProfiles
	 * @brief Contains classes for electric field profiles which may then be converted to potentials.
	 */
	namespace ElectricFieldProfiles {
		
		/// Base class for electric field profiles.
		class ElectricFieldProfile {
		protected:
			std::complex<double>* fs = nullptr;
		public:
			/**
			 * Get the electric field profile. The returned array's memory is managed internally.
			 * @return The electric field profile.
			 */
			std::complex<double> * getProfile() const {return fs;};

			/**
			 * Constructor.
			 * @param nPts Number of points in the profile.
			 */
			ElectricFieldProfile(int nPts) : fs((std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nPts)) {};

			/**
			 * Destructor.
			 */
			~ElectricFieldProfile(){if(fs) sq_free(fs);};
		};

		/// Constant field between minX and maxX of strength e.
		class ConstantFieldProfile :
			public ElectricFieldProfile {
		public:
			/**
			 * Constructor.
			 * @param nPts Number of gridpoints.
			 * @param x (in) Array of positions.
			 * @param ef Electric field strength.
			 * @param minX Minimum x.
			 * @param maxX Maximum x.
			 */
			ConstantFieldProfile(int nPts, const double * x, double ef, double minX, double maxX);
		};

		/**
		 * Field which models a geometric field enhancement near a cylinder with a provided enhancement factor.
		 * The field outside the cylinder is of the form \f$(\gamma-1)\frac{R}{x-x_0}+1\f$. 
		 * Furthermore, the field is linearly tapered to zero at the set maxX. The linear tapering begins at a position such that the field is differentiable (except for at maxX).
		 */
		class CylindricalToLinearProfile :
			public ElectricFieldProfile {
		public:
			/**
			 * Constructor.
			 * @param nPts Number of gridpoints.
			 * @param x (in) Array of positions.
			 * @param minX Minimum x to apply the field. This is assumed to be the surface of the cylinder.
			 * @param maxX Maximum x to apply the field.
			 * @param r Radius of the cylinder.
			 * @param eMax Maximum electric field strength (at minX).
			 * @param enhFact Enhancement factor.
			 */
			CylindricalToLinearProfile(int nPts, const double * x, double minX, double maxX, double r, double eMax, double enhFact);
		};

		/**
		 * Field which models a geometric field enhancement near a cylinder with a provided enhancement factor.
		 * The field outside the cylinder is of the form \f$(\gamma-1)\frac{R}{x-x_0}+1\f$.
		 * Furthermore, the field is tapered to zero at the set maxX. The tapering is applied by multiplying by a polynomial smooth step function of order 13 (7th order differentiable)
		 */
		class CylindricalToCutoffProfile :
			public ElectricFieldProfile {
		public:
			/**
			 * Constructor.
			 * @param nPts Number of gridpoints.
			 * @param x (in) Array of positions.
			 * @param minX Minimum x to apply the field. This is assumed to be the surface of the cylinder.
			 * @param maxX Maximum x to apply the field.
			 * @param r Radius of the cylinder.
			 * @param eMax Maximum electric field strength (at minX).
			 * @param enhFact Enhancement factor.
			 * @param decayLength Length over which the field is tapered to zero on the right boundary.
			 */
			CylindricalToCutoffProfile(int nPts, const double * x, double minX, double maxX, double r, double eMax, double enhFact, double decayLength);
		};

		/**
		 * Field which models the field of a surface plasmon polariton \a within a plasmonic material.
		 * This does not produce the field external to the material.
		 */
		class InternalPlasmonicFieldProfile :
			public ElectricFieldProfile {
		public:
			/**
			 * Constructor.
			 * @param nPts Number of gridpoints.
			 * @param x (in) Array of positions.
			 * @param minX Minimum x to apply the field.
			 * @param maxX Maximum x to apply the field. This is assumed to be the material surface.
			 * @param eMax Maximum \a vacuum electric field strength.
			 * @param lam Wavelength of the incident light which drove the surface plasmon.
			 * @param er Relative permittivity of the plasmonic material.
			 * @param cond Conductivity of the plasmonic material. (Not referenced.)
			 */
			InternalPlasmonicFieldProfile(int nPts, const double * x, double minX, double maxX, double eMax, double lam, std::complex<double> er, double cond);
		};

		/**
		 * Field which has an exponential dropoff profile. Effectively models a surface plasmon polariton's \a vacuum fields.
		 * The field is linearly tapered to zero at the set maxX. The linear tapering begins at a position such that the field is differentiable (except for at maxX).
		 */
		class ExponentialToLinearProfile :
			public ElectricFieldProfile {
		public:
			/**
			 * Constructor.
			 * @param nPts Number of gridpoints.
			 * @param x (in) Array of positions.
			 * @param minX Minimum x to apply the field. This is the position where the field is maximum.
			 * @param maxX Maximum x to apply the field.
			 * @param r Decay length of the exponential field. (\f$E~\exp(-\frac{x-x_0}{r})\f$)
			 * @param eMax Maximum electric field strength (at minX).
			 */
			ExponentialToLinearProfile(int nPts, const double* x, double minX, double maxX, double r, double eMax);
		};

		/**
		 * Field which is defined in a binary file.
		 * The binary file should contain the following:
		 *  - (int)\f$\times 1\f$ : Number of samples in the file \a n
		 *  - (double)\f$\times n\f$ : The positions \a x at which the field is sampled
		 *  - (double)\f$\times n\f$ : The real component of the electric field values \a E at those positions
		 *  - (double)\f$\times n\f$ : The imaginary component of the electric field values \a E at those positions
		 * To be consistent with the implementation's \a emax the maximum field in the file should have absolute value 1.
		 */
		class FileFieldProfile :
			public ElectricFieldProfile {
		public:
			/**
			 * Constructor.
			 * @param nPts Number of gridpoints.
			 * @param x (in) Array of positions.
			 * @param offset Positional offset to apply to the field.
			 * @param rightDecayPos Position at which the field starts to decay to zero on the right side.
			 * @param leftDecayPos Position at which the field starts to decay to zero on the left side.
			 * @param decayLength Length over which the field decays to zero on the right and left sides.
			 * @param emax Maximum electric field strength, assuming the peak field in the file is normalized to 1.0.
			 * @param fil File to read from.
			 */
			FileFieldProfile(int nPts, const double * x, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, const char * fil);
		};
	};

	/**
	 * @namespace Envelopes
	 * @brief Contains classes for temporal pulse envelopes.
	 */
	namespace Envelopes {
		/// Template class (no pulse).
		class Envelope {
		public:
			/**
			 * Get the envelope value at time \a t.
			 * @param t Time.
			 * @return The envelope value.
			 */
			virtual double getValue(double t) = 0;
		};

		/// Gaussian envelope.
		class GaussianEnvelope :
			public Envelope {
		private:
			double tau;
			double tmax;
		public:
			/**
			 * Constructor.
			 * @param tau Full-width-half-max power of the Gaussian.
			 * @param tmax Time at which the Gaussian is centered.
			 */
			GaussianEnvelope(double tau, double tmax);

			double getValue(double t);
		};

		/**
		 * Smoothed Gaussian envelope.
		 * This envelope smoothly transitions from zero to a Gaussian envelope using a 13th order smooth polynomial step function.
		 */
		class SmoothedInitialGaussianEnvelope :
			public Envelope {
		private:
			double tau;
			double tmax;
			double buf;
		public:
			/**
			 * Constructor.
			 * @param tau Full-width-half-max power of the Gaussian.
			 * @param tmax Time at which the Gaussian is centered.
			 * @param bufferTime Time over which the envelope transitions from zero to Gaussian.
			 */
			SmoothedInitialGaussianEnvelope(double tau, double tmax, double bufferTime);

			double getValue(double t);
		};

		/// Cosine-squared envelope.
		class CosSquaredEnvelope :
			public Envelope {
		private:
			double tau;
			double tmax;
			double a_t = 2.0 * std::acos(std::pow(0.5, 0.25)); // cos^4(a_t*(tau/2)/tau)=1/2 s.t. tau = FWHM-power
		public:
			/**
			 * Constructor.
			 * @param tau Full-width-half-max power of the cosine-squared envelope.
			 * @param tmax Time at which the envelope is centered.
			 */
			CosSquaredEnvelope(double tau, double tmax);

			double getValue(double t);
		};
	}

	/// Enum for potential complexity.
	enum PotentialComplexity{
		/// potential is a constant, only needs to be evaluated once
		STATIC, 
		/// potential changes with time but is independent of the wavefunction
		DYNAMIC, 
		/// the potential depends on the wavefunction (including density-functional potentials)
		WAVEFUNCTION_DEPENDENT 
	};

	/// Template class for potential.
	class Potential {
	public:
		virtual ~Potential() = default;

		/**
		 * Get the potential energy at time \a t. This is the potential ignoring the wavefunction or density.
		 * @param t Time.
		 * @param targ (out) Array to store the potential.
		 */
		virtual void getVBare(double t, double * targ) = 0;

		/**
		 * Get the potential energy at time \a t.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param t Time.
		 * @param targ (out) Array to store the potential.
		 */
		virtual void getV(const double * rho, const std::complex<double> * psi, double t, double * targ) = 0;

		/**
		 * Get the complexity of the potential.
		 * @return The complexity.
		 */
		virtual PotentialComplexity getComplexity() = 0;
	};

	/**
	 * Read potential energy from binary data file.
	 * The first element of the file should be an integer representing the number of points provided in the potential (n).
	 * The next n doubles are the positions that these are sampled at.
	 * The next n doubles are the (energy) potentials.
	 */
	class FilePotential :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		/**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param offset Positional offset to apply to the potential.
		 * @param fil File to read from.
		 * @param refPoint Reference point for the potential.
		 */
		FilePotential(int nPts, const double * x, double offset, const char * fil, int refPoint);

		~FilePotential();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	/**
	 * Potential which models a bias field.
	 * The field is zero outside of the bias region.
	 * The field is linearly ramped up from zero to the bias field strength over a buffer region.
	 */
	class BiasFieldPotential :
		public Potential
	{
	private:
		double * v;
		double tstart, tbuf;
		int nPts;
	public:
	    /**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param tstart Time at which the bias field starts.
		 * @param tbuf Time over which the bias field ramps up.
		 * @param xmin Minimum x to apply the field.
		 * @param xmax Maximum x to apply the field.
		 * @param xmin_buf Minimum x to start the ramp up.
		 */
		BiasFieldPotential(int nPts, const double * x, double tstart, double tbuf, double xmin, double xmax, double xmin_buf, double xmax_buf, double fieldStrength, int refPoint);
		
		~BiasFieldPotential();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> *  psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	/// Coulomb potential. The depth of the potential is capped by a separation distance of dx.
	class CoulombPotential :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
	    /**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param ne Elementary charges.
		 * @param chargePos Position of the charge.
		 * @param minX Minimum x to apply the field.
		 * @param maxX Maximum x to apply the field.
		 * @param refPoint Reference point for the potential.
		 */
		CoulombPotential(int nPts, const double * x, double ne, double chargePos, double minX, double maxX, int refPoint);
		
		~CoulombPotential();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> *  psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	/**
	 * Potential which models a finite box.
	 */
	class FiniteBox :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		/**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param left Left boundary of the box.
		 * @param right Right boundary of the box.
		 * @param vin Potential inside the box.
		 * @param refPoint Reference point for the potential.
		 */
		FiniteBox(int nPts, const double * x, double left, double right, double vin, int refPoint);
		
		~FiniteBox();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	/// Wacther's Jellium potential.
	class JelliumPotential :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
		/**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param center Center of the potential's surface.
		 * @param ef Fermi energy.
		 * @param w Work function.
		 * @param refPoint Reference point for the potential.
		 */
		JelliumPotential(int nPts, const double * x, double center, double ef, double w, int refPoint);
		
		~JelliumPotential();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	/// Jellium potential with a backing such that it smoothly returns to vacuum level on the left side.
	class JelliumPotentialBacked :
		public Potential
	{
	private:
		double * v;
		int nPts;
	public:
	    /**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param center Center of the potential's surface.
		 * @param ef Fermi energy.
		 * @param w Work function.
		 * @param backStart Position (inner side) at which the backing starts.
		 * @param backWidth Total width of the backing.
		 * @param refPoint Reference point for the potential.
		 */
		JelliumPotentialBacked(int nPts, const double * x, double center, double ef, double w, double backStart, double backWidth, int refPoint);
		
		~JelliumPotentialBacked();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	/// Shielded atomic potential, averaged across an infinite plane parallel to surface.
	class ShieldedAtomicPotential :
		public Potential {
	private:
		double * v;
		int nPts;
	public:
		/**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param center Position of the potential center.
		 * @param latticeSpacing Transverse square lattice spacing.
		 * @param zProtons Number of elementary charges.
		 * @param decayLength Decay length of the shielded potential. Use, e.g., the Thomas-Fermi length.
		 */
		ShieldedAtomicPotential(int nPts, const double * x, double center, double latticeSpacing, double zProtons, double decayLength);
		
		~ShieldedAtomicPotential();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::STATIC;};
	};

	/// Converts an electric field profile and envelope to a potential.
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
	    /**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param fieldProfile Electric field profile.
		 * @param dx Grid spacing.
		 * @param phase Phase of the field relative to the envelope center. The field is cosine-like.
		 * @param tmax Time at which the field is centered.
		 * @param lam Wavelength of the field.
		 * @param env Envelope.
		 * @param refPoint Reference point for the potential.
		 */
		ElectricFieldProfileToPotential(int nPts, ElectricFieldProfiles::ElectricFieldProfile * fieldProfile, double dx, double phase, double tmax, double lam, Envelopes::Envelope * env, int refPoint);
		
		~ElectricFieldProfileToPotential();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::DYNAMIC;};
	};

	/** 
	 * Class of potentials which generally depend on the wavefunction or density and are zero for the initial state -- they only embody the change in potential.
	 */
	class NonlinearDynamicalPotential :
		public Potential
	{
	private:
		int assembled = 0;

		/// Assemble the potential. See implementations for details on inputs.
		virtual void _assemble(const double* rho, const std::complex<double> * psi, va_list args) = 0; // child needs to implement assembly

		/**
		 * Protected method to calculate the potential at time \a t.
		 * This is the actual implementation that computes the potential based on the density and wavefunction. #getV ensures that the potential has been assembled before calling this.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param t Time.
		 * @param targ (out) Array to store the potential.
		 */
		virtual void getV_(const double* rho, const std::complex<double> * psi, double t, double * targ) = 0; // child implements unprotected potential calculation.
	public:
	    /**
		 * Assemble the potential.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param ... Additional arguments. See child classes' #_assemble documentation for specifics.
		 */
		void assemble(const double* rho, const std::complex<double> * psi, ...){
			// collect arguments
			va_list args;
			va_start(args, psi);

			// call the assembly function
			_assemble(rho, psi, args);
			
			va_end(args);

			assembled = 1; 
		};

		/**
		 * Check if the potential has been assembled.
		 * @return 1 if assembled, 0 otherwise.
		 */
		int isAssembled(){return assembled;};

		/**
		 * Get the potential energy at time \a t. This is the potential ignoring the wavefunction or density.
		 * @param t Time.
		 * @param targ (out) Array to store the potential.
		 */
		virtual void getVBare(double t, double * targ) = 0;

		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ){
			if(!isAssembled())
				throw std::runtime_error("NonlinearDynamicalPotential not assembled before evaluation.");
			getV_(rho, psi, t, targ);
		};

		virtual PotentialComplexity getComplexity() = 0;
	};

	/// Tool which integrates the current which passes through a point.
	class CurrentIntegrator
	{
	private:
		double integratedFlux, tPrev, dx;
		double * const * weights;
		int evalPoint, nPts, side;
		const int * nElec;
	public:
		/**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param dx Grid spacing.
		 * @param evalPoint Point at which to evaluate the current.
		 * @param side Side to use to evaluate the first derivative. 0 for central, 1 for right, -1 for left.
		 * @param nElec (in) Number of electrons.
		 * @param weights (in) Weights for the current calculation.
		 */
		CurrentIntegrator(int nPts, double dx, int evalPoint, int side, const int* nElec, double * const * weights);
		
		/**
		 * Integrate the current for the present time step.
		 * The time of last call is recorded internally and the current is multiplied by the difference in time.
		 * @param psi (in) Wavefunction.
		 * @param t Time.
		 */
		void integrate(const std::complex<double>* psi, double t);

		/**
		 * Get the integrated flux.
		 * @return The integrated flux.
		 */
		double getIntegratedFlux() const {return integratedFlux;};
	};

	/**
	 * Potential which models a perfectly conducting cylindrical cathode with external charges as cylindrical sheathes.
	 * The potential effectively includes the collective image charge, with charges closer to the cathode shielding the fields for those farther out.
	 * The current emitted on the right-side is integrated and used to find the total charge remaining within the cathode, regardless of the total charge bounded to the left.
	 */
	class CylindricalImageCharge :
		public NonlinearDynamicalPotential
	{
	private:
		int nPts, refPoint, posMin, posMax, surfPos;
		double dx, ef, w, rad, * origPot, * potTemp, * genTemp, * lrxr, * myRho, *nsMask, *dethin;
		const double *x;
		/**
		 * Calculate the potential at time \a t.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param cur_t Time.
		 * @param targ (out) Array to store the potential.
		 */
		void calcPot(const double* rho, const std::complex<double>* psi, double cur_t, double* targ);
		CurrentIntegrator * curInt;

		/**
		 * Assemble the potential with additional arguments.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param args Additional argument:
		 * - \a sp (int): The index of the surface position.
		 */
		void _assemble(const double* rho, const std::complex<double>* psi, va_list args);

		void getV_(const double* rho, const std::complex<double>* psi, double t, double* targ);
	public:
		/**
		 * Constructor. The surface position must be passed through #assemble. See #_assemble for details.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param dx Grid spacing.
		 * @param ef Electric field strength.
		 * @param w Work function.
		 * @param rad Radius of the cylinder.
		 * @param nElec (in) Number of electrons.
		 * @param weights (in) Weights for the density calculation. (Not referenced.)
		 * @param posMin Minimum x (index) to apply the field. This is assumed to be the surface of the cylinder.
		 * @param posMax Maximum x (index) to apply the field.
		 * @param refPoint Reference point (index) for the potential.
		 */
		CylindricalImageCharge(int nPts, const double* x, double dx, double ef, double w, double rad, const int* nElec, double * const * weights, int posMin, int posMax, int refPoint);
		
		~CylindricalImageCharge();
		void getVBare(double t, double* targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::WAVEFUNCTION_DEPENDENT;};
	};

	/**
	 * A Hartree potential which uses a planar charge geometry to the left of the surface position and a cylindrical charge geometry to the right.
	 */
	class PlanarToCylindricalHartree :
		public NonlinearDynamicalPotential
	{
	private:
		int nPts, refPoint, * nElec, posMin, posMax, surfPos;
		double dx, rad, * origPot, * potTemp, * fieldScaler, *myRho, *dethin;
		double originalCharge;
		/**
		 * Calculate the potential at time \a t.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param t Time.
		 * @param targ (out) Array to store the potential.
		 */
		void calcPot(const double* rho, const std::complex<double>* psi, double t, double* targ);
		CurrentIntegrator * curInt;
		double totalCharge;

		/**
		 * Assemble the potential with additional arguments.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param args Additional argument:
		 * - \a sp (int): The index of the surface position.
		 */
		void _assemble(const double* rho, const std::complex<double>* psi, va_list args);

		void getV_(const double* rho, const std::complex<double>* psi, double t, double* targ);
	public:
		/**
		 * Constructor. The surface position must be passed through #assemble. See #_assemble for details.
		 * @param nPts Number of points.
		 * @param x (in) Array of positions.
		 * @param dx Grid spacing.
		 * @param rad Radius of the cylinder.
		 * @param nElec (in) Number of electrons.
		 * @param weights (in) Weights for the density calculation. (Not referenced.)
		 * @param posMin Minimum x (index) to apply the field.
		 * @param posMax Maximum x (index) to apply the field.
		 * @param refPoint Reference point (index) for the potential.
		 */
		PlanarToCylindricalHartree(int nPts, const double* x, double dx, double rad, const int* nElec, double * const * weights, int posMin, int posMax, int refPoint);
		~PlanarToCylindricalHartree();
		void getVBare(double t, double* targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::WAVEFUNCTION_DEPENDENT;};
	};

	/// Types of local density approximation (LDA) functionals.
	enum class LDAFunctionalType {
		/// Slater exchange functional
		X_SLATER,
		/// PW correlation functional
		C_PW
	};
	
	class LDAFunctional :
		public NonlinearDynamicalPotential
	{
	private:
		int nPts, refPoint, * nElec;
		double * origPot, * rho, dx;
		/**
		 * Calculate the potential.
		 * @param rho (in) Density.
		 * @param targ (out) Array to store the potential.
		 */
		void calcPot(const double* rho, double* targ);

		LDAFunctionalType typ;

		/**
		 * Assemble the potential with additional arguments.
		 * @param rho (in) Density.
		 * @param psi (in) Wavefunction.
		 * @param args Additional arguments (none needed).
		 */
		void _assemble(const double* rho, const std::complex<double>* psi, va_list args);
		void getV_(const double* rho, const std::complex<double>* psi, double t, double* targ);
	public:

		/**
		 * Constructor. #assemble must be called before using the potential. See #_assemble for details.
		 * @param typ Type of LDA functional.
		 * @param nPts Number of points.
		 * @param dx Grid spacing.
		 * @param refPoint Reference point (index) for the potential.
		 */
		LDAFunctional(LDAFunctionalType typ, int nPts, double dx, int refPoint);

		~LDAFunctional();
		void getVBare(double t, double* targ);
		PotentialComplexity getComplexity(){return PotentialComplexity::WAVEFUNCTION_DEPENDENT;};
	};

	/// Combines multiple potentials into a single potential. The complexity of this potential is the most complex of its constituents.
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
		/**
		 * Constructor.
		 * @param nPts Number of points.
		 * @param numSPots Number of static potentials.
		 * @param numDPots Number of dynamic potentials.
		 * @param numWPots Number of wavefunction-dependent potentials.
		 * @param staticPots Array of static potentials.
		 * @param dynamicPots Array of dynamic potentials.
		 * @param waveFuncDependentPots Array of wavefunction-dependent potentials.
		 * @note The arrays are not copied, so they must remain valid for the lifetime of this object.
		 */
		CompositePotential(int nPts, int numSPots, int numDPots, int numWPots, Potential ** staticPots, Potential ** dynamicPots, Potential ** waveFuncDependentPots);
		
		~CompositePotential();
		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity();
	};

	/// Dynamically manages multiple potentials, combining them into a CompositePotential for evaluation.
	class PotentialManager :
		public Potential {
	private:
		bool compositeRefreshed = false;
		int nPts;
		std::vector<Potential*> staticPots, dynamicPots, waveFuncDependentPots;
		CompositePotential * pot=nullptr;
		PotentialComplexity myComplex = PotentialComplexity::STATIC;
		Potential ** spots = nullptr, ** dpots = nullptr, ** wpots = nullptr;
	public:
		/**
		 * Constructor.
		 * @param nPts Number of points.
		 */
		PotentialManager(int nPts);

		~PotentialManager(){if(pot) delete pot; if(spots) delete[] spots; if(dpots) delete[] dpots; if(wpots) delete[] wpots;};

		/**
		 * Add a potential to the manager.
		 * @param pot The potential to add.
		 * @note The potential is not copied, so it must remain valid for the lifetime of this object.
		 */
		void addPotential(Potential * pot);

		/// Refresh the composite potential if needed.
		void refreshCompositePotential();

		void getVBare(double t, double * targ);
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ);
		PotentialComplexity getComplexity(){return myComplex;};
	};

	/**
	 * Potential which combines a potential with a measurer, with the resulting potential being passed to the measurer instead of the total potential of the calculation.
	 * The measurer is called with the potential at each time step. The measurer should not be used elsewhere.
	 */
	class MeasuredPotential :
		public Potential {
	private:
		Potential * pot;
		Measurers::Measurer * meas;
		int numSteps;
		double maxT;
	public:
		/**
		 * Constructor.
		 * @param pot The potential to measure.
		 * @param meas The measurer to use.
		 * @param numSteps Number of steps to measure.
		 * @param maxT Maximum time to measure. This is used to find the present step index for the measurer.
		 */
		MeasuredPotential(Potential * pot, Measurers::Measurer * meas, int numSteps, double maxT) : pot(pot), meas(meas), numSteps(numSteps), maxT(maxT){};
		
		~MeasuredPotential(){};
		void getVBare(double t, double * targ){pot->getVBare(t, targ);};
		void getV(const double* rho, const std::complex<double> * psi, double t, double * targ){
			pot->getV(rho, psi, t, targ);
			meas->measure((int)(t/maxT*numSteps), psi, targ, t);
		};
		PotentialComplexity getComplexity(){return pot->getComplexity();};
	};
}