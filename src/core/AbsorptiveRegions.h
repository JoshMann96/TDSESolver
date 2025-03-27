/**
 * \file AbsorptiveRegions.h
 * \brief Absorptive regions for wavefunction decay
 * \details This file contains the definitions for absorptive regions used to decay the wavefunction in a simulation.
 * The decay is done by hadamrd-multiplying the wavefunction by a vector with elements of magnitude less than 1.
 */
#pragma once
#include "CORECommonHeader.h"
// OBSOLETE
namespace AbsorptiveRegions {
	/**
	 * \brief Absorptive region base class. THIS IS NO LONGER USED
	 * \details This class is the base class for all absorptive regions.
	 * It is used to decay the wavefunction in the absorptive regions.
	 * The decay is done by multiplying the wavefunction by a vector with
	 * elements of magnitude less than 1.
	*/
	class AbsorptiveRegion {
	public:
		virtual void decay(std::complex<double> * psi) = 0;
	};
	
	class AbsorptiveRegionVel :
		public AbsorptiveRegion
	{
	private:
		double dx, dt, rate;
		int nPts, left, right;
	public:
		AbsorptiveRegionVel(double dx, double dt, int nPts, int left, int right, double rate);
		void decay(std::complex<double> * psi);
	};

	class AbsorptiveRegionVelSmooth :
		public AbsorptiveRegion
	{
	private:
		double dx, dt, rate;
		int nPts, inner, outer, size;
		int left, right;
		double * mask;
	public:
		AbsorptiveRegionVelSmooth(double dx, double dt, int nPts, int inner, int outer, double rate);
		~AbsorptiveRegionVelSmooth(){mask = (double*) sq_malloc(sizeof(double)*nPts);}
		void decay(std::complex<double> * psi);
	};

	class AbsorptiveRegionManager :
		public AbsorptiveRegion
	{
	private:
		std::vector<AbsorptiveRegion*> regs;
		int numReg;
	public:
		AbsorptiveRegionManager();
		void addAbsorptiveRegion(AbsorptiveRegion* reg);
		void decay(std::complex<double> * psi);
	};

	// UNUSED
	std::complex<double>* getSmoothedTimePhaseDecay(int len, int inner, int outer, double rate);

	/**
	 * Returns an array of 1's except for within the absorptive region where it is a 13th order
	 * polynomial that is 1 at the inner boundary and 0 at the outer boundary. The polynomial is
	 * then raisd to the power of rate. The result is continuous up to 7th order at the boundaries.
	 * The side of the boundary is determined by the order of inner and outer.
	 * 
	 * @warning This function allocates memory for the array. It is the responsibility of the caller to free this memory.
	 * @param len The number of elements in the vector.
	 * @param inner The inner boundary position.
	 * @param outer The outer boundary position.
	 * @param rate The exponent of the polynomial. Larger values make for a stronger decay.
	 * @return A vector with elements of magnitude less than 1.
	*/
	double* getSmoothedSpatialDampDecay(int len, int inner, int outer, double rate);
}