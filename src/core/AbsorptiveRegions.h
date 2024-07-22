#pragma once
#include "CORECommonHeader.h"
// OBSOLETE
namespace AbsorptiveRegions {
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

	std::complex<double>* getSmoothedTimePhaseDecay(int len, int inner, int outer, double rate);
	double* getSmoothedSpatialDampDecay(int len, int inner, int outer, double rate);
}