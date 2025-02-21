#pragma once
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "SimulationManager.h"
#include "MathTools.h"

namespace WfcToRho {
	// Template function for Weight (will result in error if weight is needed and this is used).
	class Weight {
	public:
		virtual void calcWeights(int nElec, double* energies, double* weights) = 0;
	};

	// Fermi gas in slab system
	class BoundFermiGas :
		public Weight
	{
	private:
		double ef;
	public:
		BoundFermiGas(double ef) : ef(ef) {}
		void calcWeights(int nElec, double* energies, double* weights, NormalizationScheme norm);
	};

	// Fermi gas in semi-infinite system
	// Assumes incoming wavefunctions are normalized like 1.0*e^ikx-iwt
	class SemiInfiniteFermiGas :
		public Weight
	{
	private:
		double ef;
	public:
		SemiInfiniteFermiGas(double ef) : ef(ef) {}
		void calcWeights(int nElec, double* energies, double* weights, NormalizationScheme norm);
	};

	class FromDOS :
		public Weight
	{
	private:
		double ef, leff;
		boost::math::interpolators::cardinal_cubic_b_spline<double> dosISpline;
	public:
		FromDOS(double fl, double ef, double leff, const char* fil);
		void calcWeights(int nElec, double* energies, double* weights, NormalizationScheme norm);
	};

	class Density {
	public:
		virtual void calcRho(int nPts, int nElec, double dx, double* weights, std::complex<double>* psi, double* rho) = 0;
	};

	class DirectDensity :
		public Density
	{
	private:
		int first = 1;
		double* psi2=nullptr;
	public:
		~DirectDensity();
		void calcRho(int nPts, int nElec, double dx, double* weights, std::complex<double>* psi, double* rho);
	};

	class CylindricalDensity :
		public Density
	{
	private:
		Density* baseDens = nullptr;
		double center, radius, minX;
		int startIndex, endIndex;
		double* thinning=nullptr;
		int first = 1;
	public:
		CylindricalDensity(double center, double radius, double minX);
		~CylindricalDensity();
		void calcRho(int nPts, int nElec, double dx, double* weights, std::complex<double>* psi, double* rho);
		void doFirst(int nPts, double dx);

		void setBaseDens(Density* baseDens) { this->baseDens = baseDens; first = 1; };
	};

	class GaussianSmoothedDensity :
		public Density
	{
	private:
		int first = 1;
		double *psi2=nullptr, *tempRho=nullptr, sig;
		vtls::MaskConvolver<double>* conv = nullptr;
	public:
		GaussianSmoothedDensity(double sig) : sig(sig) {}
		~GaussianSmoothedDensity();
		void calcRho(int nPts, int nElec, double dx, double* weights, std::complex<double>* psi, double* rho);
	};
}