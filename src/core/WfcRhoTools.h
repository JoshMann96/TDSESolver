#pragma once
#include "CORECommonHeader.h"
#include "MathTools.h"
#include "KineticOperator.h"

namespace WfcToRho {

	enum NormalizationScheme{
		UNNORMALIZED,
		NORMALIZED
	};

	// Template function for Weight (will result in error if weight is needed and this is used).
	class Weight {
	public:
		virtual void calcWeights(int nElec, double* energies, double* weights, NormalizationScheme norm) = 0;
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
	private:
		double* psi2_work = nullptr;
	public:
		~Density() { if (psi2_work) sq_free(psi2_work); }
		static void calcRawRho(int nPts, int nElec, const double* weights, const std::complex<double>* psi, double* psi2_work, double* rho);
		void calcRho(int nPts, int nElec, double dx, const double* weights, const std::complex<double>* psi, double* rho){
			if (!psi2_work) psi2_work = (double*)sq_malloc(sizeof(double) * nPts * nElec);
			calcRawRho(nPts, nElec, weights, psi, psi2_work, rho);
			calcRho(nPts, nElec, dx, rho);
		};
		virtual void calcRho(int nPts, int nElec, double dx, double* rho) = 0; // rho is in/out (raw rho then processed rho)
	};

	class DirectDensity :
		public Density
	{
	private:
		int first = 1;
	public:
		void calcRho(int nPts, int nElec, double dx, double* rho);
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
		void calcRho(int nPts, int nElec, double dx, double* rho);
		void doFirst(int nPts, double dx);

		void setBaseDens(Density* baseDens) { this->baseDens = baseDens; first = 1; };
	};

	class GaussianSmoothedDensity :
		public Density
	{
	private:
		int first = 1;
		double *tempRho=nullptr, sig;
		vtls::MaskConvolver<double>* conv = nullptr;
	public:
		GaussianSmoothedDensity(double sig) : sig(sig) {}
		~GaussianSmoothedDensity();
		void calcRho(int nPts, int nElec, double dx, double* rho);
	};
}