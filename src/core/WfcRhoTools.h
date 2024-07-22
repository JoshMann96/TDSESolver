#pragma once
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "MathTools.h"

// Different ways of calculating weights (w*|psi(x)|^2 = rho(x) [e/m^3]) for density functional potentials
namespace WfcToRho {
	void calcEnergies(int nElec, int nPts, double dx, std::complex<double>* psi, double* totPot, KineticOperators::KineticOperator* kin, double* energies);

	// Template function for Weight (will result in error if weight is needed and this is used).
	class Weight {
	public:
		virtual void calcWeights(int nElec, double* energies, double* weights) = 0;
	};

	class FermiGasDistro :
		public Weight
	{
	private:
		double ef;
	public:
		FermiGasDistro(double ef) : ef(ef) {}
		void calcWeights(int nElec, double* energies, double* weights);
	};

	class FromDOS :
		public Weight
	{
	private:
		double ef, leff;
		boost::math::interpolators::cardinal_cubic_b_spline<double> dosISpline;
	public:
		FromDOS(double fl, double ef, double leff, const char* fil);
		void calcWeights(int nElec, double* energies, double* weights);
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