#pragma once
// Different ways of calculating weights (w*|psi(x)|^2 = rho(x) [e/m^3]) for density functional potentials
namespace WfcToRho {
	void calcEnergies(int nelec, int nPts, double dx, std::complex<double>* psi, double* totPot, KineticOperators::KineticOperator* kin, double* energies);

	// Template function for Weight (will result in error if weight is needed and this is used).
	class Weight {
	public:
		virtual void calcWeights(int nelec, double* energies, double* weights) = 0;
	};

	class FermiGasDistro :
		public Weight
	{
	private:
		double ef;
	public:
		FermiGasDistro(double ef) : ef(ef) {}
		void calcWeights(int nelec, double* energies, double* weights);
	};

	class FromDOS :
		public Weight
	{
	private:
		double ef, leff;
		boost::math::interpolators::cardinal_cubic_b_spline<double> dosISpline;
	public:
		FromDOS(double fl, double ef, double leff, const char* fil);
		void calcWeights(int nelec, double* energies, double* weights);
	};

	class Density {
	public:
		virtual void calcRho(int nPts, int nelec, double dx, double* weights, std::complex<double>* psi, double* rho) = 0;
	};

	class DirectDensity :
		public Density
	{
	private:
		int first = 1;
		double* psi2;
	public:
		void calcRho(int nPts, int nelec, double dx, double* weights, std::complex<double>* psi, double* rho);
	};

	class GaussianSmoothedDensity :
		public Density
	{
	private:
		int first = 1;
		double *psi2, *tempRho, sig;
		vtls::MaskConvolver<double>* conv;
	public:
		GaussianSmoothedDensity(double sig) : sig(sig) {}
		void calcRho(int nPts, int nelec, double dx, double* weights, std::complex<double>* psi, double* rho);
	};
}