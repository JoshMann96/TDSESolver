#pragma once

namespace KineticOperators {

	class KineticOperator
	{
	public:
		virtual std::complex<double>* getOperatorMatrix() = 0;
		virtual double evaluateKineticEnergy(std::complex<double>* psi) = 0;
		virtual void findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs) = 0;
	};

	class GenDisp_PSM :
		public KineticOperator
	{
	protected:
		GenDisp_PSM(int nPts, double dx, double dt) : nPts(nPts), dx(dx), dt(dt), osKineticEnergy(new std::complex<double>[nPts]) {};
	public:
		//Functions useful for updating potential immediately after kinetic phase for nonlinear systems
		//Half potential then full kinetic (returns in real space)
		void stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec);
		//Half potential
		void stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec);

		//Full OSFM step
		void stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec);

		std::complex<double>* getOperatorMatrix() {
			calcOpMat();
			return opMat;
		}

		void clearOpMat();

		void findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs);

		double evaluateKineticEnergy(std::complex<double>* psi) { return NULL; }

		void set_osKineticEnergy(std::complex<double>* kinIn) {
			vtls::copyArray(nPts, kinIn, osKineticEnergy); needMat = 1; firstStep = 1;
		}
	private:
		int firstStep = 1, firstStepMat = 1, needMat = 1;
		DFTI_DESCRIPTOR_HANDLE dftiHandle = 0, dftiHandleMat = 0;

		int nPts, nelec;
		std::complex<double> *osKineticPhase, * osPotentialPhase, *opMat;
		std::complex<double>* osKineticEnergy;
		double dx, dt;

		void calcOpMat();
		void initializeFFT();
		void initializeMatFFT();
	};

	class GenDisp_PSM_FreeElec :
		public GenDisp_PSM
	{
	public:
		GenDisp_PSM_FreeElec(int nPts, double dx, double dt, double m_eff);
	};
}