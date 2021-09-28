#pragma once

namespace KineticOperators {

	class KineticOperator
	{
	public:
		virtual std::complex<double>* getOperatorMatrix() = 0;
		virtual double evaluateKineticEnergy(std::complex<double>* psi) = 0;
		virtual void findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs) = 0;
	};

	class KineticOperator_PSM :
		public KineticOperator
	{
	public:
		//Functions useful for updating potential immediately after kinetic phase for nonlinear systems
		//Half potential then full kinetic (returns in real space)
		virtual void stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) = 0;
		//Half potential
		virtual void stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) = 0;

		//Full OSFM step
		virtual void stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nelec) = 0;
	};


	class GenDisp_PSM :
		public KineticOperator_PSM
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

		void clearOpMat() {
			if (opMat)
				delete[] opMat; opMat = NULL;
			needMat = 1;
		}

		void findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs);

		double evaluateKineticEnergy(std::complex<double>* psi);

		void set_osKineticEnergy(std::complex<double>* kinIn) {
			vtls::copyArray(nPts, kinIn, osKineticEnergy); needMat = 1; firstStep = 1;
		}
	private:
		int firstStep = 1, firstStepMat = 1, firstStepKin = 1, needMat = 1;
		DFTI_DESCRIPTOR_HANDLE dftiHandle = 0, dftiHandleMat = 0, dftiHandleKin = 0;

		int nPts, nelec;
		std::complex<double> *osKineticPhase, * osPotentialPhase, *opMat;
		std::complex<double>* osKineticEnergy;
		std::complex<double>* temp1, *temp2;
		double dx, dt;

		void calcOpMat();
		void initializeFFT(int nelec);
		void initializeMatFFT();
		void initializeKinFFT();
	};

	class GenDisp_PSM_FreeElec :
		public GenDisp_PSM
	{
	public:
		GenDisp_PSM_FreeElec(int nPts, double dx, double dt, double m_eff);
	};

	class GenDisp_PSM_Series :
		public GenDisp_PSM
	{
	public:
		GenDisp_PSM_Series(int nPts, double dx, double dt, int nPoly, double* polyCoeffs);
	};

	class GenDisp_PSM_MathExpr :
		public GenDisp_PSM
	{
	public:
		GenDisp_PSM_MathExpr(int nPts, double dx, double dt, std::string expr);
	};


	class NonUnifGenDisp_PSM :
		public KineticOperator_PSM
	{
	protected:
		NonUnifGenDisp_PSM(int nPts, double dx, double dt, int nDisp, int expOrder, int forceNormalization) : 
			nPts(nPts), dx(dx), dt(dt), nDisp(nDisp), expOrder(expOrder), forceNorm(forceNormalization), 
			osKineticEnergy(new std::complex<double>[nPts * nDisp]), osKineticMask(new double[nPts * nDisp]) {};
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

		void clearOpMat() {
			if (opMat)
				delete[] opMat; opMat = NULL;
			needMat = 1;
		}

		void findEigenStates(double* v, double emin, double emax, std::complex<double>** psi, int* nEigs);

		double evaluateKineticEnergy(std::complex<double>* psi);

		void set_osKineticEnergy(std::complex<double>* kinIn, double* maskIn) {
			vtls::copyArray(nPts * nDisp, kinIn, osKineticEnergy); needMat = 1; firstStep = 1;
			vtls::copyArray(nPts * nDisp, maskIn, osKineticMask);
			//take square root, as is required for this method
			for (int i = 0; i < nPts * nDisp; i++)
				osKineticEnergy[i] = std::sqrt(osKineticEnergy[i]);
		}
	private:
		int firstStep = 1, firstStepMat = 1, firstStepKin = 1, needMat = 1;
		DFTI_DESCRIPTOR_HANDLE dftiHandle = 0, dftiHandleMat = 0, dftiHandleKin = 0;

		int nPts, nelec, nDisp, expOrder, forceNorm;
		std::complex<double>* osPotentialPhase, * opMat;
		std::complex<double>* osKineticEnergy;
		std::complex<double>* tempPsi, *tempPsiCum;
		std::complex<double>* temp1, * temp2, *temp3;
		double* osKineticMask, *norms;
		double dx, dt;

		void calcOpMat();
		void initializeFFT(int nelec);
		void initializeMatFFT();
		void initializeKinFFT();
	};

	class NonUnifGenDisp_PSM_EffMassBoundary :
		public NonUnifGenDisp_PSM
	{
	public:
		//meff_r and meff_l are relative effective masses (1 for electron rest mass)
		NonUnifGenDisp_PSM_EffMassBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, double meff_l, double meff_r, double transRate, int transPos);
	};

	class NonUnifGenDisp_PSM_MathExprBoundary :
		public NonUnifGenDisp_PSM
	{
	public:
		//different regions with different dispersion relations as text, in order from left to right
		NonUnifGenDisp_PSM_MathExprBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, int nDisp, std::vector<std::string> exprs, double* transRates, int* transPoss);
	};


	class KineticOperator_FDM :
		public KineticOperator
	{

	};
}