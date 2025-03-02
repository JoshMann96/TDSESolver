#pragma once
#include "CORECommonHeader.h"
#include "MathTools.h"
#include "FDBCs.h"

namespace KineticOperators {

	class KineticOperator
	{
	public:
		virtual double evaluateKineticEnergy(std::complex<double>* psi) = 0;
		virtual void findEigenStates(double* v, double emin, double emax, std::complex<double>** states, int* nEigs) = 0;
	};

	class KineticOperator_PSM :
		public KineticOperator
	{
	protected:
		uint fftwPlanPolicy;
		KineticOperator_PSM(uint fftwPlanPolicy) : fftwPlanPolicy(fftwPlanPolicy) {};
	public:
		//Functions useful for updating potential immediately after kinetic phase for nonlinear systems
		//Half potential then full kinetic (returns in real space)
		virtual void stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) = 0;
		//Half potential
		virtual void stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) = 0;

		//Full OSFM step
		virtual void stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) = 0;
	};


	class GenDisp_PSM :
		public KineticOperator_PSM
	{
	protected:
		GenDisp_PSM(int nPts, double dx, double dt, uint fftwPlanPolicy=FFTW_PATIENT) : nPts(nPts), dx(dx), dt(dt), osKineticEnergy((std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts)), KineticOperator_PSM(fftwPlanPolicy) {};
	public:
		~GenDisp_PSM();
		//Functions useful for updating potential immediately after kinetic phase for nonlinear systems
		//Half potential then full kinetic (returns in real space)
		void stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec);
		//Half potential
		void stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec);

		//Full OSFM step
		void stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec);

		void freeOpMat() {
			if (opMat)
				sq_free(opMat); opMat = nullptr;
			needMat = 1;
		}

		void findEigenStates(double* v, double emin, double emax, std::complex<double>** states, int* nEigs);

		double evaluateKineticEnergy(std::complex<double>* psi);

		void set_osKineticEnergy(std::complex<double>* kinIn) {
			vtls::copyArray(nPts, kinIn, osKineticEnergy); needMat = 1;
		}
	private:
		int firstStepAll = 1, firstStepOne = 1, needMat = 1;
		//DFTI_DESCRIPTOR_HANDLE dftiHandle = 0, dftiHandleMat = 0, dftiHandleKin = 0;
		fftw_plan fftwAllForward=NULL, fftwAllBackward=NULL, fftwOneForward=NULL, fftwOneBackward=NULL;

		int nPts, nElec;
		std::complex<double> *osKineticPhase = nullptr, * osPotentialPhase = nullptr, *opMat = nullptr;
		std::complex<double>* osKineticEnergy = nullptr;
		double dx, dt;

		void calcOpMat();
		void initializeAllFFT(int nElec);
		void initializeOneFFT();
		void executeAllFFTForward(std::complex<double>* targ);
		void executeAllFFTBackward(std::complex<double>* targ);
		void executeOneFFTForward(std::complex<double>* targ);
		void executeOneFFTBackward(std::complex<double>* targ);
		
	};

	class GenDisp_PSM_FreeElec :
		public GenDisp_PSM
	{
	public:
		GenDisp_PSM_FreeElec(int nPts, double dx, double dt, double m_eff, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	class GenDisp_PSM_Series :
		public GenDisp_PSM
	{
	public:
		GenDisp_PSM_Series(int nPts, double dx, double dt, int nPoly, double* polyCoeffs, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	class GenDisp_PSM_MathExpr :
		public GenDisp_PSM
	{
	public:
		GenDisp_PSM_MathExpr(int nPts, double dx, double dt, std::string expr, uint fftwPlanPolicy=FFTW_PATIENT);
	};


	class NonUnifGenDisp_PSM :
		public KineticOperator_PSM
	{
	protected:
		NonUnifGenDisp_PSM(int nPts, double dx, double dt, int nDisp, int expOrder, int forceNormalization, uint fftwPlanPolicy=FFTW_PATIENT) : 
			nPts(nPts), dx(dx), dt(dt), nDisp(nDisp), expOrder(expOrder), forceNorm(forceNormalization), 
			osKineticEnergy((std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nDisp)),
			osKineticMask((double*)sq_malloc(sizeof(double)*nPts*nDisp)),
			KineticOperator_PSM(fftwPlanPolicy) {};
	public:
		~NonUnifGenDisp_PSM();
		//Functions useful for updating potential immediately after kinetic phase for nonlinear systems
		//Half potential then full kinetic (returns in real space)
		void stepOS_UW2T(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec);
		//Half potential
		void stepOS_UW(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec);

		//Full OSFM step
		void stepOS_U2TU(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec);

		void freeOpMat() {
			if (opMat)
				sq_free(opMat); opMat = NULL;
			needMat = 1;
		}

		void findEigenStates(double* v, double emin, double emax, std::complex<double>** states, int* nEigs);

		double evaluateKineticEnergy(std::complex<double>* psi);

		void set_osKineticEnergy(std::complex<double>* kinIn, double* maskIn) {
			vtls::copyArray(nPts * nDisp, kinIn, osKineticEnergy); needMat = 1; firstStepOne = 1;
			vtls::copyArray(nPts * nDisp, maskIn, osKineticMask);
			//take square root, as is required for this method
			for (int i = 0; i < nPts * nDisp; i++)
				osKineticEnergy[i] = std::sqrt(osKineticEnergy[i]);
		}
	private:
		int firstStepAll = 1, firstStepOne = 1, needMat = 1;
		//DFTI_DESCRIPTOR_HANDLE dftiHandle = 0, dftiHandleMat = 0, dftiHandleKin = 0;
		fftw_plan fftwAllForward=NULL, fftwAllBackward=NULL, fftwOneForward=NULL, fftwOneBackward=NULL;

		int nPts, nElec, nDisp, expOrder, forceNorm;
		std::complex<double>* osPotentialPhase = nullptr, * opMat = nullptr;
		std::complex<double>* osKineticEnergy = nullptr;
		std::complex<double>* tempPsi = nullptr, *tempPsiCum = nullptr;
		double* osKineticMask = nullptr, *norms = nullptr;
		double dx, dt;

		void calcOpMat();
		void initializeAllFFT(int nElec);
		void initializeOneFFT();
		void executeAllFFTForward(std::complex<double>* targ);
		void executeAllFFTBackward(std::complex<double>* targ);
		void executeOneFFTForward(std::complex<double>* targ);
		void executeOneFFTBackward(std::complex<double>* targ);
	};

	class NonUnifGenDisp_PSM_EffMassBoundary :
		public NonUnifGenDisp_PSM
	{
	public:
		//meff_r and meff_l are relative effective masses (1 for electron rest mass)
		NonUnifGenDisp_PSM_EffMassBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, double meff_l, double meff_r, double transRate, int transPos, double edgeRate, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	class NonUnifGenDisp_PSM_MathExprBoundary :
		public NonUnifGenDisp_PSM
	{
	public:
		//different regions with different dispersion relations as text, in order from left to right
		NonUnifGenDisp_PSM_MathExprBoundary(int nPts, double dx, double dt, int expOrder, int forceNormalization, int nDisp, std::vector<std::string> exprs, double* transRates, int* transPoss, uint fftwPlanPolicy=FFTW_PATIENT);
	};

	class KineticOperator_FDM :
		public KineticOperator
	{
	protected:
		FDBCs::BoundaryCondition *lbc, *rbc;
		int nPts;
	public:
		KineticOperator_FDM(int nPts, FDBCs::BoundaryCondition* lbc, FDBCs::BoundaryCondition* rbc) : nPts(nPts), lbc(lbc), rbc(rbc) {};
		virtual void stepVirtual(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) = 0; // timestep without iterating BCs
		virtual void step(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec) = 0;
		virtual void findInhomogeneousEigenStates(double* v, double* es, std::complex<double>* states, int nElec) = 0;
		void projectHistory(std::complex<double>* psi, std::complex<double>* phsL, std::complex<double>* phsR, double* v, int nElec);

		void setBC(FDBCs::BoundaryCondition* bc, FDBCs::BCSide side) { 
			switch (side) {
				case FDBCs::BCSide::LEFT:
					lbc = bc;
					break;
				case FDBCs::BCSide::RIGHT:
					rbc = bc;
					break;
			}
		};
	};

	class CrankNicolson:
		public KineticOperator_FDM
	{
		double dx, dt, m_eff;
		std::complex<double> lhsOffDiag0, lhsDiag0, rhsDiag0, rhsOffDiag, potmul; // elements of LHS tridiagonal matrix
		std::complex<double> *d, *ud, *ld;

		std::complex<double> *rbct=nullptr, *lbct=nullptr, *bct1=nullptr, *bct2=nullptr;
	
		void _step(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec, int isVirtual);
	public:
		CrankNicolson(int nPts, double dx, double dt, double m_eff, FDBCs::BoundaryCondition* leftBC, FDBCs::BoundaryCondition* rightBC);
		~CrankNicolson(){
			sq_free(d);
			sq_free(ud);
			sq_free(ld);

			if(rbct)
				sq_free(rbct);
			if(lbct)
				sq_free(lbct);
			if(bct1)
				sq_free(bct1);
			if(bct2)
				sq_free(bct2);


		};

		void step(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec){
			_step(psi0, v, spatialDamp, targ, nElec, 0);
		};
		void stepVirtual(std::complex<double>* psi0, double* v, double* spatialDamp, std::complex<double>* targ, int nElec){
			_step(psi0, v, spatialDamp, targ, nElec, 1);
		}
		void findEigenStates(double* v, double emin, double emax, std::complex<double>** states, int* nEigs);
		// iteratively solve inhomogeneous system, k0s and v0s are wavenumbers and potential values for the initial states
		void findInhomogeneousEigenStates(double* v, double* es, std::complex<double>* states, int nElec);
		double evaluateKineticEnergy(std::complex<double>* psi);

		static double wavenumberFromPhase(std::complex<double> phase, double v, double dx, double dt, double m_eff){
			dx /= PhysCon::a0;
			dt *= PhysCon::auE_ha/PhysCon::hbar;
			v  /= PhysCon::auE_ha;
			
			double cosine = 1.0 - m_eff*dx*dx*( 2.0/dt*std::tan(std::arg(phase)/2.0) - v ) ;

			if(std::abs(cosine) > 1.0)
				throw std::runtime_error("Crank-Nicolson iteration phase is too large.");
			return 1.0/dx/PhysCon::a0 * std::acos(cosine);
		};

		static double wavenumberFromEnergy(double energy, double v, double dx, double dt, double m_eff){
			energy /= PhysCon::auE_ha;
			dx /= PhysCon::a0;
			dt *= PhysCon::auE_ha/PhysCon::hbar;
			v  /= PhysCon::auE_ha;
			
			//std::cout << "energy: " << energy << std::endl;
			//std::cout << "v: " << v << std::endl;

			double cosine = 1.0 - m_eff*dx*dx*( energy - v ) ;

			//std::cout << "cosine: " << cosine << std::endl;

			if(std::abs(cosine) > 1.0)
				throw std::runtime_error("Crank-Nicolson iteration phase is too large.");
			return 1.0/dx/PhysCon::a0 * std::acos(cosine);
		};

		static std::complex<double> phaseAdvanceFromWavenumber(double k0, double v, double dx, double dt, double m_eff){
			k0 *= PhysCon::a0;
			dx /= PhysCon::a0;
			dt *= PhysCon::auE_ha/PhysCon::hbar;
			v  /= PhysCon::auE_ha;
			std::complex<double> num = 1.0-0.5*dt*PhysCon::im*( (1.0-std::cos(k0*dx))/dx/dx/m_eff + v );
			return num/std::conj(num);
		};

		static std::complex<double> phaseAdvanceFromEnergy(double e, double dt){
			std::complex<double> num = 1.0-0.5*dt*PhysCon::im*e/PhysCon::hbar;
			return num/std::conj(num);
		};
	};
}