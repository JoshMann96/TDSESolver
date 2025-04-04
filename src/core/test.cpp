#include <chrono>
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "SimulationManager.h"
#include "WfcRhoTools.h"
#include "blas.h"
#include "CuTridiagSolver.h"
#include "AbsorptiveRegions.h"

#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"

void testTransparentBCs(){
    using namespace FDBCs;
    int ne = 10;
    BoundaryCondition* bc = new UniformHDTransparentBC(1000, ne, 0.1*PhysCon::a0, 0.1*PhysCon::hbar/PhysCon::auE_ha);
    std::complex<double>* res = new std::complex<double>[ne];
    std::complex<double>* psibd = new std::complex<double>[ne];
    std::complex<double>* psiad = new std::complex<double>[ne];
    int time = 0;
    for (int i = 0; i < 100; i++){
        for (int j = 0; j < ne; j++){
            psibd[j] = std::exp(PhysCon::im*(i/100.0));
            psiad[j] = std::exp(PhysCon::im*(i/100.0+0.01));
        }
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        bc->prepareStep(psibd, psiad, i/100.0);
        bc->getRHS(psibd, psiad, i/100.0, res, ne);
        bc->finishStep(psibd, psiad, i/100.0);
        std::chrono::steady_clock::time_point ended = std::chrono::steady_clock::now();
        vtlsPrnt::printArray(ne, res);
        time += std::chrono::duration_cast<std::chrono::microseconds>(ended - begin).count();
    }
    std::cout << "calc time = " << time << " [us]" << std::endl;
    
    delete bc;
}

void testCyclicArray(){
    using namespace FDBCs;
    CyclicArray<double>* arr1 = new CyclicArray<double>(10, 2.0);
    CyclicArray<int>* arr2 = new CyclicArray<int>(10, 1);
    std::cout << arr1->inner(arr2) << std::endl;

    delete arr1;
    delete arr2;
}

void testTridiagonalAlgorithms(int nRhs=2, bool plot=true){
	// Test timing of methods for the multiplication and inversion of tridiagonal matrices
	// rhsMethod: 0 for BLAS matrix multiplication, 1 for direct treatment, 2 for direct + OMP
	// lhsMethod: 0 for zgtsv (general tridiagonal), 1 for zgtsvx (general tridiagonal with pivoting), 2 for zptsv (positive definite tridiagonal), 3 for zptsvx (positive definite tridiagonal with pivoting)

	double tmax = 400.0;
	double xmax = 200.0;

	double dx = 0.2, dt = 0.1;
    double one = 1.0;

	int nPts = std::ceil(xmax / dx);

	int nsteps = std::ceil(tmax / dt);
	int plotSteps = nsteps / 10;

	double* kins = new double[nRhs];
	for(int i = 0; i < nRhs; i++)
		kins[i] = i+1.0;

	std::cout << "nPts = " << nPts << std::endl;
	std::cout << "nsteps = " << nsteps << std::endl;

	FDBCs::BoundaryCondition* lbc = new FDBCs::UniformHDTransparentBC(1000, nRhs, dx, dt);//new FDBCs::DirichletBC((std::complex<double>)0.0);
	//FDBCs::BoundaryCondition* rbc = new FDBCs::DirichletBC((std::complex<double>)0.0);//new FDBCs::UniformHDTransparentBC(10000, nRhs, dx, dt);

	// inhomogeneous DTBC
	std::complex<double>* psibd = new std::complex<double>[nRhs];
	double* k0 = new double[nRhs];
	for (int i = 0; i < nRhs; i++){
		psibd[i] = 1.0;
		k0[i] = std::sqrt(2.0*kins[i]);
	}
	FDBCs::BoundaryCondition* rbc = new FDBCs::UniformIDTransparentBC(1000, nRhs, dx, dt, psibd, k0, 0.0);
	delete[] psibd;
	delete[] k0;

	std::complex<double> *rbct(new std::complex<double>[nRhs]), *lbct(new std::complex<double>[nRhs]), *bct1(new std::complex<double>[nRhs]), *bct2(new std::complex<double>[nRhs]);

	double* temp = (double*)sq_malloc(sizeof(double)*nPts*(nRhs+1));

	// create LHS matrix, solution vector
	std::complex<double>* d = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* ud = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::complex<double>* ld = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::complex<double>* x = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nRhs);
    std::complex<double>* v0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* v = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);

	for(int i = 0; i < nRhs; i++)
		for(int j = 0; j < nPts; j++)
			x[i*nPts+j] = 5.0*std::exp(PhysCon::im*dx*(std::sqrt(2.0*kins[i])*j))*std::exp(-dx*dx/100.0*(double)((j-nPts/2)*(j-nPts/2)));
	double norm0 = vtls::getNorm(nPts*nRhs, x, dx);

	// fill BC history
	std::complex<double> *lvs(new std::complex<double>[nRhs]), *rvs(new std::complex<double>[nRhs]);
	cblas_zcopy(nRhs, x, nPts, lvs, 1);
	cblas_zcopy(nRhs, &x[nPts-1], nPts, rvs, 1);
	std::complex<double>* phs = new std::complex<double>[nRhs];
	for(int i = 0; i < nRhs; i++)
		phs[i] = KineticOperators::CrankNicolson::phaseAdvanceFromEnergy(kins[i], dt);
	for(int i = 0; i < nRhs; i++){
		lbc->fillHistory(lvs, phs, 0.0);
		rbc->fillHistory(rvs, phs, 0.0);
	}
	delete[] kins;
	delete[] lvs;
	delete[] rvs;
	delete[] phs;

	// LHS matrix default values
	std::complex<double>* d0  = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* ud0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::complex<double>* ld0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::fill_n(d0, nPts, 1.0 + 0.5*PhysCon::im*dt/(dx*dx));
	std::fill_n(ud0, nPts-1, -0.25*PhysCon::im*dt/(dx*dx));
	std::fill_n(ld0, nPts-1, -0.25*PhysCon::im*dt/(dx*dx));
	d0[0] = lbc->getLHSEle();
	ud0[0] = lbc->getLHSAdjEle();
	d0[nPts-1] = rbc->getLHSEle();
	ld0[nPts-2] = rbc->getLHSAdjEle();

	// RHS matrix values
	std::complex<double>* rd  = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
    std::complex<double>* rd0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* x0  = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nRhs);
	std::complex<double> rhsDiag = 1.0 - 0.5*PhysCon::im*dt/(dx*dx), rhsSupDiag = 0.25*PhysCon::im*dt/(dx*dx);
    std::fill_n(rd0, nPts, rhsDiag);

	plotting::GNUPlotter* plotter;
	if(plot)
		plotter = new plotting::GNUPlotter();
				
	std::future<int> futPlot;

	// test tridiagonal inversion
	std::chrono::high_resolution_clock::time_point time0, time1;
	int copyBCTime = 0, rhsTime = 0, invTime = 0;
	for(int i = 0; i < nsteps; i++){
        for(int k = 0; k < nPts; k++)
            v0[k] = 0.0;//(nPts-1-k)*10.0/nPts;//-std::exp(-dx*dx/100.0*(double)((k-nPts/2)*(k-nPts/2)));//10.0*std::sin((30.0*i)/nsteps);//

		// status update
		if(i % plotSteps == 0){
			std::cout << "Step " << i << ": rel norm = " << vtls::getNorm(nPts*nRhs, x, dx) / norm0 <<  std::endl;
			if(plot){
				if (futPlot.valid())
					futPlot.get();
				vtls::scaMulArrayRe(nPts, 1.0, v0, temp);
				//vtls::scaMulArrayRe(nPts*nRhs, 1.0, x, &temp[nPts]);
				vtls::abs(nPts*nRhs, x, &temp[nPts]);
				futPlot = std::async(std::launch::async, [&](){
					plotter->update(nPts, nRhs+1, temp, -1.0, 10.0);
					return 0;
				});
			}
		}

		// get RHS info, update BCs
        time0 = std::chrono::high_resolution_clock::now();
		cblas_zcopy(nRhs, x, nPts, bct1, 1); // map first element of all wavefunctions to bct1
		cblas_zcopy(nRhs, &x[1], nPts, bct2, 1); // map second element of all wavefunctions to bct2

		lbc->prepareStep(bct1, bct2, std::real(v0[0]));
		lbc->getRHS(bct1, bct2, std::real(v0[0]), lbct, nRhs);
		lbc->finishStep(bct1, bct2, std::real(v0[0]));

		cblas_zcopy(nRhs, &x[nPts-1], nPts, bct1, 1); // map last element of all wavefunctions to bct1
		cblas_zcopy(nRhs, &x[nPts-2], nPts, bct2, 1); // map second to last element of all wavefunctions to bct2
		rbc->prepareStep(bct1, bct2, std::real(v0[nPts-1]));
		rbc->getRHS(bct1, bct2, std::real(v0[nPts-1]), rbct, nRhs);
		rbc->finishStep(bct1, bct2, std::real(v0[nPts-1]));

        // reset LHS matrix
		vtls::copyArray(nPts, d0, d);
		vtls::scaMulArray(nPts, 0.5*PhysCon::im*dt, v0, v);
        vtls::addArrays(nPts-2, &v[1], &d[1]); // leave out ends to leave BCs alone
		vtls::copyArray(nPts-1, ud0, ud);
		vtls::copyArray(nPts-1, ld0, ld);
        vtls::copyArray(nPts, rd0, rd);
		vtls::scaMulArray(nPts, -0.5*PhysCon::im*dt, v0, v);
		vtls::addArrays(nPts, v, rd);
		// LHS BCs
		d[0] = lbc->getLHSEle();
		ud[0] = lbc->getLHSAdjEle();
		d[nPts-1] = rbc->getLHSEle();
		ld[nPts-2] = rbc->getLHSAdjEle();

        time1 = std::chrono::high_resolution_clock::now();
        copyBCTime += std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0).count();

		// evaluate RHS of equation
		time0 = std::chrono::high_resolution_clock::now();
		vtls::copyArray(nPts*nRhs, x, x0);

		#pragma omp parallel for collapse(2)
		for(int j = 0; j < nRhs; j++)
			for(int k = 1; k < nPts-1; k++)
				x[j*nPts+k] = rd[k]*x0[j*nPts+k] + rhsSupDiag*x0[j*nPts+k-1] + rhsSupDiag*x0[j*nPts+k+1];

		/*switch(rhsMethod){
			case 0: // with submethods
				for(int j = 0; j < nRhs; j++){
					cblas_zscal(nPts, &rd, &x[j*nPts], 1);
                    cblas_zaxpy(nPts-1, &rhsSupDiag, &x0[j*nPts], 1, &x[j*nPts+1], 1); // x[:-1] += rhsSupDiag[1:]  * x0
                    cblas_zaxpy(nPts-1, &rhsSupDiag, &x0[j*nPts+1], 1, &x[j*nPts], 1); // x[1:]  += rhsSupDiag[:-1] * x0
                }
				break;
			case 1: // direct treatment
				for(int j = 0; j < nRhs; j++)
                    for(int k = 1; k < nPts-1; k++)
                        x[j*nPts+k] = rd[k]*x0[j*nPts+k] + rhsSupDiag*x0[j*nPts+k-1] + rhsSupDiag*x0[j*nPts+k+1];
                
                for(int j = 0; j < nRhs; j++){
                    x[j*nPts] = rd[0]*x0[j*nPts] + rhsSupDiag*x0[j*nPts+1];
                    x[j*nPts+nPts-1] = rd[nPts-1]*x0[j*nPts+nPts-1] + rhsSupDiag*x0[j*nPts+nPts-2];
                }
				break;
			case 2: // direct treatment with OMP
				#pragma omp parallel for collapse(2)
				for(int j = 0; j < nRhs; j++)
                    for(int k = 1; k < nPts-1; k++)
                        x[j*nPts+k] = rd[k]*x0[j*nPts+k] + rhsSupDiag*x0[j*nPts+k-1] + rhsSupDiag*x0[j*nPts+k+1];
                
                for(int j = 0; j < nRhs; j++){
                    x[j*nPts] = rd[0]*x0[j*nPts] + rhsSupDiag*x0[j*nPts+1];
                    x[j*nPts+nPts-1] = rd[nPts-1]*x0[j*nPts+nPts-1] + rhsSupDiag*x0[j*nPts+nPts-2];
                }
				break;
		}*/
		// apply RHS BCs
		cblas_zcopy(nRhs, lbct, 1, x, nPts);
		cblas_zcopy(nRhs, rbct, 1, &x[nPts-1], nPts);
		time1 = std::chrono::high_resolution_clock::now();
		rhsTime += std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0).count();

		// solve tridiagonal system: (ld, d, ud) x = rhs   ( x initially contains rhs )
		time0 = std::chrono::high_resolution_clock::now();

		int info;
		LAPACK_zgtsv(&nPts, &nRhs, reinterpret_cast<dcomplex*>(ld), reinterpret_cast<dcomplex*>(d), reinterpret_cast<dcomplex*>(ud), reinterpret_cast<dcomplex*>(x), &nPts, &info);
		/*switch(lhsMethod){
			case 0: // zgtsv
				LAPACK_zgtsv(&nPts, &nRhs, reinterpret_cast<dcomplex*>(ld), reinterpret_cast<dcomplex*>(d), reinterpret_cast<dcomplex*>(ud), reinterpret_cast<dcomplex*>(x), &nPts, &info);
				break;
			case 1: // zgtsvx
				throw std::runtime_error("zgtsvx not implemented");
				break;
			case 2: // zptsv
				throw std::runtime_error("zptsv not implemented"); // assumes positive definite, cannot be the case with TBCs
				break;
			case 3: // zptsvx
				throw std::runtime_error("zptsvx not implemented");
				break;
		}*/
		time1 = std::chrono::high_resolution_clock::now();
		invTime += std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0).count();
	}

    // print results
    std::cout << "\tCopy + BC time: " << copyBCTime/1000 << " ms" << std::endl;
	std::cout << "\tRHS time:       " << rhsTime/1000 << " ms" << std::endl;
    std::cout << "\tInversion time: " << invTime/1000 << " ms" << std::endl;
	std::cout << "\tTotal time:     " << (copyBCTime + rhsTime + invTime)/1000 << " ms" << std::endl;

	if(futPlot.valid())
		futPlot.get();

	// free memory
	delete[] rbct;
	delete[] lbct;
	delete[] bct1;
	delete[] bct2;
	sq_free(temp);
	sq_free(d);
	sq_free(ud);
	sq_free(ld);
	sq_free(x);
	sq_free(v0);
	sq_free(v);
	sq_free(d0);
	sq_free(ud0);
	sq_free(ld0);
	sq_free(rd);
	sq_free(rd0);
	sq_free(x0);
	delete lbc;
	delete rbc;

	if(plot)
		delete plotter;
}

void testCrankNicolson(){
	int nPts = 1000;
	double dx = 1e-11;
	double dt = 1e-18;
	int numSteps = 10000;
	int plotSteps = 1000;

	double* v0 = new double[nPts];
	double* xs = new double[nPts];
	double* damp = new double[nPts];

	for(int i = 0; i < nPts; i++){
		xs[i] = dx*(i-nPts/2);
		v0[i] = -PhysCon::eV*10*std::exp(-xs[i]*xs[i]/(2.0*1e-18));
		damp[i] = 1.0;
	}

	KineticOperators::CrankNicolson* cn = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, new FDBCs::DirichletBC((std::complex<double>)0.0), new FDBCs::DirichletBC((std::complex<double>)0.0));

	std::complex<double>* psi0, *psi;
	int nElec;

	cn->findEigenStates(v0, vtls::min(nPts, v0), 0.5*vtls::min(nPts, v0), &psi0, &nElec);
	psi = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nElec);

	std::cout << "nElec: " << nElec << std::endl;

	cn->setBC(new FDBCs::UniformHDTransparentBC(1000, nElec, dx, dt), FDBCs::BCSide::LEFT);
	cn->setBC(new FDBCs::UniformHDTransparentBC(1000, nElec, dx, dt), FDBCs::BCSide::RIGHT);

	plotting::GNUPlotter* plotter = new plotting::GNUPlotter();
	
	double* temp = (double*)sq_malloc(sizeof(double)*nPts*(nElec+1));
	vtls::scaMulArrayRe(nPts, 1.0, v0, temp);
	//vtls::scaMulArrayRe(nPts*nRhs, 1.0, x, &temp[nPts]);
	vtls::abs(nPts*nElec, psi0, &temp[nPts]);
	plotter->update(nPts, nElec+1, temp);

	// pause
	std::cout << "Press enter to continue..." << std::endl;
	std::cin.get();	
	
	for(int i = 0; i < nPts; i++)
		v0[i] = 1e9*PhysCon::eV*(i*dx);

	double norm0 = vtls::getNorm(nPts*nElec, psi0, dx);
	for(int i = 0; i < numSteps; i++){
		cn->step(psi0, v0, damp, psi, nElec);

		if(i % plotSteps == 0){
			std::cout << "Step " << i << ": rel norm = " << vtls::getNorm(nPts*nElec, psi, dx) / norm0 <<  std::endl;
			vtls::scaMulArrayRe(nPts, 1.0, v0, temp);
			//vtls::scaMulArrayRe(nPts*nRhs, 1.0, x, &temp[nPts]);
			vtls::abs(nPts*nElec, psi, &temp[nPts]);
			plotter->update(nPts, nElec+1, temp);
			std::cout << "Press enter to continue..." << std::endl;
			std::cin.get();
		}

		vtls::copyArray(nPts*nElec, psi, psi0);
	}
}

void testInhomogeneousEigenState(){
	int nPts = 5000;
	double dx = 0.02*PhysCon::a0;
	double dt = 0.1*PhysCon::hbar/PhysCon::auE_ha;

	double* xs = new double[nPts];
	for(int i = 0; i < nPts; i++)
		xs[i] = dx*(i-nPts/2);

	plotting::GNUPlotter* plotter = new plotting::GNUPlotter();

	SimulationManager* sm = new SimulationManager(nPts, xs[0], dx, dt);
	sm->addPotential(new Potentials::JelliumPotential(nPts, xs, 0.0, 5*PhysCon::eV, 5*PhysCon::eV, 0));
	//sm->addPotential(new Potentials::ShieldedAtomicPotential(nPts, xs, -2e-10, 4e-10, 1.5, 1e-10) );
	//sm->addPotential(new Potentials::FiniteBox(nPts, xs, -2e-9, -1e-9, -5.0*PhysCon::eV, 0));

	Measurers::Measurer* m = new Measurers::WfcRhoWeights(sm->getNElecPtr(), sm->getWeightsPtr(), "data/test");
	sm->addMeasurer(m);

	// define incoming wavefunctions
	int nElec = 50;
	std::complex<double>* psibd = new std::complex<double>[nElec];
	double* energy = new double[nElec];
	double* ks = new double[nElec];
	for (int i = 0; i < nElec; i++){
		psibd[i] = 1.0;
		energy[i] = 5.0*PhysCon::eV*std::pow((i+1.0)/nElec, 2.0);//(i+1.0)/nElec; // 0-5 eV
		//ks[i] = std::sqrt(2.0*PhysCon::me*energy[i]/PhysCon::hbar/PhysCon::hbar);
		ks[i] = KineticOperators::CrankNicolson::wavenumberFromEnergy(energy[i], 0.0, dx, dt, 1.0);
	}

	FDBCs::BoundaryCondition* rbc = new FDBCs::UniformHDTransparentBC(1000, nElec, dx, dt);
	FDBCs::BoundaryCondition* lbc = new FDBCs::UniformIDTransparentBC(1000, nElec, dx, dt, psibd, ks, 0.0);
	KineticOperators::CrankNicolson* cn = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, lbc, rbc);
	sm->setKineticOperator(cn);
	sm->setWeight(new WfcToRho::SemiInfiniteFermiGas(5.0*PhysCon::eV));
	sm->setDensity(new WfcToRho::DirectDensity());

	// plot potential
	/*double* temp = new double[nPts];
	sm->getPotPointer()->getVBare(0.0, temp);
	plotter->update(nPts, 1, temp);
	delete[] temp;*/

	sm->findInhomogeneousEigenStates(nElec, energy);

	// plot wavefunctions
	std::complex<double>* psi = sm->getPsi();
	double* temp = (double*)sq_malloc(sizeof(double)*nPts*nElec);
	vtls::normSqr(nPts*nElec, psi, temp);
	//plotter->update(nPts, nElec, temp);
	sq_free(temp);

	double* state_energies = new double[nElec];
	sm->calcEnergies(0, state_energies);
	for(int i = 0; i < nElec; i++)
		std::cout << "Expected " << i << ": " << energy[i]/PhysCon::eV << ", Got : " << state_energies[i]/PhysCon::eV << " eV" << std::endl;
	delete[] state_energies;

	// plot density and potential
	temp = (double*)sq_malloc(sizeof(double)*nPts*2);
	vtls::copyArray(nPts, sm->getRho(), temp);
	sm->getPotPointer()->getVBare(0.0, &temp[nPts]);
	vtls::scaMulArray(nPts, 1e-3/(PhysCon::eV*std::pow(PhysCon::a0,3)), &temp[nPts]);
	plotter->update(nPts, 2, temp);
	sq_free(temp);

	delete[] psibd;
	delete[] ks;
	delete[] xs;
	delete[] energy;

	// time evolution
	/*
	std::complex<double>* psi0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nElec);
	std::complex<double>* psi1 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nElec);
	std::complex<double>* initial_psi = sm->getPsi();
	vtls::copyArray(nPts*nElec, initial_psi, psi0);
	double* v = (double*)sq_malloc(sizeof(double)*nPts);
	sm->getPotPointer()->getVBare(0.0, v);
	double* damp = (double*)sq_malloc(sizeof(double)*nPts);
	std::fill_n(damp, nPts, 1.0);
	temp = (double*)sq_malloc(sizeof(double)*nPts*nElec);
	for(int j = 0; j < 100; j++){
		for(int i = 0; i < 100; i++){
			cn->step(psi0, v, damp, psi1, nElec);
			vtls::copyArray(nPts*nElec, psi1, psi0);
		}
		vtls::normSqr(nPts*nElec, psi1, temp);
		//plotter->update(nPts, nElec, temp);
		std::cout << vtls::getNorm(nPts*nElec, psi1, dx) << std::endl;
	}
	sq_free(temp);


	sq_free(psi0);
	sq_free(psi1);
	sq_free(v);
	sq_free(damp);
*/
	//delete plotter;
	delete sm;
}

void testIterationMethods(int stepType=-1, int nPts=8192){
	/* FOR RUNNING WITH WISDOM DO THIS IN MAIN
	char* wisdomFile = new char[64];
	std::snprintf(wisdomFile, 64, "fftw_nt_%04d.wisdom", omp_get_max_threads());
	fftw_init_threads();
	fftw_import_wisdom_from_filename(wisdomFile);

	//for(int i = 0; i < 3; i++)
	//	testIterationMethods(i);
	testIterationMethods();
	std::cout << "Done" << std::endl;

	fftw_export_wisdom_to_filename(wisdomFile);
	delete[] wisdomFile;
	*/

	int nSteps = 1000;
	double dx = 0.16*PhysCon::a0;
	double dt = 0.1*PhysCon::hbar/PhysCon::auE_ha;

	double* xs = new double[nPts];
	for(int i = 0; i < nPts; i++)
		xs[i] = dx*(i-nPts/2);

	SimulationManager* sm = new SimulationManager(nPts, xs[0], dx, dt);
	KineticOperators::KineticOperator* cnKin = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, new FDBCs::DirichletBC(0.0), new FDBCs::DirichletBC(0.0), true);
	KineticOperators::KineticOperator* cnKin_cpu = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, new FDBCs::DirichletBC(0.0), new FDBCs::DirichletBC(0.0), false);
	KineticOperators::KineticOperator* osKin = new KineticOperators::GenDisp_PSM_FreeElec(nPts, dx, dt, 1.0, FFTW_PATIENT);
	sm->setKineticOperator(cnKin);
	sm->setWeight(new WfcToRho::BoundFermiGas(5.0*PhysCon::eV));
	sm->setDensity(new WfcToRho::DirectDensity());
	sm->addPotential(new Potentials::FiniteBox(nPts, xs, xs[nPts/4], xs[nPts/4*3], -10.0*PhysCon::eV, 0));
	sm->addMeasurer(new Measurers::BasicMeasurers(nPts, dx, dt, "data/test/"));
	sm->addMeasurer(new Measurers::TotProb(nPts, dx, sm->getNElecPtr(), "data/test/"));
	sm->addMeasurer(new Measurers::VDProbCurrent(nPts, dx, sm->getNElecPtr(), 0, 0, "surf", "data/test/"));

	Potentials::LDAFunctional* lda_c = new Potentials::LDAFunctional(
		Potentials::LDAFunctionalType::C_PW,
		nPts, dx, 0);
	sm->addPotential(lda_c);
	Potentials::LDAFunctional* lda_x = new Potentials::LDAFunctional(
		Potentials::LDAFunctionalType::X_SLATER,
		nPts, dx, 0);
	sm->addPotential(lda_x);

	//sm->addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(nPts, nPts/10, 0, 1).get());
	//sm->addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(nPts, nPts*9/10, nPts-1, 1).get());

	if(false){
		plotting::GNUPlotter* plotter = new plotting::GNUPlotter();
		double* temp = new double[nPts];
		sm->getPotPointer()->getVBare(0.0, temp);
		plotter->update(nPts, 1, xs, temp);
		std::cout << "Press enter to continue..." << std::endl;
		std::cin.get();
		delete[] temp;
		delete plotter;
	}

	sm->addMeasurer(new Measurers::DensityPlotter(nPts, sm->getNElecPtr(), dx, xs, sm->getDensity(), sm->getWeightsPtr(), 100, false));
	sm->addMeasurer(new Measurers::PotentialPlotter(nPts, xs, 100, false));

	std::cout << "Testing the implemented iteration methods..." << std::endl;
	std::cout << "\tFinding Eigenstates..." << std::endl;
	sm->findEigenStates(-10.0*PhysCon::eV, -5.0*PhysCon::eV);

	lda_c->assemble(sm->getRho(), sm->getPsi());
	lda_x->assemble(sm->getRho(), sm->getPsi());

	// remove finite well
	sm->addPotential(new Potentials::FiniteBox(nPts, xs, xs[nPts/4], xs[nPts/4*3], 10.0*PhysCon::eV, 0));
	/*sm->addPotential(
		new Potentials::ElectricFieldProfileToPotential(
			nPts, new Potentials::ElectricFieldProfiles::ConstantFieldProfile(nPts, xs, 1e9, xs[nPts/4], xs[nPts/4*3]),
			dx, 0.0, 1.0, 800e-9, new Potentials::Envelopes::GaussianEnvelope(nSteps*dt, nSteps*dt/2.0), 0.0
		)
	);*/

	std::cout << "\tTime iterating with "  << sm->getNElec() << " wavefunctions, " << sm->getNumPoints() << " gridpoints, " << nSteps << " steps..." << std::endl;

	//std::cout << "\n\tInitializing FFTW..." << std::endl;
	//sm->runEPS_U2TU(1);

	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();


	if (stepType == -1 || stepType == 2){
		omp_set_num_threads(omp_get_max_threads());

		delete cnKin;
		cnKin = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, 
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			true);

		sm->setKineticOperator(cnKin);
		std::cout << "\tRunning FD_L GPU..." << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		sm->runCN_L(nSteps);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "\t\tTook " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
	}

	if (stepType == -1 || stepType == 3){
		omp_set_num_threads(omp_get_max_threads());

		delete cnKin_cpu;
		cnKin_cpu = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, 
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			false);

		sm->setKineticOperator(cnKin_cpu);
		std::cout << "\tRunning FD_L CPU..." << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		sm->runCN_L(nSteps);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "\t\tTook " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
	}

	if (stepType == -1 || stepType == 4){
		omp_set_num_threads(omp_get_max_threads());

		delete cnKin;
		cnKin = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, 
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			true);

		sm->setKineticOperator(cnKin);
		std::cout << "\tRunning FD_NL GPU..." << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		sm->runCN_NL(nSteps);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "\t\tTook " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
	}

	if (stepType == -1 || stepType == 5){
		omp_set_num_threads(omp_get_max_threads());

		delete cnKin_cpu;
		cnKin_cpu = new KineticOperators::CrankNicolson(nPts, dx, dt, 1.0, 
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			new FDBCs::UniformHDTransparentBC(10000, sm->getNElec(), dx, dt),
			false);

		sm->setKineticOperator(cnKin_cpu);
		std::cout << "\tRunning FD_NL CPU..." << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		sm->runCN_NL(nSteps);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "\t\tTook " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
	}

	if (stepType == -1 || stepType == 0){
		sm->setKineticOperator(osKin);
		std::cout << "\tRunning OS_U2TU..." << std::endl;
		std::cout << "\t\tInitializing FFTW..." << std::endl;
		sm->run(2); // initialize FFTW
		std::cout << "\t\tMain run..." << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		sm->run(nSteps);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "\t\tTook " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
	}

	if (stepType == -1 || stepType == 1){
		sm->setKineticOperator(osKin);
		std::cout << "\tRunning OS_UW2TUW..." << std::endl;
		std::cout << "\t\tInitializing FFTW..." << std::endl;
		sm->run(2); // initialize FFTW
		std::cout << "\t\tMain run..." << std::endl;
		t1 = std::chrono::high_resolution_clock::now();
		sm->run(nSteps);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "\t\tTook " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
	}

	std::cout << "Done!" << std::endl;

	delete sm;
	delete cnKin;
	delete cnKin_cpu;
	delete osKin;
	delete[] xs;
}

std::complex<double> randComplex(){
	return std::complex<double>(rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5);
}

void testCuTridiagSolver(){
	int n=32768, nrhs=128;
	std::complex<double> *d = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*n);
	std::complex<double> *ud = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(n-1));
	std::complex<double> *ld = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(n-1));
	std::complex<double> *x = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*n*nrhs);
	std::complex<double> *b_c = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*n*nrhs);
	std::complex<double> *b_m = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*n*nrhs);

	// fill matrix and x with random values
	for(int i = 0; i < n; i++){
		d[i] = randComplex();
		if(i < n-1){
			ud[i] = randComplex();
			ld[i] = randComplex();
		}
		for(int j = 0; j < nrhs; j++){
			x[j*n+i] = randComplex();
			//x[j*n+i] = (i == 8190) ? 1.0 : 0.0; // onehot to probe structure
		}
	}
	// fill matrix with constants, vector with onehot to probe structure
	/*for(int i = 0; i < n; i++){
		d[i] = i;
		if(i < n-1){
			ud[i] = n+i;
			ld[i] = 2*n+i;
		}
		for(int j = 0; j < nrhs; j++){
			x[j*n+i] = (i == 8191) ? 1.0 : 0.0;
		}
	}*/

	std::cout << "Testing tridiagonal matrix multiplication..." << std::endl;

	// calculate tridiagonal matrix product manually
	auto t1 = std::chrono::high_resolution_clock::now();
	for(int j = 0; j < nrhs; j++){
		b_m[j*n] = d[0]*x[j*n] + ud[0]*x[j*n+1];
		for(int i = 1; i < n-1; i++)
			b_m[j*n+i] = d[i]*x[j*n+i] + ud[i]*x[j*n+i+1] + ld[i-1]*x[j*n+i-1];
		b_m[j*n+n-1] = d[n-1]*x[j*n+n-1] + ld[n-2]*x[j*n+n-2];
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "\tManual product " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;

	// calculate using cblas
	std::complex<double>* amat = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*3*n);
	std::fill_n(amat, 3*n, 0.0);
	cblas_zcopy(n-1, ld, 1, amat+2, 3);
	cblas_zcopy(n-1, ud, 1, amat+3, 3);
	t1 = std::chrono::high_resolution_clock::now();
	cblas_zcopy(n, d, 1, amat+1, 3);
	std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);
	cblas_zgbmv(CblasColMajor, CblasNoTrans, n, n, 1, 1, &alpha, amat, 3, x, 1, &beta, b_c, 1);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "\tBLAS product took " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;

	//vtlsPrnt::printArray(n*nrhs, b_m);

	// calculate tridiagonal matrix product with cudaTridiagonalSolverSystem
	cudaTridiagonalSolverSystem* solver = new cudaTridiagonalSolverSystem(n, nrhs);
	solver->setX(x);
	solver->setOffDiag(ld, ud, cudaTridiagonalSolverSystem::RHS);
	t1 = std::chrono::high_resolution_clock::now(); // the present state will be stored on the device, setX and setOffDiag is not called repeatedly
	solver->rhsProduct(d, true);
	t2 = std::chrono::high_resolution_clock::now();
	solver->gatherRHS(b_c, true);
	std::cout << "\tCUDA product took " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;
	
	// compare results
	std::cout << "\tChecking for errors..." << std::endl;
	for(int j = 0; j < nrhs; j++)
		for(int i = 0; i < n; i++)
			if(std::abs(b_c[j*n+i] - b_m[j*n+i]) > 1e-10)
				std::cout << "\t\tMismatch at " << i << ", " << j << " : CUDA != CPU : " << b_c[j*n+i] << " != " << b_m[j*n+i] << std::endl;

	std::cout << "Testing tridiagonal matrix inversion (undoing product)..." << std::endl;
	std::complex<double>* ldt = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(n-1));
	std::complex<double>* udt = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(n-1));
	std::complex<double>* dt = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*n);
	// cpu
	t1 = std::chrono::high_resolution_clock::now();
	cblas_zcopy(n-1, ld, 1, ldt, 1);
	cblas_zcopy(n-1, ud, 1, udt, 1);
	cblas_zcopy(n, d, 1, dt, 1);
	int info;
	LAPACK_zgtsv(&n, &nrhs, reinterpret_cast<dcomplex*>(ldt), reinterpret_cast<dcomplex*>(dt), reinterpret_cast<dcomplex*>(udt), reinterpret_cast<dcomplex*>(b_m), &n, &info);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "\tLAPACK inversion took " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;

	// gpu
	solver->setOffDiag(ld, ud, cudaTridiagonalSolverSystem::LHS);
	t1 = std::chrono::high_resolution_clock::now();
	solver->solve(d, false, true);
	solver->gatherX(b_c, false);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "\tCUDA inversion took " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;
	
	// compare results, should be same as original vector
	std::cout << "\tChecking for errors..." << std::endl;
	// cpu
	for(int j = 0; j < nrhs; j++)
		for(int i = 0; i < n; i++)
			if(std::abs(x[j*n+i] - b_m[j*n+i]) > 1e-10)
				std::cout << "\t\tMismatch at " << i << ", " << j << " : CPU != EXPCTD : " << b_m[j*n+i] << " != " << x[j*n+i] << std::endl;
	// gpu
	for(int j = 0; j < nrhs; j++)
		for(int i = 0; i < n; i++)
			if(std::abs(x[j*n+i] - b_c[j*n+i]) > 1e-10)
				std::cout << "\t\tMismatch at " << i << ", " << j << " : CUDA != EXPCTD : " << b_c[j*n+i] << " != " << x[j*n+i] << std::endl;

	// test finding the density
	std::cout << "Testing calculating the density..." << std::endl;
	double* rho = (double*)sq_malloc(sizeof(double)*n);
	double* rho_c = (double*)sq_malloc(sizeof(double)*n);
	double* weights = (double*)sq_malloc(sizeof(double)*nrhs);
	std::fill_n(weights, nrhs, 1.0);

	// cpu
	t1 = std::chrono::high_resolution_clock::now();
	#pragma omp parallel for
	for(int i = 0; i < n; i++){
		rho[i] = 0;
		for(int j = 0; j < nrhs; j++)
			rho[i] += weights[j]*std::abs(x[j*n+i]*x[j*n+i]);
	}
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "\tCPU density took " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;

	// gpu
	t1 = std::chrono::high_resolution_clock::now();
	solver->calcRawRho(weights, rho_c, false);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "\tCUDA density took " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;

	// compare results
	std::cout << "\tChecking for errors..." << std::endl;
	for(int i = 0; i < n; i++)
		if(std::abs(rho[i] - rho_c[i]) > 1e-10)
			std::cout << "\t\tMismatch at " << i << " : CUDA != CPU : " << rho_c[i] << " != " << rho[i] << std::endl;

	std::cout << "Done!" << std::endl;

	delete solver;
	sq_free(d);
	sq_free(ud);
	sq_free(ld);
	sq_free(x);
	sq_free(b_c);
	sq_free(b_m);
	sq_free(amat);
	sq_free(ldt);
	sq_free(udt);
	sq_free(dt);
	sq_free(rho);
	sq_free(rho_c);
}


int main(int argc, char** argv){
	// char* wisdomFile = new char[64];
	// std::snprintf(wisdomFile, 64, "fftw_nt_%04d.wisdom", omp_get_max_threads());
	// fftw_init_threads();
	// fftw_import_wisdom_from_filename(wisdomFile);

	// testIterationMethods(-1, 2048);
	// testIterationMethods(-1, 4096);
	// testIterationMethods(-1, 8192);
	// testIterationMethods(-1, 16384);
	// testIterationMethods(-1, 32768);

	// fftw_export_wisdom_to_filename(wisdomFile);
	// delete[] wisdomFile;

	// std::cout << "Done" << std::endl;

	testIterationMethods(2, 2048);
	testIterationMethods(3, 2048);
	testIterationMethods(4, 2048);
	testIterationMethods(5, 2048);

    return 0;
}