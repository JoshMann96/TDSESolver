#include <chrono>
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "blas.h"

void testTransparentBCs(){
    using namespace FDBCs;
    int ne = 10;
    BoundaryCondition* bc = new HDTransparentBC(1000, ne, 0.1*PhysCon::a0, 0.1*PhysCon::hbar/PhysCon::auE_ha);
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

void testTridiagonalAlgorithms(int nRhs=2, bool plot=false){
	// Test timing of methods for the multiplication and inversion of tridiagonal matrices
	// rhsMethod: 0 for BLAS matrix multiplication, 1 for direct treatment, 2 for direct + OMP
	// lhsMethod: 0 for zgtsv (general tridiagonal), 1 for zgtsvx (general tridiagonal with pivoting), 2 for zptsv (positive definite tridiagonal), 3 for zptsvx (positive definite tridiagonal with pivoting)

	double tmax = 100.0;
	double xmax = 200.0;

	double dx = 0.2, dt = 0.1;
    double one = 1.0;

	int nPts = std::ceil(xmax / dx);

	int nsteps = std::ceil(tmax / dt);
	int plotSteps = nsteps / 10;

	double* kins = new double[nRhs];
	for(int i = 0; i < nRhs; i++)
		kins[i] = 4.0;

	FDBCs::BoundaryCondition* lbc = new FDBCs::HDTransparentBC(10000, nRhs, dx, dt);//new FDBCs::DirichletBC((std::complex<double>)0.0);
	FDBCs::BoundaryCondition* rbc = new FDBCs::HDTransparentBC(10000, nRhs, dx, dt);//new FDBCs::DirichletBC((std::complex<double>)0.0);
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
			x[i*nPts+j] = 5.0*std::exp(-PhysCon::im*dx*(std::sqrt(2.0*kins[i])*j))*std::exp(-dx*dx/100.0*(double)((j-nPts/2)*(j-nPts/2)));
	double norm0 = vtls::getNorm(nPts*nRhs, x, dx);

	// fill BC history
	std::complex<double> *lvs(new std::complex<double>[nRhs]), *rvs(new std::complex<double>[nRhs]);
	cblas_zcopy(nRhs, x, nPts, lvs, 1);
	cblas_zcopy(nRhs, &x[nPts-1], nPts, rvs, 1);
	for(int i = 0; i < nRhs; i++){
		lbc->fillHistory(lvs, kins);
		rbc->fillHistory(rvs, kins);
	}
	delete[] kins;
	delete[] lvs;
	delete[] rvs;

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
            v0[k] = 0.1*dt*i;//-std::exp(-dx*dx/100.0*(double)((k-nPts/2)*(k-nPts/2)));

		// status update
		if(i % plotSteps == 0){
			std::cout << "Step " << i << ": rel norm = " << vtls::getNorm(nPts*nRhs, x, dx) / norm0 <<  std::endl;
			if(plot){
				if (futPlot.valid())
					futPlot.get();
				vtls::scaMulArrayRe(nPts, 1.0, v0, temp);
				vtls::scaMulArrayRe(nPts*nRhs, 1.0, x, &temp[nPts]);
				futPlot = std::async(std::launch::async, [&](){
					plotter->update(nPts, nRhs+1, temp, -5.0, 5.0);
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

int main(int argc, char** argv){
    testTridiagonalAlgorithms();
	std::cout << "Done" << std::endl;
    return 0;
}