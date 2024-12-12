#include <chrono>
#include "CORECommonHeader.h"
#include "KineticOperator.h"
#include "blas.h"

void testTransparentBCs(){
    using namespace FDBCs;
    int ne = 10;
    UniqueBC* bc = new HDTransparentBC(1000, ne, 0.1*PhysCon::a0, 0.1*PhysCon::hbar/PhysCon::auE_ha);
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

void testTridiagonalAlgorithms(int rhsMethod, int lhsMethod, int nPts=40000, int nrhs=100, int nsteps=200){
	// Test timing of methods for the multiplication and inversion of tridiagonal matrices
	// rhsMethod: 0 for BLAS matrix multiplication, 1 for direct treatment, 2 for direct + OMP
	// lhsMethod: 0 for zgtsv (general tridiagonal), 1 for zgtsvx (general tridiagonal with pivoting), 2 for zptsv (positive definite tridiagonal), 3 for zptsvx (positive definite tridiagonal with pivoting)

	std::cout << "Multiplication: ";
	switch(rhsMethod){
		case 0:
			std::cout << "BLAS";
			break;
		case 1:
			std::cout << "Explicit";
			break;
		case 2:
			std::cout << "Explicit + OMP";
			break;
	}
	std::cout << std::endl << "Inversion:      ";
	switch(lhsMethod){
		case 0:
			std::cout << "zgtsv";
			break;
		case 1:
			std::cout << "zgtsvx";
			break;
		case 2:
			std::cout << "zptsv";
			break;
		case 3:
			std::cout << "zptsvx";
			break;
	}
	std::cout << std::endl;

	double dx = 0.1, dt = 0.05;
    double one = 1.0;

	// create LHS matrix, solution vector
	std::complex<double>* d = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* ud = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::complex<double>* ld = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::complex<double>* x = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nrhs);
    std::complex<double>* v = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	for(int i = 0; i < nPts * nrhs; i++)
		x[i] = std::complex<double>((double)rand() / RAND_MAX, (double)rand() / RAND_MAX);

	// LHS matrix default values
	std::complex<double>* d0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* ud0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::complex<double>* ld0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*(nPts-1));
	std::fill_n(d0, nPts, 1.0 + 0.5*PhysCon::im*dt/(dx*dx));
	std::fill_n(ud0, nPts-1, -0.25*PhysCon::im*dt/(dx*dx));
	std::fill_n(ld0, nPts-1, -0.25*PhysCon::im*dt/(dx*dx));

	// RHS matrix values
	std::complex<double>* rd = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
    std::complex<double>* rd0 = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts);
	std::complex<double>* rTemp = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>)*nPts*nrhs);
	std::complex<double> rhsDiag = 1.0 - 0.5*PhysCon::im*dt/(dx*dx), rhsSupDiag = 0.25*PhysCon::im*dt/(dx*dx);
    std::fill_n(rd0, nPts, rhsDiag);

	// test tridiagonal inversion
	std::chrono::high_resolution_clock::time_point time0, time1;
	int copyTime = 0, rhsTime = 0, invTime = 0;
	for(int i = 0; i < nsteps; i++){
        for(int k = 0; k < nPts; k++)
            v[k] = (double)rand() / RAND_MAX;

        // reset LHS matrix
        time0 = std::chrono::high_resolution_clock::now();
		vtls::copyArray(nPts, d0, d);
        vtls::addArrays(nPts, v, d);
		vtls::copyArray(nPts-1, ud0, ud);
		vtls::copyArray(nPts-1, ld0, ld);
        vtls::copyArray(nPts, rd0, rd);
        time1 = std::chrono::high_resolution_clock::now();
        copyTime += std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0).count();

		// evaluate RHS of equation
		time0 = std::chrono::high_resolution_clock::now();
		switch(rhsMethod){
			case 0: // with submethods
				vtls::copyArray(nPts*nrhs, x, rTemp);
				for(int j = 0; j < nrhs; j++){
					cblas_zscal(nPts, &rd, &x[j*nPts], 1);   // x = rd * x
                    cblas_zaxpy(nPts-1, &rhsSupDiag, &rTemp[j*nPts], 1, &x[j*nPts+1], 1); // x[:-1] = x[:-1] + rhsSupDiag[1:]  * rTemp
                    cblas_zaxpy(nPts-1, &rhsSupDiag, &rTemp[j*nPts+1], 1, &x[j*nPts], 1); // x[1:]  = x[1:]  + rhsSupDiag[:-1] * rTemp
                }
				break;
			case 1: // direct treatment
                vtls::copyArray(nPts*nrhs, x, rTemp);
				for(int j = 0; j < nrhs; j++)
                    for(int k = 1; k < nPts-1; k++)
                        x[j*nPts+k] = rd[k]*rTemp[j*nPts+k] + rhsSupDiag*rTemp[j*nPts+k-1] + rhsSupDiag*rTemp[j*nPts+k+1];
                
                for(int j = 0; j < nrhs; j++){
                    x[j*nPts] = rd[0]*rTemp[j*nPts] + rhsSupDiag*rTemp[j*nPts+1];
                    x[j*nPts+nPts-1] = rd[nPts-1]*rTemp[j*nPts+nPts-1] + rhsSupDiag*rTemp[j*nPts+nPts-2];
                }
				break;
			case 2: // direct treatment with OMP
				vtls::copyArray(nPts*nrhs, x, rTemp);
				#pragma omp parallel for collapse(2)
				for(int j = 0; j < nrhs; j++)
                    for(int k = 1; k < nPts-1; k++)
                        x[j*nPts+k] = rd[k]*rTemp[j*nPts+k] + rhsSupDiag*rTemp[j*nPts+k-1] + rhsSupDiag*rTemp[j*nPts+k+1];
                
                for(int j = 0; j < nrhs; j++){
                    x[j*nPts] = rd[0]*rTemp[j*nPts] + rhsSupDiag*rTemp[j*nPts+1];
                    x[j*nPts+nPts-1] = rd[nPts-1]*rTemp[j*nPts+nPts-1] + rhsSupDiag*rTemp[j*nPts+nPts-2];
                }
				break;
		}
		time1 = std::chrono::high_resolution_clock::now();
		rhsTime += std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count();

		// solve tridiagonal system
		time0 = std::chrono::high_resolution_clock::now();
		int info;
		switch(lhsMethod){
			case 0: // zgtsv
				LAPACK_zgtsv(&nPts, &nrhs, reinterpret_cast<dcomplex*>(ld), reinterpret_cast<dcomplex*>(d), reinterpret_cast<dcomplex*>(ud), reinterpret_cast<dcomplex*>(x), &nPts, &info);
				break;
			case 1: // zgtsvx
				throw std::runtime_error("zgtsvx not implemented");
				break;
			case 2: // zptsv
				throw std::runtime_error("zptsv not implemented");
				break;
			case 3: // zptsvx
				throw std::runtime_error("zptsvx not implemented");
				break;
		}
		time1 = std::chrono::high_resolution_clock::now();
		invTime += std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count();
	}

    // print results
    std::cout << "\tCopy time:      " << copyTime/1000 << " ms" << std::endl;
	std::cout << "\tRHS time:       " << rhsTime << " ms" << std::endl;
    std::cout << "\tInversion time: " << invTime << " ms" << std::endl;
}

int main(int argc, char** argv){
    for(int i = 0; i < 3; i++)
		for(int j = 0; j < 1; j++)
        	testTridiagonalAlgorithms(i, j);
	std::cout << "Done" << std::endl;
    return 0;
}