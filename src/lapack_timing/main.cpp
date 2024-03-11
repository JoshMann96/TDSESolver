#include <cstdio>
#include "CommonHeader.h"
#include "KineticOperator.h"

int main(int argc, char** arcv){
	printf("OpenBLAS lapack zhpevx_\n");
	printf("OMP Num Threads: %d\n", omp_get_max_threads());
    for(int i = 0; i < 15; i++){

		double dx = 1e-12;

        int n = std::pow(2,i);

		KineticOperators::KineticOperator * kin = new KineticOperators::GenDisp_PSM_FreeElec(n, dx, 1e-18, 1);

		double * v = new double[n];
		std::fill_n(v, n, 0.0);
		int nEigs;
		int nEigTarg = 100;

		//aim for 100 eigs
		// min energy is 0 (free electron states)
		// max energy should have 50 wavelengths within n*dx
		// \lambda = n*dx/50 -> k = 50*2*pi/(n*dx) -> maxE = 50^2 * hbar^2/2m * 4*pi^2 / (n^2 dx^2)
		double maxE = 1.0/(2.0*PhysCon::me) * std::pow(nEigTarg * PhysCon::hbar * PhysCon::pi / (n*dx),2);

		std::complex<double>* states = (std::complex<double>*)fftw_malloc(sizeof(std::complex<double>)*n*n);


		auto t1 = std::chrono::high_resolution_clock::now();

		kin->findEigenStates(v, 0, maxE, states, &nEigs);

		auto t2 = std::chrono::high_resolution_clock::now();

    	std::printf("n = %5d | nEig = %5d | t = %9.4f s\n", n, nEigs, std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0);

		std::cout << std::flush;

		fftw_free(states);
		delete[] v;
		delete kin;

    }
}