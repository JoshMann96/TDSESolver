#include "CuTridiagSolver.h"

cuTridiagSolver::cuTridiagSolver(int n, int nrhs) : n(n), nrhs(nrhs), collected(true) {
    cusparseCreate(&handle);

    cudaMalloc((void**)&cDL, n * sizeof(cuDoubleComplex));
    cudaMalloc((void**)&cD, n * sizeof(cuDoubleComplex));
    cudaMalloc((void**)&cDU, n * sizeof(cuDoubleComplex));
    cudaMalloc((void**)&cB, n * nrhs * sizeof(cuDoubleComplex));

    size_t cPBufSize;
    cusparseZgtsv2_bufferSizeExt(handle, n, nrhs, cDL, cD, cDU, cB, n, &cPBufSize);
    cudaMalloc((void**)&cPBuf, cPBufSize);
}

cuTridiagSolver::~cuTridiagSolver() {
    cusparseDestroy(handle);

    cudaFree(cDL);
    cudaFree(cD);
    cudaFree(cDU);
    cudaFree(cB);
    cudaFree(cPBuf);
}

void cuTridiagSolver::solve(std::complex<double> *DL, std::complex<double> *D, std::complex<double> *DU, std::complex<double> *x) {
    // TODO: Find ways to reduce the communication overhead.
    // e.g., DL, DU never change. x could be partly reused for the next step, storing the result of the previous step.

    cudaMemcpy(&cDL[1], DL, (n-1) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(cD, D, n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(&cDU[0], DU, (n-1) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(cB, x, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    std::complex<double> zero(0.0);
    cudaMemcpy(cDL, &zero, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(&cDU[n - 1], &zero, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);

    cusparseZgtsv2(handle, n, nrhs, cDL, cD, cDU, cB, n, cPBuf);

    cudaMemcpy(x, cB, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
}

