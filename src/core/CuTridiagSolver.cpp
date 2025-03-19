#include "CuTridiagSolver.h"
#include "cudaTools.cuh"

cudaTridiagonalSolverSystem::cudaTridiagonalSolverSystem(int n, int nrhs) : n(n), nrhs(nrhs), lhsOffdiagDefined(false), rhsOffdiagDefined(false) {
    cudaStatCheck(  cusparseCreate(&csHandle));
    cudaStatCheck(  cublasCreate_v2(&cbHandle));
    
    // LHS matrix
    cudaStatCheck(  cudaMalloc((void**)&cDL, n * sizeof(cuDoubleComplex)));
    cudaStatCheck(  cudaMalloc((void**)&cD, n * sizeof(cuDoubleComplex)));
    cudaStatCheck(  cudaMalloc((void**)&cDU, n * sizeof(cuDoubleComplex)));

    // solution vector
    cudaStatCheck(  cudaMalloc((void**)&_cX, n * nrhs * sizeof(cuDoubleComplex)));
    cudaStatCheck(  cudaMalloc((void**)&_cXV, n * nrhs * sizeof(cuDoubleComplex)));
    cX.data = _cX;
    cXV.data = _cXV;

    // density
    cudaStatCheck(  cudaMalloc((void**)&cRho, n * sizeof(double)));
    cudaStatCheck(  cudaMalloc((void**)&cWeights, nrhs * sizeof(double)));

    // workspace
    size_t cPBufSize;
    cudaStatCheck(  cusparseZgtsv2_bufferSizeExt(csHandle, n, nrhs, cDL, cD, cDU, cX.data, n, &cPBufSize));
    cudaStatCheck(  cudaMalloc((void**)&cPBuf, cPBufSize));

    // RHS matrix (banded format, tridiagonal)
    cudaStatCheck(  cudaMalloc((void**)&cRHSMat, n * 3 * sizeof(cuDoubleComplex)));
    cudaStatCheck(  cudaMemset(cRHSMat, 0, n * 3 * sizeof(cuDoubleComplex)));
    cudaStatCheck(  cudaMalloc((void**)&rhsTemp, n * sizeof(cuDoubleComplex)));
}

cudaTridiagonalSolverSystem::~cudaTridiagonalSolverSystem() {
    cudaStatCheck(  cusparseDestroy(csHandle));
    cudaStatCheck(  cublasDestroy_v2(cbHandle));

    cudaStatCheck(  cudaFree(cDL));
    cudaStatCheck(  cudaFree(cD));
    cudaStatCheck(  cudaFree(cDU));

    cudaStatCheck(  cudaFree(cRHSMat));
    cudaStatCheck(  cudaFree(rhsTemp));

    cudaStatCheck(  cudaFree(_cX));
    cudaStatCheck(  cudaFree(_cXV));
    if(tempState != nullptr)
        cudaStatCheck(  cudaFree(tempState));

    cudaStatCheck(  cudaFree(cRho));
    cudaStatCheck(  cudaFree(cWeights));

    cudaStatCheck(  cudaFree(cPBuf));
}

void cudaTridiagonalSolverSystem::setOffDiag(const std::complex<double>* DL, const std::complex<double>* DU, Side side) {
    if (side == LHS){
        cudaStatCheck(  cudaMemcpy(cDL+1, DL, (n-1) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
        cudaStatCheck(  cudaMemcpy(cDU, DU, (n-1) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

        cudaStatCheck(  cudaMemset(cDL, 0, sizeof(cuDoubleComplex))); // set first element to zero
        cudaStatCheck(  cudaMemset(cDU+n-1, 0, sizeof(cuDoubleComplex))); // set last element to zero

        lhsOffdiagDefined = true;
    }
    else { // RHS
        cudaStatCheck(  cudaMemcpy(rhsTemp, DL, (n-1) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
        cudaStatCheck(  cublasZcopy_v2(cbHandle, n-1, rhsTemp, 1, cRHSMat+2, 3));

        cudaStatCheck(  cudaMemcpy(rhsTemp, DU, (n-1) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
        cudaStatCheck(  cublasZcopy_v2(cbHandle, n-1, rhsTemp, 1, cRHSMat+3, 3));

        rhsOffdiagDefined = true;
    }
}

void cudaTridiagonalSolverSystem::solve(const std::complex<double> *DL, const std::complex<double> *D, const std::complex<double> *DU, std::complex<double> *x) {
    setOffDiag(DL, DU, LHS);
    solve(D, x);

    lhsOffdiagDefined = true;
}

void cudaTridiagonalSolverSystem::solve(const std::complex<double> *D, std::complex<double> *x) {
    if (!lhsOffdiagDefined)
        throw std::runtime_error("cudaTridiagonalSolverSystem::solve : LHS off-diagonal elements are not defined.");

    cudaStatCheck(  cudaMemcpy(cD, D, n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    cudaStatCheck(  cudaMemcpy(cX.data, x, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    
    cudaStatCheck(  cusparseZgtsv2(csHandle, n, nrhs, cDL, cD, cDU, cX.data, n, cPBuf));
    cX.status = BARE;

    cudaStatCheck(  cudaMemcpy(x, cX.data, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
}

void cudaTridiagonalSolverSystem::setX(const std::complex<double>* x, bool virt, VectorState state) {
    if (state == EMPTY)
        throw std::runtime_error("cudaTridiagonalSolverSystem::setX : Passed state must be BARE or OPERATED.");

    CudaVector& myX = virt ? cXV : cX;
    cudaStatCheck(  cudaMemcpy(myX.data, x, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    myX.status = state;
}

void cudaTridiagonalSolverSystem::gatherX(std::complex<double>* x, bool virt) {
    CudaVector& myX = virt ? cXV : cX;
    if (myX.status != BARE)
        throw std::runtime_error("cudaTridiagonalSolverSystem::gatherX : Stored state is not BARE.");

    cudaStatCheck(  cudaMemcpy(x, myX.data, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
}

void cudaTridiagonalSolverSystem::gatherRHS(std::complex<double>* x, bool virt) {
    CudaVector& myX = virt ? cXV : cX;
    if (myX.status != OPERATED)
        throw std::runtime_error("cudaTridiagonalSolverSystem::gatherRHS : Stored state is not OPERATED.");

    cudaStatCheck(  cudaMemcpy(x, myX.data, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
}

void cudaTridiagonalSolverSystem::rhsProduct(const std::complex<double>* D, bool destVirt, bool sourceVirt) {
    CudaVector& sourceX = sourceVirt ? cXV : cX;
    if (sourceX.status != BARE)
        throw std::runtime_error("cudaTridiagonalSolverSystem::rhsProduct : Source state is not BARE.");
    CudaVector& destX = destVirt ? cXV : cX;

    cudaStatCheck(  cudaMemcpy(rhsTemp, D, n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    cudaStatCheck(  cublasZcopy_v2(cbHandle, n, rhsTemp, 1, cRHSMat+1, 3));

    cuDoubleComplex one = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex zero = make_cuDoubleComplex(0.0, 0.0);
    // gbmv cannot work in-place, so if the source and dest are the same use a temporary buffer
    if(destVirt == sourceVirt){
        for(int i = 0; i < nrhs; i++){
            cudaStatCheck(  cublasZgbmv_v2(cbHandle, CUBLAS_OP_N, n, n, 1, 1, &one, cRHSMat, 3, sourceX.data+(i*n), 1, &zero, rhsTemp, 1));
            cudaStatCheck(  cublasZcopy_v2(cbHandle, n, rhsTemp, 1, destX.data+(i*n), 1));
        }
    }
    else{
        for(int i = 0; i < nrhs; i++)
            cudaStatCheck(  cublasZgbmv_v2(cbHandle, CUBLAS_OP_N, n, n, 1, 1, &one, cRHSMat, 3, sourceX.data+(i*n), 1, &zero, destX.data+(i*n), 1));
    }

    // enforce BCs if set
    if(bdyRHSL != nullptr)
        cudaStatCheck(  cublasZcopy_v2(cbHandle, nrhs, bdyRHSL, 1, destX.data, n));
    if(bdyRHSR != nullptr)
        cudaStatCheck(  cublasZcopy_v2(cbHandle, nrhs, bdyRHSR, 1, destX.data+(n-1), n));

    destX.status = OPERATED;
}

void cudaTridiagonalSolverSystem::solve(const std::complex<double> *D, bool destVirt, bool sourceVirt){
    if (!lhsOffdiagDefined)
        throw std::runtime_error("cudaTridiagonalSolverSystem::solve : LHS off-diagonal elements are not defined.");

    CudaVector& sourceX = sourceVirt ? cXV : cX;
    if (sourceX.status != OPERATED)
        throw std::runtime_error("cudaTridiagonalSolverSystem::solve : Source state is not OPERATED.");
    CudaVector& destX = destVirt ? cXV : cX;

    cudaStatCheck(  cudaMemcpy(cD, D, n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    if (destVirt == sourceVirt)
        cudaStatCheck(  cusparseZgtsv2(csHandle, n, nrhs, cDL, cD, cDU, sourceX.data, n, cPBuf));
    else{
        if(tempState == nullptr)
            cudaStatCheck(  cudaMalloc((void**)&tempState, n * nrhs * sizeof(cuDoubleComplex)));
        cudaStatCheck(  cudaMemcpy(tempState, sourceX.data, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice));
        cudaStatCheck(  cusparseZgtsv2(csHandle, n, nrhs, cDL, cD, cDU, tempState, n, cPBuf));
        cudaStatCheck(  cudaMemcpy(destX.data, tempState, n * nrhs * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice));
    }

    destX.status = BARE;
}

void cudaTridiagonalSolverSystem::setBdyCond(std::complex<double> DOv, const std::complex<double>* RHSv, Side side){
    cuDoubleComplex* myDO = side == LHS ? cDU : cDL+(n-1);
    if (side == LHS && bdyRHSL == nullptr)
        cudaStatCheck(  cudaMalloc((void**)&bdyRHSL, sizeof(cuDoubleComplex)*nrhs));
    if (side == RHS && bdyRHSR == nullptr)
        cudaStatCheck(  cudaMalloc((void**)&bdyRHSR, sizeof(cuDoubleComplex)*nrhs));
    cuDoubleComplex* myRHS = side == LHS ? bdyRHSL : bdyRHSR;
    
    cudaStatCheck(  cudaMemcpy(myDO, &DOv, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
    cudaStatCheck(  cudaMemcpy(myRHS, RHSv, sizeof(cuDoubleComplex)*nrhs, cudaMemcpyHostToDevice));
}

void cudaTridiagonalSolverSystem::resetBdyCond(Side side){
    if (side == LHS && bdyRHSL != nullptr) {
        cudaStatCheck(  cudaFree(bdyRHSL));
        bdyRHSL = nullptr;
    } else if (side == RHS && bdyRHSR != nullptr) {
        cudaStatCheck(  cudaFree(bdyRHSR));
        bdyRHSR = nullptr;
    }
}

void cudaTridiagonalSolverSystem::calcRawRho(const double* weights, double* rho, bool virt) {
    CudaVector& myX = virt ? cXV : cX;
    if (myX.status != BARE)
        throw std::runtime_error("cudaTridiagonalSolverSystem::calcRawRho : Stored state is not BARE.");
    
    cudaStatCheck(  cudaMemcpy(cWeights, weights, nrhs * sizeof(double), cudaMemcpyHostToDevice));
    cudaStatCheck(  cudaDensity(cWeights, myX.data, cRho, n, nrhs));
    cudaStatCheck(  cudaMemcpy(rho, cRho, n * sizeof(double), cudaMemcpyDeviceToHost)); // TEMPORARY FOR TESTING
}