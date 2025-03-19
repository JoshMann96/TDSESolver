#pragma once

#include <cusparse.h>
#include <cublas.h>
#include <cuda_runtime.h>
#include <iostream>
#include <complex>

// manages the solution of a tridiagonal system of equations on the GPU
// aims to minimize communication between CPU and GPU
class cudaTridiagonalSolverSystem {
public:
    enum VectorState {EMPTY, BARE, OPERATED}; // EMPTY: vector not initialized, BARE: state is set (it stores a wavefunction), OPERATED: state has been operated on (it stores the RHS of the linear equation)
private:
    struct CudaVector {
        cuDoubleComplex *data = nullptr;
        VectorState status = EMPTY;
    } cX, cXV; // solution vector and virtual solution vector

    int n, nrhs;
    cusparseHandle_t csHandle;
    cublasHandle_t cbHandle;
    cuDoubleComplex *cDL, *cD, *cDU; // LHS matrix definition
    cuDoubleComplex *_cX, *_cXV, *cPBuf, *tempState = nullptr; // solution vector (data) and workspace
    double *cRho, *cWeights; // density

    cuDoubleComplex *cRHSMat, *rhsTemp; // RHS matrix values

    bool lhsOffdiagDefined, rhsOffdiagDefined;

    cuDoubleComplex *bdyRHSL = nullptr, *bdyRHSR = nullptr; // boundary condition values for the LHS and RHS

    void cudaStatCheck(cudaError_t status){
        if (status != cudaSuccess)
            throw std::runtime_error("cudaTridiagonalSolverSystem : cuda error with status: " + std::string(cudaGetErrorString(status)));
    }
    void cudaStatCheck(cublasStatus_t status){
        if (status != CUBLAS_STATUS_SUCCESS)
            throw std::runtime_error("cudaTridiagonalSolverSystem : cublas error with status: " + std::string(cublasGetStatusString(status)));
    }
    void cudaStatCheck(cusparseStatus_t status){
        if (status != CUSPARSE_STATUS_SUCCESS)
            throw std::runtime_error("cudaTridiagonalSolverSystem : cusparse error with status: " + std::string(cusparseGetErrorString(status)));
    }

public:
    enum Side {LHS, RHS};

    cudaTridiagonalSolverSystem(int n, int nrhs);
    ~cudaTridiagonalSolverSystem();
    
    // descriptions assume system is A X = B X0
    // virtual states which preserve the original system are also supported

    void solve(const std::complex<double> *DL, const std::complex<double> *D, const std::complex<double> *DU, std::complex<double> *x); // solve a tridiagonal system, x contains X0 and is overwritten with X
    void setOffDiag(const std::complex<double>* DL, const std::complex<double>* DU, Side side); // sets the internal off-diagonal elements of the tridiagonal system, side determines whether it is A or B being set
    void solve(const std::complex<double> *D, std::complex<double> *x); // solve a tridiagonal system with the LHS diagonals already set, x contains X0 and is overwritten with X
    void setX(const std::complex<double>* x, bool virt = false, VectorState state = BARE); // sets the solution vector on the GPU, virt = true sets the virtual solution
    void gatherX(std::complex<double>* x, bool virt = false); // gathers the solution vector from the GPU to the CPU, virt = true returns the virtual solution
    void gatherRHS(std::complex<double>* x, bool virt = false); // gathers the RHS vector from the GPU to the CPU, virt = true returns the virtual RHS
    void rhsProduct(const std::complex<double>* D, bool destVirt = false, bool sourceVirt = false); // computes the product B X0 and stores it in the RHS, destVirt determines whether the result (B X0) is stored in the regular or virtual state, and sourceVirt determines whether the source (X0) is taken from the regular or virtual state
    void setBdyCond(std::complex<double> DOv, const std::complex<double>* RHSv, Side side); // sets the boundary conditions for the LHS or RHS of the physical system
    void resetBdyCond(Side side); // resets the boundary conditions for the LHS or RHS of the physical system

    void solve(const std::complex<double> *D, bool destVirt = false, bool sourceVirt = false); // solves a tridiagonal system with the LHS diagonals already set, destVirt determines whether the solution X is stored in the regular or virtual state, and sourceVirt determines whether the source (B X0) is taken from the regular or virtual state
    void calcRawRho(const double* weights, double* rho, bool virt = false); // calculates the density of the state on the GPU and returns it to the CPU, virt determines whether the density is calculated from the regular or virtual state
};
