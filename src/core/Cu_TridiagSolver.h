/**
 * @file CuTridiagSolver.h
 * @brief CUDA tridiagonal solver system
 * @details This file contains the definition of the cudaTridiagonalSolverSystem class, which is used to solve tridiagonal systems of the form \f$A X = B X_0\f$ on the GPU.
 * The class minimizes communication between the CPU and GPU, allowing for efficient solving of tridiagonal systems.
 * It supports virtual states which preserve the original system, useful for accurate iterations for nonlinear systems.
 */
#pragma once

#include <cusparse.h>
#include <cublas.h>
#include <cuda_runtime.h>
#include <iostream>
#include <complex>

/**
 * CUDA tridiagonal solver system
 * 
 * This class is used to solve tridiagonal systems of the form \f$A X = B X_0\f$ on the GPU while minimizing communication between the CPU and GPU.
 * The present state is stored on the device and its status may be EMPTY (uninitialized), BARE (initialized), or OPERATED (operated on by a matrix).
 * Supports virtual states which preserve the original system, useful for accurate iterations for nonlinear systems.
*/
class cudaTridiagonalSolverSystem {
public:
    /// The state of a wavefunction vector stored on the GPU.
    enum VectorState {
        /// The vector is not initialized or has been reset.
        EMPTY, 
        /// The vector is initialized and contains a valid wavefunction.
        BARE, 
        /// The vector has been operated on and contains the result of a linear product (RHS of the equation).
        OPERATED
    };
private:
    struct CudaVector {
        cuDoubleComplex *data = nullptr;
        VectorState status = EMPTY;
    } cX, cXV; // solution vector and virtual solution vector

    size_t n, nrhs;
    cusparseHandle_t csHandle;
    cublasHandle_t cbHandle;
    cuDoubleComplex *cDL, *cD, *cDU; // LHS matrix definition
    cuDoubleComplex *_cX, *_cXV, *cPBuf, *tempState = nullptr; // solution vector (data) and workspace
    double *cRho, *cCur, *cWeights, *cVec=nullptr; // density

    cuDoubleComplex *cRHSMat, *rhsTemp; // RHS matrix values

    bool lhsOffdiagDefined, rhsOffdiagDefined;

    cuDoubleComplex *bdyRHSL = nullptr, *bdyRHSR = nullptr; // boundary condition values for the LHS and RHS

    static void cudaStatCheck(cudaError_t status){
        #ifdef USE_CUDA
        if (status != cudaSuccess)
            throw std::runtime_error("cudaTridiagonalSolverSystem : cuda error with status: " + std::string(cudaGetErrorString(status)));
        #endif
    }
    static void cudaStatCheck(cublasStatus_t status){
        if (status != CUBLAS_STATUS_SUCCESS)
            throw std::runtime_error("cudaTridiagonalSolverSystem : cublas error with status: " + std::string(cublasGetStatusString(status)));
    }
    static void cudaStatCheck(cusparseStatus_t status){
        if (status != CUSPARSE_STATUS_SUCCESS)
            throw std::runtime_error("cudaTridiagonalSolverSystem : cusparse error with status: " + std::string(cusparseGetErrorString(status)));
    }

public:
    enum Side {LHS, RHS};

    /**
     * Initializes a new cudaTridiagonalSolverSystem.
     * @param n number of gridpoints in the system.
     * @param nrhs number of right hand sides (Kohn-Sham orbitals).
    */
    cudaTridiagonalSolverSystem(size_t n, size_t nrhs);

    ~cudaTridiagonalSolverSystem();

    /**
     * Solves a tridiagonal system of the form \f$A X = X_0\f$. This is intended for solving without maintaining the system on the GPU.
     * The resulting state is stored in the regular state.
     * @param DL (in) lower diagonal of \f$A\f$, \a n-1 elements.
     * @param D (in) main diagonal of \f$A\f$, \a n elements.
     * @param DU (in) upper diagonal of \f$A\f$, \a n elements.
     * @param x (in/out) solution vector, contains \f$X_0\f$ and is overwritten with \f$X\f$, \a n * \a nrhs elements.
     */
    void solve(const std::complex<double> *DL, const std::complex<double> *D, const std::complex<double> *DU, std::complex<double> *x);
    
    /**
     * Sets the internal off-diagonal elements of the tridiagonal system. Side determines whether it is \f$A\f$ or \f$B\f$ being set.
     * @param DL (in) lower diagonal of \f$A\f$ or \f$B\f$, \a n-1 elements.
     * @param DU (in) upper diagonal of \f$A\f$ or \f$B\f$, \a n-1 elements.
     * @param side whether to set the system's LHS (\f$A\f$) or RHS (\f$B\f$).
     */
    void setOffDiag(const std::complex<double>* DL, const std::complex<double>* DU, Side side);
    
    /**
     * Solves a tridiagonal system of the form \f$A X = X_0\f$. This is intended for solving with the off-diagonal components of A already set on the GPU.
     * The resulting state is stored in the regular state.
     * @param D (in) main diagonal of \f$A\f$, \a n elements.
     * @param x (in/out) solution vector, contains \f$X_0\f$ and is overwritten with \f$X\f$, \a n*nrhs elements.
     * @throws std::runtime_error if the off-diagonal elements of \f$A\f$ have not been set.
     */
    void solve(const std::complex<double> *D, std::complex<double> *x);
    
    /**
     * Sets the state vector on the GPU.
     * @param x (in) state vector, \a n*nrhs elements.
     * @param virt true sets the virtual state, false sets the regular state.
     * @param state vector state to set, BARE, or OPERATED.
     * @throws std::runtime_error if the vector state is not BARE or OPERATED.
     */
    void setX(const std::complex<double>* x, bool virt = false, VectorState state = BARE);
    
    /**
     * Gathers the state vector from the GPU to the CPU.
     * @param x (out) state vector, \a n*nrhs elements.
     * @param virt true gathers the virtual state, false gathers the regular state.
     * @throws std::runtime_error if the vector state is not BARE.
     */
    void gatherX(std::complex<double>* x, bool virt = false);
    
    /**
     * Gathers the RHS vector from the GPU to the CPU.
     * @param x (out) RHS vector, \a n*nrhs elements.
     * @param virt true gathers the virtual RHS, false gathers the regular RHS.
     * @throws std::runtime_error if the vector state is not OPERATED.
     */
    void gatherRHS(std::complex<double>* x, bool virt = false);
    
    /**
     * Computes the product \f$B X_0\f$.
     * @param D (in) main diagonal of \f$B\f$, n elements.
     * @param destVirt true stores the result in the virtual state, false stores the result in the regular state.
     * @param sourceVirt true uses the virtual state as \f$X_0\f$, false uses the regular state as \f$X_0\f$.
     * @throws std::runtime_error if the source vector state is not BARE.
     * @throws std::runtime_error if the off-diagonal elements of \f$B\f$ have not been set.
     */
    void rhsProduct(const std::complex<double>* D, bool destVirt = false, bool sourceVirt = false); // computes the product B X0 and stores it in the RHS, destVirt determines whether the result (B X0) is stored in the regular or virtual state, and sourceVirt determines whether the source (X0) is taken from the regular or virtual state
    
    /**
     * Sets the boundary conditions for the left or right side of the \a physical system.
     * The diagonal component of the boundary condition must be passed by the call to solve.
     * @param DOv off-diagonal element of \f$A\f$.
     * @param RHSv the right-hand side value of the operated vector.
     * @param side whether to set the LHS or RHS boundary conditions.
     */
    void setBdyCond(std::complex<double> DOv, const std::complex<double>* RHSv, Side side); // sets the boundary conditions for the LHS or RHS of the physical system
    
    /**
     * Resets the boundary conditions for the left or right side of the \a physical system.
     * This is necessary to remove boundary conditions that are no longer needed.
     * @param side whether to reset the LHS or RHS boundary conditions.
     */
    void resetBdyCond(Side side); // resets the boundary conditions for the LHS or RHS of the physical system

    /**
     * Solves a tridiagonal system of the form \f$A X = B X_0\f$. 
     * The off-diagonal elements of \f$A\f$ must have been set previously, and rhsProduct must have been called to compute \f$B X_0\f$.
     * @param D (in) main diagonal of \f$A\f$, \a n elements.
     * @param destVirt true stores the result in the virtual state, false stores the result in the regular state.
     * @param sourceVirt true uses the virtual state as \f$X_0\f$, false uses the regular state as \f$X_0\f$.
     * @throws std::runtime_error if the source vector state is not OPERATED (did yuo call rhsProduct?).
     * @throws std::runtime_error if the off-diagonal elements of \f$A\f$ have not been set.
     */
    void solve(const std::complex<double> *D, bool destVirt = false, bool sourceVirt = false);
    
    /**
     * Calculates the density according to weights on the GPU and returns it to the CPU.
     * This reduces the amount of communication overhead between the CPU and GPU.
     * @param weights (in) weights of the states, \a nrhs elements.
     * @param rho (out) density, n elements.
     * @param virt true calculates the density from the virtual state, false calculates the density from the regular state.
     */
    void calcRawRho(const double* weights, double* rho, bool virt = false);

    /**
     * Calculates the current density according to weights on the GPU and returns it to the CPU.
     * This reduces the amount of communication overhead between the CPU and GPU.
     * @param weights (in) weights of the states, \a nrhs elements.
     * @param cur (out) current density, n elements.
     * @param prefactor includes the prefactor for the current density: \f$ \hbar/(m \Delta x) \f$
     * @param virt true calculates the current density from the virtual state, false calculates the current density from the regular state.
     */
    void calcRawCur(const double* weights, double* cur, double prefactor, bool virt = false);

    /**
     * Calculates the Hadamard product of the current state vector with a given vector.
     * This is an element-wise multiplication of the two vectors, repeating for each substate vector.
     * @param vec (in) The vector to multiply with, must be of size \a n.
     * @param virt If true, performs the operation and stores the result on the virtual state; otherwise, operates and stores on the regular state.
     */
    void vectorHadamardProduct(const double* vec, bool virt = false);

    // TODO:
    // Permit gathering only part of the state to minimize communication further
    // This will need some modification of the Measurers classes
    // Add method to Measurer requiring each to declare which points are needed.
    // Allow measurers to return some sort of Kernel, so that whatever they need to calculate can be done efficiently on the GPU.
};