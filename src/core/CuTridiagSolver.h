#pragma once

#include <cusparse.h>
#include <cuda_runtime.h>
#include <iostream>
#include <complex>

class cuTridiagSolver {
private:
    int n, nrhs;
    cusparseHandle_t handle;
    cuDoubleComplex *cDL, *cD, *cDU, *cB, *cPBuf;
    bool collected;
public:
    cuTridiagSolver(int n, int nrhs);
    ~cuTridiagSolver();
    void solve(std::complex<double> *DL, std::complex<double> *D, std::complex<double> *DU, std::complex<double> *x);
};
