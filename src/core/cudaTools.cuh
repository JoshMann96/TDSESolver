#pragma once
#include <cuComplex.h>

cudaError_t cudaNormSquare(const cuDoubleComplex* x, double* result, int n);

cudaError_t cudaDensity(const double* weights, const cuDoubleComplex* x, double* result, int n, int nrhs);