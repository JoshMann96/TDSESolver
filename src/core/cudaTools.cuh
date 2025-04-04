/**
 * \file cudaTools.cuh
 * \brief CUDA tools for computing norms and densities.
 */
#pragma once
#include <cuComplex.h>

/**
 * Compute the square of the norm of a vector via CUDA. All arrays are assumed to be on the device.
 * @param x (in) The vector.
 * @param result (out) The result.
 * @param n The number of elements in the vector.
 * @return cudaError code.
 */
cudaError_t cudaNormSquare(const cuDoubleComplex* x, double* result, int n);

/**
 * Compute the density of a vector via CUDA. All arrays are assumed to be on the device.
 * @param weights (in) The weights.
 * @param x (in) Wavefunctions.
 * @param result (out) Result.
 * @param n Number of gridpoints.
 * @param nrhs Number of right-hand sides.
 * @return cudaError code.
 */
cudaError_t cudaDensity(const double* weights, const cuDoubleComplex* x, double* result, int n, int nrhs);

/**
 * Compute the Hadamard product of a vector and a matrix via CUDA. All arrays are assumed to be on the device.
 * @param vec (in) The which will multiply the elemtns of x. Size must be n.
 * @param x (in/out) The vectors to be modified. It will be multiplied element-wise by vec. Size must be n*nrhs.
 * @param n The number of elements in each vector (length of vec and rows of x).
 * @param nrhs The number of right-hand sides (columns of x).
 * @return cudaError code.
 */
cudaError_t cudaHadamard(const double* vec, cuDoubleComplex* x, int n, int nrhs);