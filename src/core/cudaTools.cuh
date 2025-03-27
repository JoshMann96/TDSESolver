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