#include "cudaTools.cuh"
#include <cuda_runtime.h>
#include <iostream>

__global__ void _cudaNormSquare(const cuDoubleComplex* x, double* result, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n)
        result[idx] = x[idx].x * x[idx].x + x[idx].y * x[idx].y;
}

cudaError_t cudaNormSquare(const cuDoubleComplex* x, double* result, size_t n) {
    size_t threadsPerBlock = 256;
    size_t blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
    _cudaNormSquare<<<blocksPerGrid, threadsPerBlock>>>(x, result, n);
    return cudaGetLastError();
}


__global__ void _cudaDensity(const double* weights, const cuDoubleComplex* x, double* result, size_t n, size_t nrhs) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n){
        for(size_t i = 0; i < nrhs; i++){
            result[idx] += weights[i] * (x[i*n + idx].x * x[i*n + idx].x + x[i*n + idx].y * x[i*n + idx].y);
        }
    }
}

cudaError_t cudaDensity(const double* weights, const cuDoubleComplex* x, double* result, size_t n, size_t nrhs) {
    cudaMemset(result, 0, n * sizeof(double));
    size_t threadsPerBlock = 256;
    size_t blocksPerGrid = (n*nrhs + threadsPerBlock - 1) / threadsPerBlock;
    _cudaDensity<<<blocksPerGrid, threadsPerBlock>>>(weights, x, result, n, nrhs);
    return cudaGetLastError();
}

__global__ void _cudaHadamard(const double* vec, cuDoubleComplex* x, size_t n, size_t nrhs) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        for (size_t i = 0; i < nrhs; i++) {
            size_t index = i * n + idx;
            x[index].x *= vec[idx];
            x[index].y *= vec[idx];
        }
    }
}

cudaError_t cudaHadamard(const double* vec, cuDoubleComplex* x, size_t n, size_t nrhs) {
    size_t threadsPerBlock = 256;
    size_t blocksPerGrid = (n * nrhs + threadsPerBlock - 1) / threadsPerBlock;

    _cudaHadamard<<<blocksPerGrid, threadsPerBlock>>>(vec, x, n, nrhs);
    return cudaGetLastError();
}