#include "cudaTools.cuh"
#include <cuda_runtime.h>

__global__ void _cudaNormSquare(const cuDoubleComplex* x, double* result, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n)
        result[idx] = x[idx].x * x[idx].x + x[idx].y * x[idx].y;
}

cudaError_t cudaNormSquare(const cuDoubleComplex* x, double* result, int n) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
    _cudaNormSquare<<<blocksPerGrid, threadsPerBlock>>>(x, result, n);
    return cudaGetLastError();
}


__global__ void _cudaDensity(const double* weights, const cuDoubleComplex* x, double* result, int n, int nrhs) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n){
        for(int i = 0; i < nrhs; i++){
            result[idx] += weights[i] * (x[i*n + idx].x * x[i*n + idx].x + x[i*n + idx].y * x[i*n + idx].y);
        }
    }
}

cudaError_t cudaDensity(const double* weights, const cuDoubleComplex* x, double* result, int n, int nrhs) {
    cudaMemset(result, 0, n * sizeof(double));
    int threadsPerBlock = 256;
    int blocksPerGrid = (n*nrhs + threadsPerBlock - 1) / threadsPerBlock;
    _cudaDensity<<<blocksPerGrid, threadsPerBlock>>>(weights, x, result, n, nrhs);
    return cudaGetLastError();
}

__global__ void _cudaHadamard(const double* vec, cuDoubleComplex* x, int n, int nrhs) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        for (int i = 0; i < nrhs; i++) {
            int index = i * n + idx;
            x[index].x *= vec[idx];
            x[index].y *= vec[idx];
        }
    }
}

cudaError_t cudaHadamard(const double* vec, cuDoubleComplex* x, int n, int nrhs) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (n * nrhs + threadsPerBlock - 1) / threadsPerBlock;
    _cudaHadamard<<<blocksPerGrid, threadsPerBlock>>>(vec, x, n, nrhs);
    return cudaGetLastError();
}