#pragma once

#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <climits>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <chrono>
#define NOMINMAX
#include <iomanip>
#include <future>
#include <cstring>
#include <string>
#include <stdexcept>

#include <omp.h>

#ifdef USE_RESTRICT
#else
#define __restrict
#endif

#include <boost/unordered_map.hpp>
#include <boost/assign/list_of.hpp>

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/tools/rational.hpp>

#include <fftw3.h>

inline std::mutex mtx;

inline void* thread_safe_malloc(size_t size){
    mtx.lock();
    //auto res = fftw_malloc(size);
    void* res;
    posix_memalign(&res, 128, size);
    if(res == NULL){
        std::cerr << "Memory allocation failed" << std::endl;
        exit(1);
    }
    mtx.unlock();
    return res;
}

inline void thread_safe_free(void* ptr){
    mtx.lock();
    //fftw_free(ptr);
    free(ptr);
    mtx.unlock();
}

#define sq_malloc thread_safe_malloc
#define sq_free thread_safe_free