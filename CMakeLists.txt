cmake_minimum_required(VERSION 3.19)

project(TDSESolveLinux)

set(CMAKE_CXX_STANDARD 17)

message("C   Compiler: ${CMAKE_C_COMPILER}")
message("CXX Compiler: ${CMAKE_CXX_COMPILER}")
message("FTN Compiler: ${CMAKE_Fortran_COMPILER}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


### MODIFY FOR YOUR SYSTEM HERE ###

#OpenBLAS config
set(OBLAS_ROOT "/home/jmann")
set(FFTW_ROOT "/home/jmann")

# Boost config
set(BOOST_INC "/usr/include/boost")

#PYBIND11
set(Python_ROOT_DIR ${CMAKE_SOURCE_DIR}/.venv)
set(pybind11_DIR ${CMAKE_SOURCE_DIR}/.venv/lib/python3.11/site-packages/pybind11/share/cmake/pybind11/)

###################################


#pybind11
set(PYTHON_EXECUTABLE ${Python_ROOT_DIR}/bin/python3)
find_package(pybind11 REQUIRED)

#OMP
find_package(OpenMP REQUIRED)

#MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#BOOST
include_directories(${BOOST_INC})
set(BOOST_LINK_FLAGS "-lboost_system -lboost_filesystem")

#FFTW config
set(FFTW_INC "${FFTW_ROOT}/include")
set(FFTW_LIB "${FFTW_ROOT}/lib")
include_directories(${FFTW_INC})
set(FFTW_LIBS "-L${FFTW_LIB} -lfftw3_omp -lfftw3 -lm")

#OpenBLAS config
set(OBLAS_INC "${OBLAS_ROOT}/include")
set(OBLAS_LIB "${OBLAS_ROOT}/lib")
find_library(OBLAS_LIBS
  NAMES openblas libopenblas
  REQUIRED
  HINTS ${OBLAS_LIB}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lstdc++ -O3 -march=znver2 -mtune=native -fPIC -fopenmp")

add_subdirectory(src)
