cmake_minimum_required(VERSION 3.19)

project(TDSESolver)

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

# Boost config
set(BOOST_INC "/usr/include/boost")

# CPU type
set(CPU_MARCH "znver2")

###################################

#PYBIND11
set(Python_ROOT_DIR ${CMAKE_SOURCE_DIR}/.venv)
set(pybind11_DIR ${CMAKE_SOURCE_DIR}/.venv/lib/python3.$ENV{PYTHON_SUBVERSION_NUMBER}/site-packages/pybind11/share/cmake/pybind11/)

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
find_package(PkgConfig REQUIRED)     
pkg_check_modules(FFTW IMPORTED_TARGET REQUIRED fftw3)
          
if(NOT FFTW_ROOT AND DEFINED ENV{FFTWDIR})
    set(FFTW_ROOT $ENV{FFTWDIR})
endif()
          
find_library(
    FFTW_DOUBLE_OMP_LIB
    NAMES "fftw3_omp"
    PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
)
          
if (FFTW_DOUBLE_OMP_LIB)
    set(FFTW_DOUBLE_OMP_LIB_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_OMP_LIB})
    add_library(FFTW::DoubleOMP INTERFACE IMPORTED)
    set_target_properties(FFTW::DoubleOMP
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_OMP_LIB}"
    )
else()
    set(FFTW_DOUBLE_OMP_LIB_FOUND FALSE)
endif()

#set(FFTW_INC "${FFTW_ROOT}/include")
#set(FFTW_LIB "${FFTW_ROOT}/lib")
#include_directories(${FFTW_INC})
#set(FFTW_LIBS "-L${FFTW_LIB} -lfftw3_omp -lfftw3 -lm")

#OpenBLAS config
set(OBLAS_INC "${OBLAS_ROOT}/include")
set(OBLAS_LIB "${OBLAS_ROOT}/lib")
find_library(OBLAS_LIBS
  NAMES openblas libopenblas
  REQUIRED
  HINTS ${OBLAS_LIB}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lstdc++ -O3 -march=${CPU_MARCH} -mtune=native -fPIC -fopenmp")

add_subdirectory(src)
