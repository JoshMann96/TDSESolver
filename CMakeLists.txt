cmake_minimum_required(VERSION 3.13)

#aocc environment variables
#set(ENV{LIBRARY_PATH} "/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/lib:/usr/lib32:$ENV{LIBRARY_PATH}")
#set(ENV{LD_LIBRARY_PATH} "/usr/lib64:/usr/lib/x86_64-linux-gnu:/usr/lib:/usr/lib32:$ENV{LD_LIBRARY_PATH}")
#set(ENV{PATH} "/home/jmann/aocc-compiler-4.1.0/bin:$ENV{PATH}")
#set(ENV{PATH} "/home/jmann/aocc-compiler-4.1.0/share/opt-viewer:$ENV{PATH}")
#set(ENV{LIBRARY_PATH} "/home/jmann/aocc-compiler-4.1.0/lib:/home/jmann/aocc-compiler-4.1.0/lib32:$ENV{LIBRARY_PATH}")
#set(ENV{LD_LIBRARY_PATH} "/home/jmann/aocc-compiler-4.1.0/ompd:/home/jmann/aocc-compiler-4.1.0/lib:/home/jmann/aocc-compiler-4.1.0/lib32:$ENV{LD_LIBRARY_PATH}")

#IF(DEFINED $ENV{C_INCLUDE_PATH})
#    set(ENV{C_INCLUDE_PATH} "$ENV{C_INCLUDE_PATH}:/home/jmann/aocc-compiler-4.1.0/include")
#ELSE()
#    set(ENV{C_INCLUDE_PATH} "/home/jmann/aocc-compiler-4.1.0/include")
#ENDIF()

#IF(DEFINED $ENV{CPLUS_INCLUDE_PATH})
#    set(ENV{CPLUS_INCLUDE_PATH} "$ENV{CPLUS_INCLUDE_PATH}:/home/jmann/aocc-compiler-4.1.0/include")
#ELSE()
#    set(ENV{CPLUS_INCLUDE_PATH} "/home/jmann/aocc-compiler-4.1.0/include")
#ENDIF()

#aocl environment variables
#set(ENV{AOCL_ROOT} "/home/jmann/aocl/4.1.0/aocc")
#set(ENV{C_INCLUDE_PATH} "/home/jmann/aocl/4.1.0/aocc/include:$ENV{C_INCLUDE_PATH}")
#set(ENV{CPLUS_INCLUDE_PATH} "/home/jmann/aocl/4.1.0/aocc/include:$ENV{CPLUS_INCLUDE_PATH}")
#set(ENV{LD_LIBRARY_PATH} "/home/jmann/aocl/4.1.0/aocc/lib:$ENV{LD_LIBRARY_PATH}")
#set(ENV{LIBRARY_PATH} "/home/jmann/aocl/4.1.0/aocc/lib:$ENV{LIBRARY_PATH}")

project(TDSESolveLinux)

#set compiler
set(CMAKE_C_COMPILER "/home/jmann/aocc-compiler-4.1.0/bin/clang")
set(CMAKE_CXX_COMPILER "/home/jmann/aocc-compiler-4.1.0/bin/clang++")
set(CMAKE_Fortran_COMPILER "/home/jmann/aocc-compiler-4.1.0/bin/flang")

message("C   Compiler: ${CMAKE_C_COMPILER}")
message("CXX Compiler: ${CMAKE_CXX_COMPILER}")
message("FTN Compiler: ${CMAKE_Fortran_COMPILER}")

# Compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -std=c++1y -O3 -march=znver2 -mtune=native -DNDEBUG -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp")

# OpenMP Configuration (USING AOCL BELOW)
#find_package(OpenMP REQUIRED)

# Boost Configuration
set(BOOST_INC "/usr/include/boost")
include_directories(${BOOST_INC})
set(BOOST_LINK_FLAGS "-lboost_system -lboost_filesystem")

#AOCC config
set(AOCC_ROOT "/home/jmann/aocc")
set(AOCC_INC "${AOCC_ROOT}/include")
set(AOCC_LIB "${AOCC_ROOT}/lib")
set(AOCC_LIBS "-L${AOCL_LIB} -lomp -lflang")
include_directories(${AOCC_INC})
#find_library(FLANG_LIB      NAMES flang     HINTS ${AOCC_LIB})
#find_library(OMP_LIB        NAMES omp      HINTS ${AOCC_LIB})

#message("Using flang: ${FLANG_LIB}")
#message("Using omp: ${OMP_LIB}")

#AOCL config
set(AOCL_ROOT "/home/jmann/aocl/4.1.0/aocc")
set(AOCL_INC "${AOCL_ROOT}/include")
set(AOCL_LIB "${AOCL_ROOT}/lib")
include_directories(${AOCL_INC})
set(AOCL_LIBS "-L${AOCL_LIB} -lblis-mt -lflame -lfftw3_omp -lfftw3 -lalm -lm")
#find_library(BLIS_LIB       NAMES libblis-mt  HINTS ${AOCL_LIB})
#find_library(PTHREAD_LIB    NAMES pthread   HINTS ${AOCL_LIB})
#find_library(FFTW3_OMP_LIB  NAMES fftw3_omp     HINTS ${AOCL_LIB})
#find_library(FFTW3_LIB      NAMES fftw3         HINTS ${AOCL_LIB})
#find_library(FLAME_LIB      NAMES flame         HINTS ${AOCL_LIB})
#find_library(ALM_LIB        NAMES alm           HINTS ${AOCL_LIB})
#find_library(M_LIB          NAMES m)

#message("Using blis: ${BLIS_LIB}")
#message("Using FFTW3_OMP: ${FFTW3_OMP_LIB}")
#message("Using FFTW3: ${FFTW3_LIB}")
#message("Using libFLAME: ${FLAME_LIB}")
#message("Using libALM: ${ALM_LIB}")
#message("Using libM: ${M_LIB}")

#set(AOCL_LINK_FLAGS "-lm -lblis -lpthread -lfftw3 -lstdc++ -lflame")

# FFTW3 Config
#set(FFTW_INC "/usr/include")
#set(FFTW_LIB "/usr/lib")
#include_directories(${FFTW_INC})
#link_directories(${FFTW_LIB})
#set(FFTW_LINK_FLAGS "-lfftw3 -lm")

# Intel MKL Configuration
#set(MKL_INC "/usr/include/mkl")
#set(MKL_LIB "/usr/lib/x86_64-linux-gnu/mkl/lib")
#include_directories(${MKL_INC})
#link_directories(${MKL_LIB})
#set(MKL_LINK_FLAGS "-lmkl_intel_thread -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl")

#MPI Configuration
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Add source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB IMPLEMENTS "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Add the executable
add_executable(TDSESolveLinux ${SOURCES} ${IMPLEMENTS} ${HEADERS})

# Include directories
#target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${BOOST_INC_ROOT} ${MKL_INC_ROOT} ${MPI_INCLUDE_PATH})
target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${BOOST_INC_ROOT} ${MPI_INCLUDE_PATH} ${AOCL_INC})

# Link libraries
#target_link_libraries(TDSESolveLinux PRIVATE OpenMP::OpenMP_CXX MPI::MPI_CXX ${MKL_LIBS} ${BOOST_LIBS})
target_link_libraries(TDSESolveLinux PRIVATE MPI::MPI_CXX ${AOCC_LIBS} ${BOOST_LINK_FLAGS} ${AOCL_LIBS})