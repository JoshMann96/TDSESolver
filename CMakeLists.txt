cmake_minimum_required(VERSION 3.13)

project(TDSESolveLinux)

#set compiler
set(CMAKE_C_COMPILER "/home/jmann/bin/clang")
set(CMAKE_CXX_COMPILER "/home/jmann/bin/clang++")
set(CMAKE_Fortran_COMPILER "/home/jmann/bin/flang")

message("C   Compiler: ${CMAKE_C_COMPILER}")
message("CXX Compiler: ${CMAKE_CXX_COMPILER}")
message("FTN Compiler: ${CMAKE_Fortran_COMPILER}")

# Compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=znver2 -mtune=native -DNDEBUG -fopenmp")

# Boost Configuration
set(BOOST_INC "/usr/include/boost")
include_directories(${BOOST_INC})
set(BOOST_LINK_FLAGS "-lboost_system -lboost_filesystem")

#AOCC config
set(AOCC_ROOT "/home/jmann")
set(AOCC_INC "${AOCC_ROOT}/include")
set(AOCC_LIB "${AOCC_ROOT}/lib")
set(AOCC_LIBS "-L${AOCC_LIB} -lomp -lflang")
include_directories(${AOCC_INC})

#AOCL config
set(AOCL_ROOT "/home/jmann")
set(AOCL_INC "${AOCL_ROOT}/include")
set(AOCL_LIB "${AOCL_ROOT}/lib")
set(BLIS_INC "${AOCL_INC}/blis")
include_directories(${AOCL_INC})
#include_directories(${BLIS_INC})
#set(AOCL_LIBS "-L${AOCL_LIB} -lflame -lfftw3_omp -lfftw3 -lblis-mt  -lalm -lm -lstdc++ -lpthread")
set(AOCL_LIBS "-L${AOCL_LIB} -lfftw3_omp -lfftw3 -lalm -lm -lstdc++ -lpthread")

#OpenBLAS config
set(OBLAS_ROOT "/home/jmann")
set(OBLAS_INC "${OBLAS_ROOT}/include")
set(OBLAS_LIB "${OBLAS_ROOT}/lib")
set(OBLAS_LIBS "-L${OBLAS_LIB} -lopenblas -lflang -lpthread")

#MPI Configuration
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Add source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB IMPLEMENTS "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Add the executable
add_executable(TDSESolveLinux ${SOURCES} ${IMPLEMENTS} ${HEADERS})

#add_executable(Timings
#    "${CMAKE_CURRENT_SOURCE_DIR}/../lapack_timing/main.cpp" 
#    "${CMAKE_CURRENT_SOURCE_DIR}/KineticOperator.cpp" 
#    "${CMAKE_CURRENT_SOURCE_DIR}/MathTools.cpp"
#    "${CMAKE_CURRENT_SOURCE_DIR}/MathTools.tpp"
#    "${CMAKE_CURRENT_SOURCE_DIR}/KineticOperator.h"
#    "${CMAKE_CURRENT_SOURCE_DIR}/CommonHeader.h"
#    "${CMAKE_CURRENT_SOURCE_DIR}/MathTools.h"
#    "${CMAKE_CURRENT_SOURCE_DIR}/PhysCon.h")

# Include directories
#target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${BOOST_INC_ROOT} ${MKL_INC_ROOT} ${MPI_INCLUDE_PATH})
target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${MPI_INCLUDE_PATH} ${AOCL_INC} ${BLIS_INC} ${BOOST_INC_ROOT} ${OBLAS_INC})
#target_include_directories(Timings PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${MPI_INCLUDE_PATH} ${AOCL_INC} ${BLIS_INC} ${BOOST_INC_ROOT} ${OBLAS_INC})

# Link libraries
#target_link_libraries(TDSESolveLinux PRIVATE OpenMP::OpenMP_CXX MPI::MPI_CXX ${MKL_LIBS} ${BOOST_LIBS})
target_link_libraries(TDSESolveLinux PRIVATE MPI::MPI_CXX ${AOCC_LIBS} ${BOOST_LINK_FLAGS} ${AOCL_LIBS} ${OBLAS_LIBS})
#target_link_libraries(Timings PRIVATE MPI::MPI_CXX ${AOCC_LIBS} ${BOOST_LINK_FLAGS} ${AOCL_LIBS} ${OBLAS_LIBS})