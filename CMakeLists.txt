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

#OpenBLAS config
set(OBLAS_ROOT "/home/jmann")

#AOCL config
set(AOCL_ROOT "/home/jmann")

# Boost Configuration
set(BOOST_INC "/usr/include/boost")

#AOCC config
set(AOCC_ROOT "/home/jmann")

add_subdirectory(src)

# Include directories
#target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${BOOST_INC_ROOT} ${MKL_INC_ROOT} ${MPI_INCLUDE_PATH})
#target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${MPI_INCLUDE_PATH} ${AOCL_INC} ${BLIS_INC} ${BOOST_INC_ROOT} ${OBLAS_INC})
#target_include_directories(Timings PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${MPI_INCLUDE_PATH} ${AOCL_INC} ${BLIS_INC} ${BOOST_INC_ROOT} ${OBLAS_INC})

# Link libraries
#target_link_libraries(TDSESolveLinux PRIVATE OpenMP::OpenMP_CXX MPI::MPI_CXX ${MKL_LIBS} ${BOOST_LIBS})
#target_link_libraries(TDSESolveLinux PRIVATE MPI::MPI_CXX ${AOCC_LIBS} ${BOOST_LINK_FLAGS} ${AOCL_LIBS} ${OBLAS_LIBS})
#target_link_libraries(Timings PRIVATE MPI::MPI_CXX ${AOCC_LIBS} ${BOOST_LINK_FLAGS} ${AOCL_LIBS} ${OBLAS_LIBS})