cmake_minimum_required(VERSION 3.13)
project(TDSESolveLinux)

# Compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -DNDEBUG")

# Boost Configuration
set(BOOST_INC "/usr/include/boost")
include_directories(${BOOST_INC})
set(BOOST_LIBS "-lboost_system -lboost_filesystem")

# Intel MKL Configuration
set(MKL_INC "/usr/include/mkl")
set(MKL_LIB "/usr/lib/x86_64-linux-gnu/mkl/lib")
include_directories(${MKL_INC})
link_directories(${MKL_LIB})
set(MKL_LIBS "-Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl")

# OpenMP Configuration
find_package(OpenMP REQUIRED)


# Add source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Add the executable
add_executable(TDSESolveLinux ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${BOOST_INC_ROOT} ${MKL_INC_ROOT})

# Link libraries
target_link_libraries(TDSESolveLinux PRIVATE OpenMP::OpenMP_CXX ${MKL_LIBS} ${BOOST_LIBS})