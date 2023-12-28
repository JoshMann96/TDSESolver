cmake_minimum_required(VERSION 3.13)
project(TDSESolveLinux)

# Compiler options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -DNDEBUG")

# Boost Configuration
set(BOOST_INC_ROOT "/usr/include/boost")
include_directories(${BOOST_INC_ROOT})

# Intel MKL Configuration
set(MKLROOT "/opt/intel/oneapi/mkl/latest")
include_directories(${MKLROOT}/include)
link_directories(${MKLROOT}/lib/intel64)
set(MKL_LIBS "-Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl")

# OpenMP Configuration
find_package(OpenMP REQUIRED)


# Add source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Add the executable
add_executable(TDSESolveLinux ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${BOOST_INC_ROOT})

# Link directories
link_directories(${Boost_LIBRARY_DIRS})

# Link libraries
target_link_libraries(TDSESolveLinux PRIVATE ${MKL_LIBS} OpenMP::OpenMP_CXX)