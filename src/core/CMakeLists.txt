set(EXTERN_DIR "${CMAKE_SOURCE_DIR}/src/extern")

set(TEST_SRC "${CMAKE_CURRENT_SOURCE_DIR}/test.cpp")
get_filename_component(TEST_SRC_ABS ${TEST_SRC} ABSOLUTE)

file(GLOB CUDA_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Cu_*")
message(STATUS "CUDA_SRC: ${CUDA_SRC}")

#cuda stuff
if(USE_CUDA)
    file(GLOB CUSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cu")
    file(GLOB CUHDR "${CMAKE_CURRENT_SOURCE_DIR}/*.cuh")
    #set(CMAKE_CUDA_ARCHITECTURES 50 52 60 61 70 72 75 80 86 87 89 90 100 101 120)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets")
    add_library(myCuda SHARED ${CUSRC} ${CUHDR})
    target_compile_options(myCuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
    target_compile_features(myCuda PRIVATE cxx_std_17 cuda_std_17)
    set_target_properties(myCuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
else()
    set(CUSRC "")
    set(CUHDR "")
endif()

#cpu stuff
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
list(REMOVE_ITEM SRC ${TEST_SRC_ABS})
if(NOT USE_CUDA)
    list(REMOVE_ITEM SRC ${CUDA_SRC})
endif()
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
if(NOT USE_CUDA)
    list(REMOVE_ITEM HDR ${CUSRC})
endif()
file(GLOB EXTRN "${EXTERN_DIR}/*.h")

add_library(core SHARED ${SRC} ${IMP} ${HRP} ${HDR} ${EXTRN} ${CUHDR})
target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" Boost::boost PkgConfig::FFTW "${EXTERN_DIR}")
if(USE_CUDA)
    target_link_libraries(core PUBLIC myCuda MPI::MPI_CXX OpenMP::OpenMP_CXX ${Boost_LIBRARIES} PkgConfig::FFTW ${FFTW_DOUBLE_OMP_LIB} ${BLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_LIBRARIES} ${CUDA_cublas_LIBRARY})
else()
    target_link_libraries(core PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX ${Boost_LIBRARIES} PkgConfig::FFTW ${FFTW_DOUBLE_OMP_LIB} ${BLAS_LIBRARIES})
endif()


add_executable(test ${TEST_SRC})
target_include_directories(test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" Boost::boost PkgConfig::FFTW "${EXTERN_DIR}")
target_link_libraries(test PRIVATE core)