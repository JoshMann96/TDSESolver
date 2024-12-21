set(EXTERN_DIR "${CMAKE_SOURCE_DIR}/src/extern")
set(TEST_SRC "${CMAKE_CURRENT_SOURCE_DIR}/test.cpp")
get_filename_component(TEST_SRC_ABS ${TEST_SRC} ABSOLUTE)

file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
list(REMOVE_ITEM SRC ${TEST_SRC_ABS})
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB EXTRN "${EXTERN_DIR}/*.h")

add_library(core SHARED ${SRC} ${IMP} ${HRP} ${HDR} ${EXTRN})
target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" Boost::boost PkgConfig::FFTW "${EXTERN_DIR}")
target_link_libraries(core PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX ${Boost_LIBRARIES} PkgConfig::FFTW ${FFTW_DOUBLE_OMP_LIB} ${BLAS_LIBRARIES})

add_executable(test ${TEST_SRC})
target_include_directories(test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" Boost::boost PkgConfig::FFTW "${EXTERN_DIR}")
target_link_libraries(test PRIVATE core)