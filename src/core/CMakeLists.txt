#Boost
include_directories(${BOOST_INC})
set(BOOST_LINK_FLAGS "-lboost_system -lboost_filesystem")

#MPI Configuration
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#AOCC
set(AOCC_INC "${AOCC_ROOT}/include")
set(AOCC_LIB "${AOCC_ROOT}/lib")
set(AOCC_LIBS "-L${AOCC_LIB} -lomp -lflang -lstdc++")
include_directories(${AOCC_INC})

#AOCL config
set(AOCL_INC "${AOCL_ROOT}/include")
set(AOCL_LIB "${AOCL_ROOT}/lib")
set(BLIS_INC "${AOCL_INC}/blis")
include_directories(${AOCL_INC})
set(AOCL_LIBS "-L${AOCL_LIB} -lfftw3_omp -lfftw3 -lm -lpthread")

#OpenBLAS config
set(OBLAS_INC "${OBLAS_ROOT}/include")
set(OBLAS_LIB "${OBLAS_ROOT}/lib")
set(OBLAS_LIBS "-L${OBLAS_LIB} -lopenblas -lflang -lpthread")

#files
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(core_lib STATIC ${SRC} ${IMP} ${HDR})
target_include_directories(core_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" ${AOCL_INC} ${BLIS_INC} ${BOOST_INC_ROOT} ${OBLAS_INC})
target_link_libraries(core_lib PRIVATE MPI::MPI_CXX ${AOCC_LIBS} ${BOOST_LINK_FLAGS} ${AOCL_LIBS} ${OBLAS_LIBS})