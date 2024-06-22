file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(core SHARED ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" ${OBLAS_INC} ${BOOST_INC} PkgConfig::FFTW)
target_link_libraries(core PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX ${BOOST_LINK_FLAGS} PkgConfig::FFTW ${FFTW_DOUBLE_OMP_LIB} ${OBLAS_LIBS})