file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(core SHARED ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" ${FFTW_INC} ${BOOST_INC_ROOT} ${OBLAS_INC})
target_link_libraries(core PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX ${BOOST_LINK_FLAGS} ${FFTW_LIBS} ${OBLAS_LIBS})