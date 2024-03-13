# Compiler options
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=znver2 -mtune=native -DNDEBUG -fopenmp -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=znver2 -mtune=native -fopenmp -fPIC")

#files
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(core STATIC ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" ${AOCL_INC} ${BLIS_INC} ${BOOST_INC_ROOT} ${OBLAS_INC})
target_link_libraries(core PRIVATE MPI::MPI_CXX ${AOCC_LIBS} ${BOOST_LINK_FLAGS} ${AOCL_LIBS} ${OBLAS_LIBS})