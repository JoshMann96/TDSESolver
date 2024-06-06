file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_executable(TDSESolveLinux ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(TDSESolveLinux PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${AOCL_INC} ${BOOST_INC_ROOT})
target_link_libraries(TDSESolveLinux PRIVATE MPI::MPI_CXX ${BOOST_LINK_FLAGS} ${AOCL_LIBS} core)