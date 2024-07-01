file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_executable(TDSESolver ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(TDSESolver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" ${BOOST_INC_ROOT})
target_link_libraries(TDSESolver PRIVATE core)