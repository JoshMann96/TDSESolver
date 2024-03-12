#files
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(TDSEpy SHARED ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(TDSEpy PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" ${MPI_INCLUDE_PATH} ${BOOST_INC_ROOT})
target_link_libraries(TDSEpy ${BOOST_LINK_FLAGS} core)

#generate Python bindings scripts
execute_process(COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/generatePythonBindingsFile.sh ${CMAKE_CURRENT_SOURCE_DIR}/Bindings.cpp ${CMAKE_CURRENT_SOURCE_DIR}/TDSEpyBindings.py)
#install Python bindings scripts
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/TDSEpyBindings.py ${CMAKE_BINARY_DIR}/TDSEpyBindings.py COPYONLY)