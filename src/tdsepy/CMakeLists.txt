# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=znver2 -mtune=native -DNDEBUG -fopenmp -fPIC")

#files
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(TDSEpy SHARED ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(TDSEpy PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(TDSEpy core)

#generate Python bindings scripts
execute_process(COMMAND bash generatePythonBindingsFile.sh Bindings.cpp TDSEpyBindings.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#install Python bindings scripts
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/TDSEpyBindings.py ${CMAKE_BINARY_DIR}/TDSEpyBindings.py COPYONLY)