# Compiler options
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=znver2 -mtune=native -DNDEBUG -fopenmp -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=znver2 -mtune=native -fopenmp -fPIC")

set(BUILD_PY_COMMAND bash generatePythonBindingsFile.sh Bindings.cpp TDSEpyBindings.py)

#files
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_library(TDSEpy SHARED ${SRC} ${IMP} ${HRP} ${HDR})
target_include_directories(TDSEpy PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(TDSEpy core)

# Define the bash script path

add_custom_target(buildPy
    COMMAND ${BUILD_PY_COMMAND}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building python scripts"
)

add_custom_target(installPy
    DEPENDS buildPy
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/TDSEpyBindings.py ${CMAKE_BINARY_DIR}
    COMMENT "Installing python scripts"
)

add_dependencies(installPy buildPy) #must make py code before installing

add_dependencies(TDSEpy installPy) #must install

#generate Python bindings scripts
#execute_process(COMMAND bash generatePythonBindingsFile.sh Bindings.cpp TDSEpyBindings.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#install Python bindings scripts
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/TDSEpyBindings.py ${CMAKE_BINARY_DIR}/TDSEpyBindings.py COPYONLY)
#message complete
#message("-- Built Python bindings wrapper") 