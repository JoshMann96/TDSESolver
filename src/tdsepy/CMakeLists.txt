# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=znver2 -mtune=native -DNDEBUG -fopenmp -fPIC")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -glldb -fopenmp -fPIC")
#set(CMAKE_BUILD_TYPE Debug)

#PYBIND11
set(pybind11_DIR ${CMAKE_SOURCE_DIR}/.venv/lib/python3.10/site-packages/pybind11/share/cmake/pybind11/)
find_package(pybind11 REQUIRED)

#files
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

pybind11_add_module(tdsepy SHARED ${SRC} ${IMP} ${HRP} ${HDR} )
target_link_libraries(tdsepy PRIVATE core pybind11::pybind11)

# Use stubgen to create .pyi files to sit alongside the just-built python module
set(Stubgen_Executable "${CMAKE_SOURCE_DIR}/.venv/bin/stubgen")
if(EXISTS "${Stubgen_Executable}")
    add_custom_command(TARGET tdsepy POST_BUILD
        COMMAND ${Stubgen_Executable} --include-docstrings -m tdsepy -o .
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMENT "Use stubgen to create .pyi for statement completion")
else()
    message( SEND_ERROR "${Stubgen_Executable} does not exist -- .pyi file not generated. Ensure that mypy is installed.")
endif()

# Copy tdsepy_helpers to build dir
add_custom_target(tdsepy_helpers ALL DEPENDS tdsepy)

add_custom_command(TARGET tdsepy_helpers POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/helpers ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tdsepy/helpers)

# Copy setup.py to build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/setup.py COPYONLY)