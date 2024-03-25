# Compiler options
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=znver2 -mtune=native -DNDEBUG -fopenmp -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -glldb -fopenmp -fPIC")
set(CMAKE_BUILD_TYPE Debug)

#PYBIND11
find_package(pybind11)

#files
file(GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HRP "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
file(GLOB IMP "${CMAKE_CURRENT_SOURCE_DIR}/*.tpp")
file(GLOB HDR "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

pybind11_add_module(tdsepy MODULE ${SRC} ${IMP} ${HRP} ${HDR})
target_link_libraries(tdsepy PRIVATE core pybind11::pybind11)

# Use stubgen to create .pyi files to sit alongside the just-built python module
set(Stubgen_Executable "${CMAKE_SOURCE_DIR}/.venv/bin/stubgen")
if(EXISTS "${Stubgen_Executable}")
    add_custom_command(TARGET tdsepy POST_BUILD
        COMMAND ${Stubgen_Executable} --include-docstrings -p tdsepy -o .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Use stubgen to create .pyi for statement completion")
else()
    message( SEND_ERROR "${Stubgen_Executable} does not exist -- .pyi file not generated. Ensure that mypy is installed.")
endif()