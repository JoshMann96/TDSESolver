#!/bin/bash
set -e

# This script is used to build the project and create a Python virtual environment (venv).
# To install to a custom venv, cd to build/lib/tdsepy and run `pip install -e .`

# This project uses CMake's find_package functionality for the following required packages
# OpenMP
# Boost
# FFTW3 & FFTW3_OMP (Double) [-DFFTW_ROOT=path]
# OpenBLAS (highly recommend OMP version, only the library is needed) [-DBLAS_HINTS=path/to/lib]

# You can also specify the following optional parameters:
# Cuda                      [-DUSE_CUDA=ON|OFF]
# Disable Python bindings   [--disablePython]
# Debugging comp. flags     [-DDEBUG=ON|OFF]

original_params=("$@")
echo "${original_params[@]}"

PYTHON_BINDINGS=TRUE

while test $# -gt 0
do
    case "$1" in
        --disablePython) 
            PYTHON_BINDINGS=FALSE
            echo "Disabling Python bindings"
            ;;
    esac
    shift
done

if [ "$PYTHON_BINDINGS" = TRUE ]
    then
    python3 -m venv .venv
    source .venv/bin/activate

    python_version=$(python3 --version)
    version_numbers=(${python_version//./ })
    export PYTHON_SUBVERSION_NUMBER=${version_numbers[2]}

    pip3 install --upgrade pip
    pip3 install -r requirements.txt
    fi

if [ -d build ];
then echo "build folder already exists."
else
mkdir build
fi

filtered_params=()
for param in "${original_params[@]}"; do
    if [[ $param == -D* ]]; then
        filtered_params+=("$param")
    fi
done

cmake -S . -B build "${filtered_params[@]}" -DPYTHON_BINDINGS=$PYTHON_BINDINGS
cd build
make -j 8

if [ "$PYTHON_BINDINGS" = TRUE ]
    then
    cd lib/tdsepy
    pip3 install .

    echo ""
    echo "BUILD COMPLETE"
    echo "To install to a custom venv"
    echo "    1. source the desired venv"
    echo "    2. cd to build/lib/tdsepy"
    echo "    3. run 'pip install .'"

    deactivate

    cd ../..

    fi
cd ..