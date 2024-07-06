#!/bin/bash

# This script is used to build the project and create a Python virtual environment (venv).
# To install to a custom venv, cd to build/lib/tdsepy and run `pip install -e .`

# This project uses CMake's find_package functionality for the following (required) packages
# OpenMP
# MPI
# Boost
# FFTW3 & FFTW3_OMP (Double) [-DFFTW_ROOT=path]
# OpenBLAS (highly recommend OMP version, only the library is needed) [-DBLAS_HINTS=path/to/lib]

if [ -d .venv ]; 
then echo ".venv folder already exists. Deleting..." & rm -r .venv
fi

python3 -m venv .venv
source .venv/bin/activate

python_version=$(python3 --version)
version_numbers=(${python_version//./ })
export PYTHON_SUBVERSION_NUMBER=${version_numbers[2]}

pip3 install --upgrade pip
pip3 install -r requirements.txt

if [ -d build ];
then echo "build folder already exists."
else
mkdir build
fi

cmake -S . -B build "$@"
cd build
make -j 8

cd lib/tdsepy
pip3 install -e .

echo ""
echo "BUILD COMPLETE"
echo "To install to a custom venv"
echo "    1. source the desired venv"
echo "    2. cd to build/lib/tdsepy"
echo "    3. run 'pip install -e .'"

deactivate