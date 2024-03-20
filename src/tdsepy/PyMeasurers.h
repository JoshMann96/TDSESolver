#pragma once
#include "PyCommon.h"
#include "Measurers.h"
#include "PySimulation.h"

void init_Measurers(py::module &m);

using namespace Measurers;

class PyConstant
    : public DoubleConst {
        public:
        PyConstant(double c, std::string filName, std::string fol)
            : DoubleConst(c, filName.c_str(), fol.c_str()){}
    };

class PyBasic
    : public BasicMeasurers{
        public:
        PyBasic(PySimulation* sim, std::string fol)
            : BasicMeasurers("_py_", sim->getNumPoints(), sim->getDX(), sim->getDT(), sim->getXPtr(), fol.c_str()){}
    };