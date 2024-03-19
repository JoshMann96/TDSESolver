#pragma once
#include "PyCommon.h"
#include "Measurers.h"

void init_Measurers(py::module &m);

using namespace Measurers;

class PyConstant
    : DoubleConst {
        public:
        PyConstant(double c, std::string filName, std::string fol)
            : DoubleConst(c, filName.c_str(), fol.c_str()){}
    };

