#pragma once
#include "PyCommon.h"
#include "PySimulation.h"
#include "Potentials.h"

// simplify instantiation of potentials
class PyFilePotential
    : public Potentials::FilePotential{
        public:
        PyFilePotential(PySimulation * sim, double offset, std::string fil, double refPoint)
            : Potentials::FilePotential(sim->getNumPoints(), sim->getXPtr(), offset, fil.c_str(), sim->findXIdx(refPoint)){}
    };

void init_Potentials(py::module &m);