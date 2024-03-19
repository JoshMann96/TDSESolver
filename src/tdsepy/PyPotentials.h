#pragma once
#include "PyCommon.h"
#include "PySimulation.h"
#include "Potentials.h"

// simplify instantiation of potentials
class PyFilePotential
    : public Potentials::FilePotential{
        private:
        int nPts;
        public:
        PyFilePotential(PySimulation * sim, double offset, std::string fil, double refPoint)
            : nPts(sim->getNumPoints()), Potentials::FilePotential(sim->getNumPoints(), sim->getXPtr(), offset, fil.c_str(), sim->findXIdx(refPoint)){}
        std::vector<double> getValue();
    };

void init_Potentials(py::module &m);