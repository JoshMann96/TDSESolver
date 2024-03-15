#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "Potentials.h"
#include "Measurers.h"
#include "SimulationManager.h"
#include "MathTools.h"

// couple particular grid with SimulationManager for ease of use
class PySimulation 
    : public SimulationManager {
    private:
        double* x;
        int nPts;
    public:
        PySimulation(double xmin, double xmax, double dx, double dt, double maxT)
            : nPts((int)((xmax-xmin)/dx)), 
              SimulationManager(nPts, dx, dt, maxT, NULL){
                this->x = new double[nPts];
                for(int i = 0; i < nPts; i++)
                    this->x[i] = i*dx + xmin;
        }
        ~PySimulation(){delete[] x;}
        double* getXPtr(){return x;}
        std::vector<double> getX(){return std::vector<double>(x, x + nPts);}
        int findXIdx(double xp){return vtls::findValue(nPts, x, xp);}
    };

// simplify instantiation of potentials
class PyFilePotential
    : public Potentials::FilePotential{
        public:
        PyFilePotential(PySimulation sim, double offset, std::string fil, double refPoint)
            : Potentials::FilePotential(sim.getNumPoints(), sim.getXPtr(), offset, fil.c_str(), sim.findXIdx(refPoint)){}
    };

namespace py = pybind11;
using namespace py::literals; //easy args

PYBIND11_MODULE(tdsepy, m) {
    py::class_<PySimulation>(m, "Simulation")
        .def(py::init<double, double, double, double, double>(), R"V0G0N(
            Manages TDSE simulations.

            Parameters
            ----------
            xmin : float
                Left boundary position.
            xmax : float
                Right boundary position.
            dx : float
                Spatial step size.
            dt : float
                Temporal step size.
            maxT : float
                Maximum time in simulation (end time).

            Returns
            -------
            Simulation)V0G0N",
            "xmin"_a, "xmax"_a, "dx"_a, "dt"_a, "maxT"_a)
        .def("getX", &PySimulation::getX)
        .def("getDX", &PySimulation::getDX);
    
    py::class_<Potentials::Potential>(m, "Potential");

    py::class_<PyFilePotential, Potentials::Potential>(m, "FilePotential")
        .def(py::init<PySimulation, double, const std::string, double>(), R"V0G0N(
            Includes a static potential as defined in a binary file.
            See documentation for appropriate file format.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            offset : float
                Translational offset with respect to input data.
            fil : str
                File path with data.
            refPoint : float
                Potential reference point.

            Returns
            -------
            FilePotential)V0G0N",
            "sim"_a, "offset"_a, "fil"_a, "refPoint"_a);
}