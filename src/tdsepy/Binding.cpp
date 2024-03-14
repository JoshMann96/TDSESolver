#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "Potentials.h"
#include "Measurers.h"
#include "SimulationManager.h"
#include "MathTools.h"

// couple particular grid with SimulationManager for ease of use
class Simulation 
    : public SimulationManager {
    private:
        double* x;
        int nPts;
    public:
        Simulation(double xmin, double xmax, double dx, double dt, double maxT)
            : nPts((int)((xmax-xmin)/dx)), 
              SimulationManager(nPts, dx, dt, maxT, NULL){
                this->x = new double[nPts];
                for(int i = 0; i < nPts; i++)
                    this->x[i] = i*dx + xmin;
        }
        ~Simulation(){delete[] x;}
        void* getXPtr(){return x;}
        std::vector<double> getX(){return std::vector<double>(x, x + nPts);}
        int checkPtr(void* ix){return reinterpret_cast<double*>(ix) == x;}
    };

namespace py = pybind11;
using namespace pybind11::literals; //easy args

PYBIND11_MODULE(tdsepy, m) {
    py::class_<Simulation>(m, "Simulation")
        .def(py::init<double, double, double, double, double>(), R"V0G0N(
            Creates a TDSE Simulation instance.

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
            Simulation
                Instance of TDSE Simulator.)V0G0N",
            "xmin"_a, "xmax"_a, "dx"_a, "dt"_a, "maxT"_a)
        .def("getX", &Simulation::getX)
        .def("getXPtr", &Simulation::getXPtr)
        .def("getDX", &Simulation::getDX)
        .def("checkPtr", &Simulation::checkPtr);
}