#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "SimulationManager.h"
#include "MathTools.h"

// couple particular grid with SimulationManager for ease of use
class Sim 
    : public SimulationManager {
    private:
        double* x;
        int nPts;
    public:
        Sim(double xmin, double xmax, double dx, double dt, double maxT)
            : nPts((int)((xmax-xmin)/dx)), 
              SimulationManager(nPts, dx, dt, maxT, NULL){
                this->x = new double[nPts];
                for(int i = 0; i < nPts; i++)
                    this->x[i] = i*dx + xmin;
        }
        ~Sim(){delete[] x;}
        double* getXPtr(){return x;}
        std::vector<double> getX(){return std::vector<double>(x, x + nPts);}
    };

namespace py = pybind11;
using namespace pybind11::literals; //easy args

PYBIND11_MODULE(tdsepy, m) {
    py::class_<Sim>(m, "Sim")
        .def(py::init<double, double, double, double, double>(), 
            "Initialize a Simulation instance.",
            "xmin"_a, "xmax"_a, "dx"_a, "dt"_a, "maxT"_a)
        .def("getX", &Sim::getX)
        .def("getXPtr", &Sim::getXPtr)
        .def("getDX", &Sim::getDX);
}