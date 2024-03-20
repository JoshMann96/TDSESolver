#include "PyMeasurers.h"

using namespace Measurers;

void init_Measurers(py::module &m) {
    py::class_<Measurer>(m, "Measurer");

    py::class_<PyConstant, Measurer>(m, "Constant")
        .def(py::init<double, std::string, std::string>(), R"V0G0N(
            Records a constant.

            Parameters
            ----------
            c : float
                Value to be recorded.
            fileName : str
                Name of output file (no extension).
            fol : str
                Directory to contain file.

            Returns
            -------
            Constant)V0G0N",
            "c"_a, "fileName"_a, "fol"_a);
    
    py::class_<PyBasic, Measurer>(m, "Basic")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records basic simulation parameters:
             - nPts     - Number of points in grid
             - NSteps   - Number of time steps
             - DX       - Spatial step size
             - DT       - Temporal step size

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            Basic)V0G0N",
            "sim"_a, "fol"_a);
}