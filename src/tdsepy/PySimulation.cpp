#include "PySimulation.h"

void init_Simulation(py::module &m) {
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
        .def("getDX", &PySimulation::getDX)
        .def("addPot", &PySimulation::addPotential, R"V0G0N(
            Adds potential to the simulation.

            Parameters
            ----------
            pot : Potential
                Potential to be added.)V0G0N",
            "pot"_a);
}