#include "PyPotentials.h"

void init_Potentials(py::module &m) {
    py::class_<Potentials::Potential>(m, "Potential");

    py::class_<PyFilePotential, Potentials::Potential>(m, "FilePotential")
        .def(py::init<PySimulation*, double, const std::string, double>(), R"V0G0N(
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