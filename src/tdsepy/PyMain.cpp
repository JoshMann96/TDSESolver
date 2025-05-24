#include "PyCommon.h"
#include "PySimulation.h"
#include "PyKinetics.h"
#include "PyDensities.h"
#include "PyPotentials.h"
#include "PyMeasurers.h"

PYBIND11_MODULE(tdsepy, m) {
    py::module_ ms = m.def_submodule("Simulation", "Simulation for TDSE solver.");
    py::module_ mk = m.def_submodule("Kinetics", "Kinetic operators for TDSE solver.");
    py::module_ md = m.def_submodule("Densities", "Density calculators for TDSE solver.");
    py::module_ mp = m.def_submodule("Potentials", "Potentials for TDSE solver.");
    py::module_ mm = m.def_submodule("Measurers", "Measurers for TDSE solver.");

    init_Simulation(ms);
    init_Kinetics(mk);
    init_Densities(md);
    init_Potentials(mp);
    init_Measurers(mm);
}