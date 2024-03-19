#include "PyCommon.h"
#include "PySimulation.h"
#include "PyPotentials.h"
#include "PyMeasurers.h"

PYBIND11_MODULE(tdsepy, m) {
    py::module_ mp = m.def_submodule("Potentials", "Potentials for TDSE solver.");
    py::module_ mm = m.def_submodule("Measurers", "Measurers for TDSE solver.");

    init_Simulation(m);
    init_Potentials(mp);
    init_Measurers(mm);
    //init_Measurers(mm);
}