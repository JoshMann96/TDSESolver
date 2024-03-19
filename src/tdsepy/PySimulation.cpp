#include "PySimulation.h"
#include "WfcRhoTools.h"

void init_Simulation(py::module &m) {
    py::class_<WfcToRho::Weight>(m, "Weight");

    py::class_<WfcToRho::FermiGasDistro, WfcToRho::Weight>(m, "FermiGasDistro")
        .def(py::init<double>(), R"V0G0N(
            Uses 3-D Fermi gas distribution to convert 1-D wavefunctions to an effective 3-D density.
            Mapping uses wavefunction initial eigenstates.

            Parameters
            ----------
            ef : float
                Fermi energy.

            Returns
            -------
            FermiGasDistro)V0G0N",
            "ef"_a);

    py::class_<WfcToRho::Density>(m, "Density");

    py::class_<WfcToRho::DirectDensity, WfcToRho::Density>(m, "DirectDensity")
        .def(py::init<>(), R"V0G0N(
            Uses no preprocessing in calculating the final density.
            Density = sum over states (weight x psi*psi)

            Returns
            -------
            DirectDensity)V0G0N");

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
            "pot"_a)
        .def("addMeas", &PySimulation::addMeasurer, R"V0G0N(
            Adds measurer to the simulation.

            Parameters
            ----------
            meas : Measurer
                Measurer to be added.)V0G0N",
            "meas"_a)
        .def("setDens", &PySimulation::setDens, R"V0G0N(
            Sets density calculator for simulation.

            Parameters
            ----------
            dens : Density
                Density calculator to be used.)V0G0N",
            "dens"_a)
        .def("setWght", &PySimulation::setWght, R"V0G0N(
            Sets state weight calculator for simulation.

            Parameters
            ----------
            wght : Weight
                Weight calculator to be used.)V0G0N",
            "wght"_a);
}