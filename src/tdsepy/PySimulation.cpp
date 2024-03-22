#include "PySimulation.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"
#include <memory>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>

void init_Simulation(py::module &m) {

// WEIGHT CALCULATIONS

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

// DENSITY PROCESSING

    py::class_<WfcToRho::Density>(m, "Density");

    py::class_<WfcToRho::DirectDensity, WfcToRho::Density>(m, "DirectDensity")
        .def(py::init<>(), R"V0G0N(
            Uses no preprocessing in calculating the final density.
            Density = sum over states (weight x psi*psi)

            Returns
            -------
            DirectDensity)V0G0N");

// KINETIC OPERATORS

    py::class_<KineticOperators::KineticOperator_PSM>(m, "KineticOperator_PSM");

    py::class_<KineticOperators::KineticOperator_FDM>(m, "KineticOperator_FDM");

    py::class_<KineticOperators::GenDisp_PSM_FreeElec, KineticOperators::KineticOperator_PSM>(m, "PSM_FreeElec")
        .def(py::init([](PySimulation* sim, double meff){
            return std::unique_ptr<KineticOperators::GenDisp_PSM_FreeElec>(new KineticOperators::GenDisp_PSM_FreeElec(sim->getNumPoints(), sim->getDX(), sim->getDT(), meff));
        }), R"V0G0N(
            Free electron dispersion relation with uniform effective mass using pseudospectral derivatives.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            meff : float
                Effective mass (1.0 = free electron).

            Returns
            -------
            PSM_FreeElec)V0G0N",
            "sim"_a, "meff"_a);

// SIMULATION

    py::class_<PySimulation>(m, "Simulation")
        .def(py::init<double, double, double, double, double, std::function<void(int)>>(), R"V0G0N(
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
            callback : function
                Callback -- reports percentage complete of time-stepping runs.

            Returns
            -------
            Simulation)V0G0N",
            "xmin"_a, "xmax"_a, "dx"_a, "dt"_a, "maxT"_a, "callback"_a)
        .def("__enter__", &PySimulation::enter)
        .def("__exit__", &PySimulation::exit)
        .def("getX", &PySimulation::getX)
        .def("getDX", &PySimulation::getDX)
        .def("addPot", &PySimulation::addPotential, py::keep_alive<1,2>(), R"V0G0N(
            Adds potential to the simulation.

            Parameters
            ----------
            pot : Potential
                Potential to be added.)V0G0N",
            "pot"_a)
        .def("addMeas", &PySimulation::addMeasurer, py::keep_alive<1,2>(), R"V0G0N(
            Adds measurer to the simulation.

            Parameters
            ----------
            meas : Measurer
                Measurer to be added.)V0G0N",
            "meas"_a)
        .def("setDens", &PySimulation::setDens, py::keep_alive<1,2>(), R"V0G0N(
            Sets density calculator for simulation.

            Parameters
            ----------
            dens : Density
                Density calculator to be used.)V0G0N",
            "dens"_a)
        .def("setWght", &PySimulation::setWght, py::keep_alive<1,2>(), R"V0G0N(
            Sets state weight calculator for simulation.

            Parameters
            ----------
            wght : Weight
                Weight calculator to be used.)V0G0N",
            "wght"_a)
        .def("setKin", py::overload_cast<KineticOperators::KineticOperator_PSM*>(&PySimulation::setKineticOperator_PSM), py::keep_alive<1,2>(), R"V0G0N(
            Sets kinetic operator for simulation.

            Parameters
            ----------
            nkin : KineticOperator
                Kinetic operator to be used.)V0G0N",
            "nkin"_a)
        .def("setKin", py::overload_cast<KineticOperators::KineticOperator_FDM*>(&PySimulation::setKineticOperator_FDM), py::keep_alive<1,2>(), R"V0G0N(
            Sets kinetic operator for simulation.

            Parameters
            ----------
            nkin : KineticOperator
                Kinetic operator to be used.)V0G0N",
            "nkin"_a)
        .def("addLeftAbsBdy", &PySimulation::addLeftAbsBdy, R"V0G0N(
            Adds absorptive boundary to left side of simulation.
            Decay is applied by multiplying states near boundary by a number of magnitude less than one.
            Over a chage of time Dt a stationary wavefunction will be
                psi = psi_0 sigma(x)^(rate * Dt)
            With sigma(x) a polynomial smooth function, 1 on the inner boundary and 0 on the outer boundary.

            Parameters
            ----------
            rate : float
                Decay rate.
            width : float
                Width of boundary)V0G0N",
            "rate"_a, "width"_a)
        .def("addRightAbsBdy", &PySimulation::addLeftAbsBdy, R"V0G0N(
            Adds absorptive boundary to right side of simulation.
            Decay is applied by multiplying states near boundary by a number of magnitude less than one.
            Over a chage of time Dt a stationary wavefunction will be
                psi = psi_0 sigma(x)^(rate * Dt)
            With sigma(x) a polynomial smooth function, 1 on the inner boundary and 0 on the outer boundary.

            Parameters
            ----------
            rate : float
                Decay rate.
            width : float
                Width of boundary.)V0G0N",
            "rate"_a, "width"_a)
        .def("finishInit", &PySimulation::finishInitialization, R"V0G0N(
            Finishes initialization of simulation.)V0G0N")
        .def("eigenSolve", &PySimulation::findEigenStates, R"V0G0N(
            Finds eigenstates of current system, without self-consistent potentials.

            Parameters
            ----------
            minE : float
                Eigenvalue lower bound.
            maxE : float
                Eigenvalue upper bound)V0G0N",
            "minE"_a, "maxE"_a)
        .def("runOS_U2TU", &PySimulation::runOS_U2TU, R"V0G0N(
            Runs simulation using operator splitting method. Potential is not updated between kinetic operator propagation steps.)V0G0N")
        .def("runOS_UW2TUW", &PySimulation::runOS_UW2TUW, R"V0G0N(
            Runs simulation using operator splitting method. Potential is updated between kinetic operator propagation steps.)V0G0N");
}