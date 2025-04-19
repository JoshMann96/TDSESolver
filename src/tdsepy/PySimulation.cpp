#include "PySimulation.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"
#include <memory>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pytypes.h>

void init_Simulation(py::module &m) {
    py::class_<PySimulation>(m, "Simulation")
        .def(py::init<double, double, double, double, std::function<void(int)>>(), R"V0G0N(
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
            callback : function
                Callback -- reports percentage complete of time-stepping runs.

            Returns
            -------
            Simulation)V0G0N",
            "xmin"_a, "xmax"_a, "dx"_a, "dt"_a, "callback"_a)
        .def("getXVec", &PySimulation::getXVec)
        .def("getDX", &PySimulation::getDX)
        .def("findXIdx", &PySimulation::findXIdx, R"V0G0N(
            Finds index of position in grid.

            Parameters
            ----------
            xp : float
                Position.

            Returns
            -------
            uint)V0G0N",
            "xp"_a)
        .def("addPot", &PySimulation::addPotential, py::keep_alive<1,2>(), R"V0G0N(
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
            "meas"_a) // no need for keep_alive as Measurer has py::nodelete, the MeasurementManager will take care of the memory
        .def("setDens", &PySimulation::setDensity, py::keep_alive<1,2>(), R"V0G0N(
            Sets density calculator for simulation.

            Parameters
            ----------
            dens : Density
                Density calculator to be used.)V0G0N",
            "dens"_a)
        .def("setWght", &PySimulation::setWeight, py::keep_alive<1,2>(), R"V0G0N(
            Sets state weight calculator for simulation.

            Parameters
            ----------
            wght : Weight
                Weight calculator to be used.)V0G0N",
            "wght"_a)
        .def("setKin", py::overload_cast<KineticOperators::KineticOperator*>(&PySimulation::setKineticOperator), py::keep_alive<1,2>(), R"V0G0N(
            Sets kinetic operator for simulation.

            Parameters
            ----------
            nkin : KineticOperator
                Kinetic operator to be used.)V0G0N",
            "nkin"_a)
        .def("addLeftAbsBdy", &PySimulation::addLeftAbsBdy, R"V0G0N(
            Adds absorptive boundary to left side of simulation.
            Decay is applied by multiplying states near boundary by a number of magnitude less than one.
            Over a change of time Dt a stationary wavefunction will be
                psi = psi_0 sigma(x)^(rate * Dt)
            With sigma(x) a polynomial smooth function, 1 on the inner boundary and 0 on the outer boundary.

            Parameters
            ----------
            rate : float
                Decay rate.
            width : float
                Width of boundary)V0G0N",
            "rate"_a, "width"_a)
        .def("addRightAbsBdy", &PySimulation::addRightAbsBdy, R"V0G0N(
            Adds absorptive boundary to right side of simulation.
            Decay is applied by multiplying states near boundary by a number of magnitude less than one.
            Over a change of time Dt a stationary wavefunction will be
                psi = psi_0 sigma(x)^(rate * Dt)
            With sigma(x) a polynomial smooth function, 1 on the inner boundary and 0 on the outer boundary.

            Parameters
            ----------
            rate : float
                Decay rate.
            width : float
                Width of boundary.)V0G0N",
            "rate"_a, "width"_a)
        .def("findEigenStates", &PySimulation::findEigenStates, R"V0G0N(
            Finds eigenstates of current system, without self-consistent potentials.

            Parameters
            ----------
            minE : float
                Eigenvalue lower bound.
            maxE : float
                Eigenvalue upper bound)V0G0N",
            "minE"_a, "maxE"_a)
        .def("findInhomogeneousEigenStates", &PySimulation::findInhomogeneousEigenStates, R"V0G0N(
            Finds eigenstates of current system, with inhomogeneous boundary conditions.
            This is only intended to work for finite difference schemes with supported boundary conditions.
            The resulting states should be orthogonal eigenstates of the open system.

            Parameters
            ----------
            nElec : uint
                Number of electrons in the system.
            energies : float array
                Eigenstate energies. The boundary conditions must be consistent with these energies.)V0G0N",
            "nElec"_a, "energies"_a)
        .def("run", &PySimulation::run, R"V0G0N(
            Runs simulation using the appropriate method for the kinetic operator and potential.)V0G0N",
            "nSteps"_a)
        .def("runEPS_U2TU", &PySimulation::runEPS_U2TU, R"V0G0N(
            Runs simulation using operator splitting method. Potential is not updated between kinetic operator propagation steps.)V0G0N",
            "nSteps"_a)
        .def("runEPS_UW2TUW", &PySimulation::runEPS_UW2TUW, R"V0G0N(
            Runs simulation using operator splitting method. Potential is updated between kinetic operator propagation steps.)V0G0N",
            "nSteps"_a)
        .def("runCN_L", &PySimulation::runCN_L, R"V0G0N(
            Runs simulation using Crank-Nicolson method. Potential is not updated between kinetic operator propagation steps.)V0G0N",
            "nSteps"_a)
        .def("runCN_NL", &PySimulation::runCN_NL, R"V0G0N(
            Runs simulation using Crank-Nicolson method. Potential is updated between kinetic operator propagation steps.)V0G0N",
            "nSteps"_a)
        .def("getElectricalCentroidSurface", &PySimulation::findElectricalSurfaceCentroidRule, R"V0G0N(
            Finds the index of the electrical surface using the centroid rule.

            Parameters
            ----------
            minPos : float
                Minimum position for region including centroid of response.
            maxPos : float
                Maximum position for region including centroid of response.

            Returns
            -------
            uint)V0G0N",
            "minPos"_a, "maxPos"_a);
}