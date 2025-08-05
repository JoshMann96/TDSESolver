#include "PySimulation.h"
#include "KineticOperator.h"
#include "Densities.h"
#include <memory>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pytypes.h>

void init_Simulation(py::module &m) {
    py::class_<SimulationManager>(m, "SimulationManager", R"V0G0N(
        DO NOT USE DIRECTLY. C++ class for managing TDSE simulations.
    )V0G0N");

    py::enum_<Densities::NormalizationScheme>(m, "NormalizationScheme")
        .value("UNNORMALIZED", Densities::UNNORMALIZED)
        .value("NORMALIZED", Densities::NORMALIZED);

    py::class_<PySimulation, SimulationManager>(m, "Simulation")
        .def(py::init<double, double, double, double, const std::optional<std::function<void(double)>>, std::optional<size_t>, std::optional<std::string>>(), R"V0G0N(
            Manages TDSE simulations.
            This constructor initializes the simulation with a specified spatial range, step size, and time step.

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
                Callback function. Takes in a float between 0 and 1.0 for the progress of the present calculation. Default is None (no callback).
            numCallbackCalls : int
                Number of times to call the callback function. Default is 101.
                The callback function will be called with doubles ranging from 0 to 1.0, inclusive.
            fftwWisdomPrefix : str
                Prefix for the FFTW wisdom file. Default is an empty string. Wisdom filenames are already unique for each distinct HOST_NAME.

            Returns
            -------
            Simulation)V0G0N",
            "xmin"_a, "xmax"_a, "dx"_a, "dt"_a, "callback"_a = py::none(), "numCallbackCalls"_a = 101, "fftwWisdomPrefix"_a = "")
        .def(py::init<size_t, double, double, double, const std::optional<std::function<void(double)>>, std::optional<size_t>, std::optional<std::string>>(), R"V0G0N(
            Manages TDSE simulations.
            This constructor initializes the simulation with a specified number of grid points, left boundary position, step size, and time step.

            Parameters
            ----------
            nPts : int
                Number of grid points in the simulation.
            xmin : float
                Left boundary position.
            dx : float
                Spatial step size.
            dt : float
                Temporal step size.
            callback : function
                Callback function. Takes in a float between 0 and 1.0 for the progress of the present calculation. Default is None (no callback).
            numCallbackCalls : int
                Number of times to call the callback function. Default is 101.
                The callback function will be called with doubles ranging from 0 to 1.0, inclusive.
            fftwWisdomPrefix : str
                Prefix for the FFTW wisdom file. Default is an empty string. Wisdom filenames are already unique for each distinct HOST_NAME.

            Returns
            -------
            Simulation)V0G0N",
            "nPts"_a, "xmin"_a, "dx"_a, "dt"_a, "callback"_a = py::none(), "numCallbackCalls"_a = 101, "fftwWisdomPrefix"_a = "")
        .def(py::init<double, double, size_t, double, const std::optional<std::function<void(double)>>, std::optional<size_t>, std::optional<std::string>>(), R"V0G0N(
            Manages TDSE simulations.
            This constructor initializes the simulation with a specified spatial range, number of grid points, and time step.

            Parameters
            ----------
            xmin : float
                Left boundary position.
            xmax : float
                Right boundary position.
            nPts : int
                Number of grid points in the simulation.
            dt : float
                Temporal step size.
            callback : function
                Callback function. Takes in a float between 0 and 1.0 for the progress of the present calculation. Default is None (no callback).
            numCallbackCalls : int
                Number of times to call the callback function. Default is 101.
                The callback function will be called with doubles ranging from 0 to 1.0, inclusive.
            fftwWisdomPrefix : str
                Prefix for the FFTW wisdom file. Default is an empty string. Wisdom filenames are already unique for each distinct HOST_NAME.

            Returns
            -------
            Simulation)V0G0N",
            "xmin"_a, "xmax"_a, "nPts"_a, "dt"_a, "callback"_a = py::none(), "numCallbackCalls"_a = 101, "fftwWisdomPrefix"_a = "")
        .def("getX", &PySimulation::getX)
        .def("getDX", &PySimulation::getDX)
        .def("getDT", &PySimulation::getDT)
        .def("getNumPoints", &PySimulation::getNumPoints, R"V0G0N(
            Returns the number of grid points in the simulation.

            Returns
            -------
            uint : Number of grid points.)V0G0N")
        .def("getNumStates", &PySimulation::getNElec, R"V0G0N(
            Returns the number of states in the simulation.

            Returns
            -------
            uint : Number of states.)V0G0N")
        .def("getRho", &PySimulation::getRho, R"V0G0N(
            Returns the electron density of the simulation.

            Returns
            -------
            float array : Electron density at each grid point.)V0G0N")
        .def("getCur", &PySimulation::getCur, R"V0G0N(
            Returns the current density of the simulation.

            Returns
            -------
            float array : Current density at each grid point.)V0G0N")
        .def("getPsi", &PySimulation::getPsi, R"V0G0N(
            Returns the wavefunction of the simulation.
            If the wavefunction is not yet set an error will be raised.

            Returns
            -------
            complex array : Wavefunction at each grid point for each state.
            The array is of size nPts * nStates, where nPts is the number of grid points and nStates is the number of states.)V0G0N")
        .def("getWeights", &PySimulation::getWeightValues, R"V0G0N(
            Returns the weights of the states in the simulation.
            If the weights are not yet calculated, they will be calculated using the current wavefunction.

            Returns
            -------
            float array : Weights of each state.)V0G0N")
        .def("getV", &PySimulation::getV, R"V0G0N(
            Returns the potential of the simulation.
            If the potential is not yet calculated, it will be calculated using the current density and wavefunction.
            If the wavefunction is not initialized, it will calculate the bare potential.

            Returns
            -------
            float array : Potential at each grid point.)V0G0N")
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
        .def("setPsi", &PySimulation::setPsi, R"V0G0N(
            Sets the wavefunction for the simulation.
            The wavefunction must be of size nPts, where nPts is the number of grid points.

            Parameters
            ----------
            psi : complex array
                Wavefunction to be set.
            nElec : uint
                Number of electrons in the system. Default is 1.   
            norm : NormalizationScheme
                Normalization scheme to use for the wavefunction. Default is UNNORMALIZED.
                If NORMALIZED, the wavefunction will be normalized such that the integral of |psi|^2 over the entire simulation space equals 1.
                If UNNORMALIZED, the wavefunction is unnormalized, typically for inhomogeneous/open systems.)V0G0N",
            "psi"_a, "nElec"_a = 1, "norm"_a = Densities::UNNORMALIZED)
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
        .def("findGroundState", &PySimulation::findGroundState, R"V0G0N(
            Finds the ground state of the system using the given max energy and max number of states.
            Nonlinear potentials assume a neutral charge distribution -- this function does not find a self-consistent solution.

            Parameters
            ----------
            maxStates : uint
                Maximum number of states to be found.
            emax : float
                Maximum energy of the eigenstates to be found.)V0G0N",
            "maxStates"_a, "emax"_a)
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
            Runs simulation using the appropriate method for the kinetic operator and potential.
            
            Parameters
            ----------
            nSteps : uint
                Number of steps to run the simulation for.
            scfIts : int, optional
                Number of self-consistent field iterations to perform. Default is 8.
                If 'scfIts = -m < 0' then a polynomial extrapolation of order (m-1) is used for implicitly needed potentials (Crank-Nicolson only).
                For other values, the SCF logic is detailed below.
            scfTol : float, optional
                Tolerance for the self-consistent field iterations. Default is 1e-6.
            
            SCF is only performed for nonlinear Crank-Nicolson calculations. The logic is as follows:
            - If `scfIts = 0` then no SCF iterations are performed.
            - If `scfTol > 0.0` (default behavior), SCF iterations continue until $\frac{\Delta t}{\hbar}\max_j{|V_j'-V_j|} < scfTol$ or if `scfIts` is reached.
            - If `scfTol = 0.0` and `scfIts = 0`, no SCF iterations are performed.
            - If `scfTol = 0.0` and `scfIts != 0`, SCF is performed for `scfIts` iterations.
            )V0G0N",
            "nSteps"_a, "scfIts"_a = 8, "scfTol"_a = 1e-6)
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
            "nSteps"_a, "scfIts"_a = 0, "scfTol"_a = 1e-6)
        .def("runCN_P", &PySimulation::runCN_P, R"V0G0N(
            Runs simulation using Crank-Nicolson method with polynomial extrapolation for potential estimation.)V0G0N",
            "nSteps"_a, "order"_a = 3)
        .def("setNumCallbackCalls", &PySimulation::setNumCallbackCalls, R"V0G0N(
            Sets the number of times throughout a run that the callback function will be called.

            Parameters
            ----------
            nCalls : int
                Number of times to call the callback function.)V0G0N",
            "nCalls"_a)
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