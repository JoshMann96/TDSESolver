#include "PySimulation.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"
#include <memory>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pytypes.h>

void init_Simulation(py::module &m) {

// WEIGHT CALCULATIONS

    py::class_<WfcToRho::Weight>(m, "Weight");

    py::class_<WfcToRho::BoundFermiGas, WfcToRho::Weight>(m, "BoundFermiGas")
        .def(py::init<double>(), R"V0G0N(
            Uses 3-D Fermi gas distribution at zero temperature to convert 1-D wavefunctions to an effective 3-D density.
            Mapping uses wavefunction initial eigenstates assuming they are bound and normalized.

            Parameters
            ----------
            ef : float
                Fermi energy.

            Returns
            -------
            BoundFermiGas)V0G0N",
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

    py::class_<WfcToRho::GaussianSmoothedDensity, WfcToRho::Density>(m, "GaussianSmoothedDensity")
        .def(py::init<double>(), R"V0G0N(
            Uses a Gaussian smoothing function to calculate the final density.
            Density = sum over states (weight x psi*psi) * Gaussian
                (* = convolution)

            Parameters
            ----------
            sigma : float
                Standard deviation of Gaussian. Typically the inverse of the Thomas-Fermi wavenumber.

            Returns
            -------
            GaussianSmoothedDensity)V0G0N",
            "sigma"_a);
    
    py::class_<WfcToRho::CylindricalDensity, WfcToRho::Density>(m, "CylindricalDensity")
        .def(py::init<double, double, double>(), R"V0G0N(
            Applies geometric dispersion to density assuming a cylindrical geometry (azimuthal symmetry) with some definite radius.


            Parameters
            ----------
            center : float
                Center of azimuthal symmetry.
            radius : float 
                Radius of azimuthal symmetry.
            minX : float
                Left boundary position.

            Returns
            -------
            CylindricalDensity)V0G0N",
            "center"_a, "radius"_a, "minX"_a);
            
// FDBCs

    py::class_<FDBCs::BoundaryCondition>(m, "BoundaryCondition");

    py::class_<FDBCs::CommonBC, FDBCs::BoundaryCondition>(m, "CommonBC");

    py::class_<FDBCs::TimeIndependentBC, FDBCs::CommonBC>(m, "TimeIndependentBC");

    // TODO: implement PyNeumannBC to gather correct dx (or do something like PSM_FreeElec below for brevity)
    //       also need to add enum Side to Python wrapper
    /*py::class_<FDBCs::NeumannBC, FDBCs::TimeIndependentBC>(m, "NeumannBC")
        .def(py::init<std::complex<double>>(), R"V0G0N(
            Neumann boundary condition.

            Parameters
            ----------
            coeff : complex
                Coefficient for Neumann BC.

            Returns
            -------
            NeumannBC)V0G0N",
            "coeff"_a);
    
    py::class_<FDBCs::DirichletBC, FDBCs::TimeIndependentBC>(m, "DirichletBC")
        .def(py::init<std::complex<double>>(), R"V0G0N(
            Dirichlet boundary condition.

            Parameters
            ----------
            coeff : complex
                Coefficient for Dirichlet BC.

            Returns
            -------
            DirichletBC)V0G0N",
            "coeff"_a);*/

    //TODO: Other BCs

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
    
    py::class_<KineticOperators::CrankNicolson, KineticOperators::KineticOperator_FDM>(m, "CrankNicolson")
        .def(py::init([](PySimulation* sim, double meff, FDBCs::BoundaryCondition* leftBC, FDBCs::BoundaryCondition* rightBC, bool useCuda){
            return std::unique_ptr<KineticOperators::CrankNicolson>(new KineticOperators::CrankNicolson(sim->getNumPoints(), sim->getDX(), sim->getDT(), meff, leftBC, rightBC, useCuda));
        }), py::keep_alive<1,4>(), py::keep_alive<1,5>(),
        R"V0G0N(
            Crank-Nicolson method using finite difference derivatives.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            meff : float
                Effective mass (1.0 = free electron).
            leftBC : BoundaryCondition
                Left boundary condition.
            rightBC : BoundaryCondition
                Right boundary condition.
            useCuda : bool
                Whether to use the CUDA solver for the tridiagonal system.

            Returns
            -------
            CrankNicolson)V0G0N",
            "sim"_a, "meff"_a, "leftBC"_a, "rightBC"_a, "useCuda"_a);

// SIMULATION

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
            int)V0G0N",
            "xp"_a)
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
        .def("addRightAbsBdy", &PySimulation::addRightAbsBdy, R"V0G0N(
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
            nElec : int
                Number of electrons in the system.
            energies : list
                Eigenstate energies. The boundary conditions must be consistent with these energies.)V0G0N",
            "nElec"_a, "energies"_a)
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
        .def("run", &PySimulation::run, R"V0G0N(
            Runs simulation using the appropriate method for the kinetic operator and potential.)V0G0N",
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
            int)V0G0N",
            "minPos"_a, "maxPos"_a);
}