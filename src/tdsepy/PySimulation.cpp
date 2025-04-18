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

    py::class_<WfcToRho::UniformWeight, WfcToRho::Weight>(m, "UniformWeight")
        .def(py::init<double>(), R"V0G0N(
            Sets the weights to be a constant value for all states.

            Parameters
            ----------
            weight : float
                Weight to be assigned to each state.

            Returns
            -------
            UniformWeight)V0G0N",
            "weight"_a);

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
    
    py::class_<WfcToRho::SemiInfiniteFermiGas, WfcToRho::Weight>(m, "SemiInfiniteFermiGas")
        .def(py::init<double>(), R"V0G0N(
            Uses 3-D Fermi gas distribution at zero temperature to convert 1-D wavefunctions to an effective 3-D density.
            Mapping uses wavefunction initial eigenstates assuming they are eigenstates of the open system, in contact with a zero-temperature free electron gas (FEG).
            The boundary conditions on the FEG side must be inhomogeneous and have an incoming component of magnitude 1.

            Parameters
            ----------
            ef : float
                Fermi energy.

            Returns
            -------
            SemiInfiniteFermiGas)V0G0N",
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

    py::enum_<FDBCs::BCSide>(m, "BCSide")
        .value("LEFT", FDBCs::BCSide::LEFT)
        .value("RIGHT", FDBCs::BCSide::RIGHT);

    py::class_<FDBCs::BoundaryCondition>(m, "BoundaryCondition");

    py::class_<FDBCs::CommonBC, FDBCs::BoundaryCondition>(m, "CommonBC");

    py::class_<FDBCs::TimeIndependentBC, FDBCs::CommonBC>(m, "TimeIndependentBC");

    py::class_<FDBCs::NeumannBC, FDBCs::TimeIndependentBC>(m, "NeumannBC")
        .def(py::init([](PySimulation* sim, std::complex<double> bdDer, FDBCs::BCSide side){
            return std::unique_ptr<FDBCs::NeumannBC>(new FDBCs::NeumannBC(bdDer, sim->getDX(), side));
        }), R"V0G0N(
            Neumann boundary condition.

            Parameters
            ----------
            sim : Simulation
                Associated simulation (only used for grabbing dx, ownership is not given to the Simulation).
            bdDer : complex
                Derivative for Neumann BC.
            side : BCSide
                Side of the boundary condition. Determines the sign of the derivative.

            Returns
            -------
            NeumannBC)V0G0N",
            "sim"_a, "bdDer"_a, "side"_a);
    
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
            "coeff"_a);

    py::class_<FDBCs::UniformHDTransparentBC, FDBCs::BoundaryCondition>(m, "UniformHDTransparentBC")
        .def(py::init([](PySimulation* sim, size_t order, size_t nElec){
            return std::unique_ptr<FDBCs::UniformHDTransparentBC>(new FDBCs::UniformHDTransparentBC(order, nElec, sim->getDX(), sim->getDT()));
        }), R"V0G0N(
            Uniform homogeneous discrete transparent boundary condition.
                The grid the simulation is on must be uniform.
                The potential must be constant at the boundary.
                The boundary condition is homogeneous -- it has no sources.
            For more information, see https://doi.org/10.4310/CMS.2003.v1.n3.a7

            Parameters
            ----------
            sim : Simulation
                Associated simulation (only used for grabbing dx, ownership is not given to the Simulation).
            order : uint
                Order of the boundary condition. This determines the length of the truncated history to use.
                It is recommended to have this be much larger than the number of time steps in the period associated with the loest relevant kinetic energy:
                order >> h / ( dt * E_min )
            nElec : uint
                Number of electrons. This must be known at the time of construction, so either the eigenstates must be found first with Dirichlet BCs 
                or an inhomogeneous BC (with a set number of wavefunctions) must be additionally used.

            Returns
            -------
            UniformHDTransparentBC)V0G0N",
            "sim"_a, "order"_a, "nElec"_a);

    py::class_<FDBCs::UniformIDTransparentBC, FDBCs::UniformHDTransparentBC>(m, "UniformIDTransparentBC")
        .def(py::init([](PySimulation* sim, size_t order, size_t nElec, std::complex<double>* psibd, double* energies, double m_eff, FDBCs::BCSide side){
            // get boundary potential
            double* v = (double*) sq_malloc(sizeof(double) * nElec);
            sim->getPotPointer()->getVBare(0.0, v);
            double vb;
            if (side == FDBCs::BCSide::LEFT)
                vb = v[0];
            else
                vb = v[sim->getNumPoints()-1];
            sq_free(v);

            // get k0
            double* k0s = (double*) sq_malloc(sizeof(double) * nElec);
            for (size_t i = 0; i < nElec; i++)
                k0s[i] = KineticOperators::CrankNicolson::wavenumberFromEnergy(energies[i], vb, sim->getDX(), sim->getDT(), m_eff);

            std::unique_ptr<FDBCs::UniformIDTransparentBC> res(new FDBCs::UniformIDTransparentBC(order, nElec, sim->getDX(), sim->getDT(), psibd, k0s, vb));

            sq_free(k0s);

            return res;
        }), R"V0G0N(
            Uniform inhomogeneous discrete transparent boundary condition.
                The grid the simulation is on must be uniform.
                The boundary condition is inhomogeneous -- it has a source of the form psibd e^{i k0 x} phi(k0)^(t/dt)
                    phi(k0) is the one-step complex phase advance representing the discrete dispersion relation of the Crank-Nicolson method (in Hartree atomic units): 
                        phi(k) = 1 - 0.5i * dt ( (1 - cos(k * dx)) / (dx^2 * m_eff) + vb )
            For more information, see https://doi.org/10.4310/CMS.2003.v1.n3.a7

            Construction calls Potential::getVBare to get the potential at the boundary. All initially non-zero potentials should be added to the Simulation prior to constructing this.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            order : uint
                Order of the boundary condition. This determines the length of the truncated history to use.
                It is recommended to have this be much larger than the number of time steps in the period associated with the loest relevant kinetic energy:
                order >> h / ( dt * E_min )
            nElec : uint
                Number of electrons.
            psibd : complex array
                Boundary wavefunction, nElec elements.
                It is recommended that the magnitude of the values are 1.0.
            energies : float array
                Energies of the eigenstates, nElec elements.
            m_eff : float
                Effective mass of the electrons. 1.0 = free electron mass.
            side : BCSide
                Side of the boundary condition. Used for finding the initial potential.
            
            Returns
            -------
            UniformIDTransparentBC)V0G0N",
            "sim"_a, "order"_a, "nElec"_a, "psibd"_a, "energies"_a, "m_eff"_a, "side"_a);
    

// KINETIC OPERATORS

    py::class_<KineticOperators::KineticOperator>(m, "KineticOperator");

    py::class_<KineticOperators::KineticOperator_PSM, KineticOperators::KineticOperator>(m, "KineticOperator_PSM");

    py::class_<KineticOperators::KineticOperator_FDM, KineticOperators::KineticOperator>(m, "KineticOperator_FDM");

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
            uint)V0G0N",
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
            energies : list
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