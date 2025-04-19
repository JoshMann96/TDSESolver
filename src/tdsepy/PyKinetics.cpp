#include "PyKinetics.h"
#include "KineticOperator.h"
#include "WfcRhoTools.h"
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pytypes.h>
#include "PySimulation.h"

void init_Kinetics(py::module &m) {
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
        .def(py::init([](PySimulation* sim, size_t order, size_t nElec, 
            py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> psibd, 
            py::array_t<double, py::array::c_style | py::array::forcecast> energies, 
            double m_eff, double vb){
                
            vtlsPrnt::printArray(nElec, energies.data());

            // get k0
            double* k0s = (double*) sq_malloc(sizeof(double) * nElec);
            for (size_t i = 0; i < nElec; i++)
                k0s[i] = KineticOperators::CrankNicolson::wavenumberFromEnergy(energies.at(i), vb, sim->getDX(), sim->getDT(), m_eff);

            std::unique_ptr<FDBCs::UniformIDTransparentBC> res(new FDBCs::UniformIDTransparentBC(order, nElec, sim->getDX(), sim->getDT(), psibd.data(), k0s, vb));

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
            vb : float
                Potential at the boundary.
            
            Returns
            -------
            UniformIDTransparentBC)V0G0N",
            "sim"_a, "order"_a, "nElec"_a, "psibd"_a, "energies"_a, "m_eff"_a, "vb"_a);
    

// KINETIC OPERATORS

    py::class_<KineticOperators::KineticOperator>(m, "KineticOperator");

    py::class_<KineticOperators::KineticOperator_PSM, KineticOperators::KineticOperator>(m, "KineticOperator_PSM");

    py::class_<KineticOperators::KineticOperator_FDM, KineticOperators::KineticOperator>(m, "KineticOperator_FDM")
        .def("setBC", &KineticOperators::KineticOperator_FDM::setBC, py::keep_alive<1,2>(), R"V0G0N(
            Sets the left or right boundary condition for the system.
            This may be done at any time.

            Parameters
            ----------
            bc : BoundaryCondition
                The boundary condition to set.
            side : BCSide
                The side of the system to set the boundary condition for (left or right).)V0G0N",
            "bc"_a, "side"_a);

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
    }