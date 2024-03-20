#include "PyMeasurers.h"

using namespace Measurers;

void init_Measurers(py::module &m) {
    py::class_<Measurer>(m, "Measurer");

    py::class_<PyConstant, Measurer>(m, "Constant")
        .def(py::init<double, std::string, std::string>(), R"V0G0N(
            Records a constant.

            Parameters
            ----------
            c : float
                Value to be recorded.
            fileName : str
                Name of output file (no extension).
            fol : str
                Directory to contain file.

            Returns
            -------
            Constant)V0G0N",
            "c"_a, "fileName"_a, "fol"_a);
    
    py::class_<PyBasic, Measurer>(m, "Basic")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records basic simulation parameters.
            Grid size and step size in space and time:
             - nPts     - Number of points in grid
             - NSteps   - Number of time steps
             - DX       - Spatial step size
             - DT       - Temporal step size

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            Basic)V0G0N",
            "sim"_a, "fol"_a);
    
    py::class_<PyXS, Measurer>(m, "Xs")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records spatial grid.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            Xs)V0G0N",
            "sim"_a, "fol"_a);
    
    py::class_<PyTS, Measurer>(m, "Ts")
        .def(py::init<std::string>(), R"V0G0N(
            Records time steps.

            Parameters
            ----------
            fol : str
                Directory to contain file.

            Returns
            -------
            Ts)V0G0N",
            "fol"_a);
    
    py::class_<PyNElec, Measurer>(m, "NElec")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records the number of states in the simulation.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            NElec)V0G0N",
            "sim"_a, "fol"_a);

    py::class_<PyWeights, Measurer>(m, "Weights")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records the weights for each state to map from 1-D to 3-D densities.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            Weights)V0G0N",
            "sim"_a, "fol"_a);

    py::class_<PyOrigPot, Measurer>(m, "OrigPot")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records initial potential.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            OrigPot)V0G0N",
            "sim"_a, "fol"_a);

    py::class_<PyPsi2t, Measurer>(m, "Psi2t")
        .def(py::init<PySimulation*, int, int, std::string>(), R"V0G0N(
            Records wavefunction probability densities, downsampling to nx spatial points and nt temporal points.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            nx : int
                Number of spatial points to sample.
            nt : int
                Number of temporal points to sample.
            fol : str
                Directory to contain file.

            Returns
            -------
            Psi2t)V0G0N",
            "sim"_a, "nx"_a, "nt"_a, "fol"_a);

    py::class_<PyVfunct, Measurer>(m, "Vfunct")
        .def(py::init<PySimulation*, int, int, std::string>(), R"V0G0N(
            Records potential, downsampling to nx spatial points and nt temporal points.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            nx : int
                Number of spatial points to sample.
            nt : int
                Number of temporal points to sample.
            fol : str
                Directory to contain file.

            Returns
            -------
            Vfunct)V0G0N",
            "sim"_a, "nx"_a, "nt"_a, "fol"_a);

    py::class_<PyExpectE, Measurer>(m, "ExpectE")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records expectation value of Hamiltonian at each time step.
            --- WARNING ---
            Uses free-electron 3-point stencil for kinetic energy evaluation.
            This should be changed to use the KineticOperator
            See implementation of ExpectE0
            ---------------

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            ExpectE)V0G0N",
            "sim"_a, "fol"_a);
    
    py::class_<PyExpectE0, Measurer>(m, "ExpectE0")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records expectation value of Hamiltonian at first time step.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            ExpectE0)V0G0N",
            "sim"_a, "fol"_a);
    
    py::class_<PyExpectX, Measurer>(m, "ExpectX")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records expectation value of position for each state.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            ExpectX)V0G0N",
            "sim"_a, "fol"_a);

    py::class_<PyExpectP, Measurer>(m, "ExpectP")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records expectation value of momentum for each state. Note: computationally expective.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            ExpectP)V0G0N",
            "sim"_a, "fol"_a);

    py::class_<PyExpectA, Measurer>(m, "ExpectA")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records expectation value of acceleration for each tate.
            More precisely, the negative gradient of the potential divided by the electron mass.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            ExpectA)V0G0N",
            "sim"_a, "fol"_a);

    py::class_<PyTotProb, Measurer>(m, "TotProb")
        .def(py::init<PySimulation*, std::string>(), R"V0G0N(
            Records integrated probability of each state.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fol : str
                Directory to contain file.

            Returns
            -------
            TotProb)V0G0N",
            "sim"_a, "fol"_a);

    py::class_<PyVDProbCurrent, Measurer>(m, "VDProbCurrent")
        .def(py::init<PySimulation*, double, int, std::string, std::string>(), R"V0G0N(
            Virtual detector which records probability current at a position over time for each state.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            vdPos : float
                Position of virtual detector.
            vdNum : int
                Index of virtual detector.
            name : str
                Name of virtual detector (4 characters)
            fol : str
                Directory to contain file.

            Returns
            -------
            VDProbCurrent)V0G0N",
            "sim"_a, "vdPos"_a, "vdNum"_a, "name"_a, "fol"_a);

    py::class_<PyVDPsi, Measurer>(m, "VDPsi")
        .def(py::init<PySimulation*, double, int, std::string, std::string>(), R"V0G0N(
            Virtual detector which records the wavefunction's complex value at a position over time for each state.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            vdPos : float
                Position of virtual detector.
            vdNum : int
                Index of virtual detector.
            name : str
                Name of virtual detector (4 characters)
            fol : str
                Directory to contain file.

            Returns
            -------
            VDPsi)V0G0N",
            "sim"_a, "vdPos"_a, "vdNum"_a, "name"_a, "fol"_a);

    py::class_<PyVDPot, Measurer>(m, "VDPot")
        .def(py::init<PySimulation*, double, int, std::string, std::string>(), R"V0G0N(
            Virtual detector which records the potential at a position over time for each state.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            vdPos : float
                Position of virtual detector.
            vdNum : int
                Index of virtual detector.
            name : str
                Name of virtual detector (4 characters)
            fol : str
                Directory to contain file.

            Returns
            -------
            VDPot)V0G0N",
            "sim"_a, "vdPos"_a, "vdNum"_a, "name"_a, "fol"_a);

    py::class_<PyVDFluxSpec, Measurer>(m, "VDPsi")
        .def(py::init<PySimulation*, double, int, int, double, std::string, std::string>(), R"V0G0N(
            Virtual detector which measures the bidirectional flux spectrum of the state passing through a point for each state.
            Useful for obtaining electron emission spectra without saving the entire wavefunction history.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            vdPos : float
                Position of virtual detector.
            vdNum : int
                Index of virtual detector.
            nSamp : int
                Number of energy samples.
            emax : float
                Maximum energy of spectrum.
            name : str
                Name of virtual detector (4 characters)
            fol : str
                Directory to contain file.

            Returns
            -------
            VDPsi)V0G0N",
            "sim"_a, "vdPos"_a, "vdNum"_a, "nSamp"_a, "emax"_a, "name"_a, "fol"_a);

    py::class_<PyPsiT, Measurer>(m, "PsiT")
        .def(py::init<PySimulation*, double, int, std::string, std::string>(), R"V0G0N(
            Virtual detector which records the wavefunctions at a set time.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            meaT : float
                Measurement time.
            vdNum : int
                Index of virtual detector.
            name : str
                Name of virtual detector (4 characters)
            fol : str
                Directory to contain file.

            Returns
            -------
            PsiT)V0G0N",
            "sim"_a, "meaT"_a, "vdNum"_a, "name"_a, "fol"_a);

    py::class_<PyPotT, Measurer>(m, "PotT")
        .def(py::init<PySimulation*, double, int, std::string, std::string>(), R"V0G0N(
            Virtual detector which records the potential at a set time.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            meaT : float
                Measurement time.
            vdNum : int
                Index of virtual detector.
            name : str
                Name of virtual detector (4 characters)
            fol : str
                Directory to contain file.

            Returns
            -------
            PotT)V0G0N",
            "sim"_a, "meaT"_a, "vdNum"_a, "name"_a, "fol"_a);
}
