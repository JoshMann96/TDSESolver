#include "PyMeasurers.h"

using namespace Measurers;

void init_Measurers(py::module &m) {
    py::class_<Measurer, std::unique_ptr<Measurer,py::nodelete>>(m, "Measurer");

    py::class_<DoubleConst, Measurer, std::unique_ptr<DoubleConst, py::nodelete>>(m, "Constant")
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
    
    py::class_<BasicMeasurers, Measurer, std::unique_ptr<BasicMeasurers, py::nodelete>>(m, "Basic")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<BasicMeasurers, py::nodelete>(new BasicMeasurers(
                sim->getNumPoints(), sim->getDX(), sim->getDT(), fol
            ));
        }), R"V0G0N(
            Records basic simulation parameters.
            Grid size and step size in space and time:
             - NPts     - Number of points in grid
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
    
    py::class_<XS, Measurer, std::unique_ptr<XS, py::nodelete>>(m, "Xs")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<XS, py::nodelete>(new XS(
                sim->getNumPoints(), sim->getX(), fol
            ));
        }), R"V0G0N(
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
    
    py::class_<TS, Measurer, std::unique_ptr<TS, py::nodelete>>(m, "Ts")
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

    py::class_<OrigPot, Measurer, std::unique_ptr<OrigPot, py::nodelete>>(m, "OrigPot")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<OrigPot, py::nodelete>(new OrigPot(
                sim->getNumPoints(), fol
            ));
        }), R"V0G0N(
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
    
    py::class_<NElec, Measurer, std::unique_ptr<NElec, py::nodelete>>(m, "NElec")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<NElec, py::nodelete>(new NElec(
                sim->getNElecPtr(), fol
            ));
        }), R"V0G0N(
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

    py::class_<WfcRhoWeights, Measurer, std::unique_ptr<WfcRhoWeights, py::nodelete>>(m, "Weights")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<WfcRhoWeights, py::nodelete>(new WfcRhoWeights(
                sim->getNElecPtr(), sim->getWeightsPtr(), fol
            ));
        }), R"V0G0N(
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

    py::class_<Psi2t, Measurer, std::unique_ptr<Psi2t, py::nodelete>>(m, "Psi2t")
        .def(py::init([](PySimulation* sim, size_t nx, size_t nt, size_t numSteps, std::string fol){
            return std::unique_ptr<Psi2t, py::nodelete>(new Psi2t(
                sim->getNumPoints(), nx, nt, numSteps, sim->getX(), sim->getNElecPtr(), fol
            ));
        }), R"V0G0N(
            Records wavefunction probability densities, downsampling to nx spatial points and nt temporal points.
            Parsing the output requires that Measurers NElec and Weights also be added to the simulation.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            nx : uint
                Number of spatial points to sample.
            nt : uint
                Number of temporal points to sample.
            numSteps : uint
                Number of time steps in the full calculation.
            fol : str
                Directory to contain file.

            Returns
            -------
            Psi2t)V0G0N",
            "sim"_a, "nx"_a, "nt"_a, "numSteps"_a, "fol"_a);

    py::class_<Vfunct, Measurer, std::unique_ptr<Vfunct, py::nodelete>>(m, "Vfunct")
        .def(py::init([](PySimulation* sim, size_t nx, size_t nt, size_t numSteps, std::string fol, int idx){
            return std::unique_ptr<Vfunct, py::nodelete>(new Vfunct(
                idx, sim->getNumPoints(), nx, nt, numSteps, numSteps*sim->getDT(), sim->getX(), fol
            ));
        }), R"V0G0N(
            Records potential, downsampling to nx spatial points and nt temporal points.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            nx : uint
                Number of spatial points to sample.
            nt : uint
                Number of temporal points to sample.
            numSteps : uint
                Number of time steps in the full calculation.
            fol : str
                Directory to contain file.
            idx: int
                Index of record, to prepend the output file. < 0 for no index, >= 0 for listed index. Defaults to -1.

            Returns
            -------
            Vfunct)V0G0N",
            "sim"_a, "nx"_a, "nt"_a, "numSteps"_a, "fol"_a, "idx"_a=-1);

    py::class_<ExpectE, Measurer, std::unique_ptr<ExpectE, py::nodelete>>(m, "ExpectE")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<ExpectE, py::nodelete>(new ExpectE(
                sim->getNumPoints(), sim->getDX(), sim->getNElecPtr(), fol, sim->getKin()
            ));
        }), R"V0G0N(
            Records expectation value of Hamiltonian at each time step.
            
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
    
    py::class_<ExpectE0, Measurer, std::unique_ptr<ExpectE0, py::nodelete>>(m, "ExpectE0")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<ExpectE0, py::nodelete>(new ExpectE0(
                sim->getNumPoints(), sim->getDX(), sim->getNElecPtr(), fol, sim->getKin()
            ));
        }), R"V0G0N(
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
    
    py::class_<ExpectX, Measurer, std::unique_ptr<ExpectX, py::nodelete>>(m, "ExpectX")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<ExpectX, py::nodelete>(new ExpectX(
                sim->getNumPoints(), sim->getX(), sim->getDX(), sim->getNElecPtr(), fol
            ));
        }), R"V0G0N(
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

    py::class_<ExpectP, Measurer, std::unique_ptr<ExpectP, py::nodelete>>(m, "ExpectP")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<ExpectP, py::nodelete>(new ExpectP(
                sim->getNumPoints(), sim->getDX(), sim->getNElecPtr(), fol
            ));
        }), R"V0G0N(
            Records expectation value of momentum for each state. Note: computationally expensive.

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

    py::class_<ExpectA, Measurer, std::unique_ptr<ExpectA, py::nodelete>>(m, "ExpectA")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<ExpectA, py::nodelete>(new ExpectA(
                sim->getNumPoints(), sim->getDX(), sim->getNElecPtr(), fol
            ));
        }), R"V0G0N(
            Records expectation value of acceleration for each state.
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

    py::class_<TotProb, Measurer, std::unique_ptr<TotProb, py::nodelete>>(m, "TotProb")
        .def(py::init([](PySimulation* sim, std::string fol){
            return std::unique_ptr<TotProb, py::nodelete>(new TotProb(
                sim->getNumPoints(), sim->getDX(), sim->getNElecPtr(), fol
            ));
        }), R"V0G0N(
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

    py::class_<VDProbCurrent, Measurer, std::unique_ptr<VDProbCurrent, py::nodelete>>(m, "VDProbCurrent")
        .def(py::init([](PySimulation* sim, double vdPos, int vdNum, std::string name, std::string fol){
            return std::unique_ptr<VDProbCurrent, py::nodelete>(new VDProbCurrent(
                sim->getNumPoints(), sim->getDX(), sim->getNElecPtr(), sim->findXIdx(vdPos), vdNum, name, fol
            ));
        }), R"V0G0N(
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

    py::class_<VDPsi, Measurer, std::unique_ptr<VDPsi, py::nodelete>>(m, "VDPsi")
        .def(py::init([](PySimulation* sim, double vdPos, int vdNum, std::string name, std::string fol){
            return std::unique_ptr<VDPsi, py::nodelete>(new VDPsi(
                sim->getNElecPtr(), sim->findXIdx(vdPos), vdNum, name, fol
            ));
        }), R"V0G0N(
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

    py::class_<VDPot, Measurer, std::unique_ptr<VDPot, py::nodelete>>(m, "VDPot")
        .def(py::init([](PySimulation* sim, double vdPos, int vdNum, std::string name, std::string fol){
            return std::unique_ptr<VDPot, py::nodelete>(new VDPot(
                sim->findXIdx(vdPos), vdNum, name, fol
            ));
        }), R"V0G0N(
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

    py::class_<VDFluxSpec, Measurer, std::unique_ptr<VDFluxSpec, py::nodelete>>(m, "VDFluxSpec")
        .def(py::init([](PySimulation* sim, double vdPos, int vdNum, size_t nSamp, double emax, double maxT, std::string name, std::string fol){
            return std::unique_ptr<VDFluxSpec, py::nodelete>(new VDFluxSpec(
                sim->getNumPoints(), sim->findXIdx(vdPos), vdNum, sim->getNElecPtr(), nSamp, emax, maxT, name, fol
            ));
        }), R"V0G0N(
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
            nSamp : uint
                Number of energy samples.
            emax : float
                Maximum energy of spectrum.
            maxT : float
                Maximum time of simulation.
            name : str
                Name of virtual detector (4 characters)
            fol : str
                Directory to contain file.

            Returns
            -------
            VDFluxSpec)V0G0N",
            "sim"_a, "vdPos"_a, "vdNum"_a, "nSamp"_a, "emax"_a, "maxT"_a, "name"_a, "fol"_a);

    py::class_<PsiT, Measurer, std::unique_ptr<PsiT, py::nodelete>>(m, "PsiT")
        .def(py::init([](PySimulation* sim, double meaT, int vdNum, std::string name, std::string fol){
            return std::unique_ptr<PsiT, py::nodelete>(new PsiT(
                sim->getNumPoints(), meaT, sim->getNElecPtr(), vdNum, name, fol
            ));
        }), R"V0G0N(
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

    py::class_<PotT, Measurer, std::unique_ptr<PotT, py::nodelete>>(m, "PotT")
        .def(py::init([](PySimulation* sim, double meaT, int vdNum, std::string name, std::string fol){
            return std::unique_ptr<PotT, py::nodelete>(new PotT(
                sim->getNumPoints(), meaT, vdNum, name, fol
            ));
        }), R"V0G0N(
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
