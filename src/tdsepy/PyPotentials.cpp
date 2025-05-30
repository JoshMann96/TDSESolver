#include "PyPotentials.h"
#include "Potentials.h"
#include "Measurers.h"
#include "PySimulation.h"

/*std::vector<double> PyFilePotential::getValue(){
    double* temp = new double[nPts];
    FilePotential::getV(0.0, temp, nullptr);
    std::vector out(temp, temp+nPts);
    return out;
}*/

using namespace Potentials;

void init_Potentials(py::module &m) {

// Envelopes

    py::class_<Envelopes::Envelope>(m, "Envelope");

    py::class_<Envelopes::SmoothedInitialGaussianEnvelope, Envelopes::Envelope>(m, "GaussianEnvelope")
        .def(py::init<double, double, double>(), R"V0G0N(
            Gaussian laser envelope.

            Parameters
            ----------
            tau : float
                Full-width half-max power.
            tmax : float
                Envelope center in time.
            bufferTime: float
                Time length of polynomial-smoothing factor such that the profile is zero at t=0.

            Returns
            -------
            GaussianEnvelope)V0G0N",
            "tau"_a, "tmax"_a, "bufferTime"_a);

    py::class_<Envelopes::CosSquaredEnvelope, Envelopes::Envelope>(m, "CosSquaredEnvelope")
        .def(py::init<double, double>(), R"V0G0N(
            Cosine-squared laser envelope.

            Parameters
            ----------
            tau : float
                Full-width half-max power.
            tmax : float
                Envelope center in time.

            Returns
            -------
            CosSquaredEnvelope)V0G0N",
            "tau"_a, "tmax"_a);

// Field Profiles

    py::class_<ElectricFieldProfiles::ElectricFieldProfile>(m, "FieldProfile");

    py::class_<ElectricFieldProfiles::FileFieldProfile, ElectricFieldProfiles::ElectricFieldProfile>(m, "FileFieldProfile")
        .def(py::init([](SimulationManager* sim, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, std::string fil){
            return std::unique_ptr<ElectricFieldProfiles::FileFieldProfile>(new ElectricFieldProfiles::FileFieldProfile(
                sim->getNumPoints(), sim->getX(), offset, rightDecayPos, leftDecayPos, decayLength, emax, fil
            ));
        }), R"V0G0N(
            Spatial laser field profile as defined in a file.
            See documentation for appropriate file format.
            The field may be further confined by the "decay" parameters.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            offset : float
                Translational offset with respect to input data.
            rightDecayPos : float
                Right-side position to begin decaying field to zero.
            leftDecayPos : float
                Left-side position to begin decaying field to zero.
            decayLength : float
                Lengthscale over which field is decayed.
            emax : float
                Maximum field strength.
            fil : str
                File path with data.

            Returns
            -------
            FileFieldProfile)V0G0N",
            "sim"_a, "offset"_a, "rightDecayPos"_a, "leftDecayPos"_a, "decayLength"_a, "emax"_a, "fil"_a);

// Potentials

    py::class_<Potential>(m, "Potential");

    py::class_<FilePotential, Potential>(m, "FilePotential")
        .def(py::init([](SimulationManager* sim, double offset, const std::string fil, double refPoint){
            return std::unique_ptr<FilePotential>(new FilePotential(
                sim->getNumPoints(), sim->getX(), offset, fil, sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Static potential as defined in a binary file.
            See documentation for appropriate file format.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            offset : float
                Translational offset with respect to input data.
            fil : str
                File path with data.
            refPoint : float
                Potential reference point.

            Returns
            -------
            FilePotential)V0G0N",
            "sim"_a, "offset"_a, "fil"_a, "refPoint"_a);


    py::class_<JelliumPotentialBacked, Potential>(m, "JelliumPotentialBacked")
        .def(py::init([](SimulationManager* sim, double center, double ef, double w, double backStart, double backWidth, double refPoint){
            return std::unique_ptr<JelliumPotentialBacked>(new JelliumPotentialBacked(
                sim->getNumPoints(), sim->getX(), center, ef, w, backStart, backWidth, sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Static Jellium slab potential. The right side of the potential is the Jellium surface while the left side is a polynomial-smooth backing.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            center : float
                Center point of surface sigmoid function.
            ef : float
                Fermi energy.
            w : float
                Work function.
            backStart : float
                Start of polynomial-smooth rear potential.
            backWidth : float
                Width of polynomial-smooth rear potential.
            refPoint : float
                Potential reference point.

            Returns
            -------
            JelliumPotentialBacked)V0G0N",
            "sim"_a, "center"_a, "ef"_a, "w"_a, "backStart"_a, "backWidth"_a, "refPoint"_a);

    py::class_<JelliumPotential, Potential>(m, "JelliumPotential")
        .def(py::init([](SimulationManager* sim, double center, double ef, double w, double refPoint){
            return std::unique_ptr<JelliumPotential>(new JelliumPotential(
                sim->getNumPoints(), sim->getX(), center, ef, w, sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Static semiinfinite Jellium potential. The right side of the potential is the Jellium surface, and to the left is within the material.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            center : float
                Center point of surface sigmoid function.
            ef : float
                Fermi energy.
            w : float
                Work function.
            refPoint : float
                Potential reference point.

            Returns
            -------
            JelliumPotentialBacked)V0G0N",
            "sim"_a, "center"_a, "ef"_a, "w"_a, "refPoint"_a);

    py::class_<FiniteBox, Potential>(m, "FiniteBox")
        .def(py::init([](SimulationManager* sim, double left, double right, double vin, double refPoint){
            return std::unique_ptr<FiniteBox>(new FiniteBox(
                sim->getNumPoints(), sim->getX(), left, right, vin, sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Static finite box potential.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            left : float
                Left-side position of box.
            right : float
                Right-side position of box.
            vin : float
                Potential inside the box.
            refPoint : float
                Potential reference point.

            Returns
            -------
            FiniteBox)V0G0N",
            "sim"_a, "left"_a, "right"_a, "vin"_a, "refPoint"_a);

    py::class_<ElectricFieldProfileToPotential, Potential>(m, "PulsePotential")
        .def(py::init([](SimulationManager* sim, ElectricFieldProfiles::ElectricFieldProfile* fieldProfile, Envelopes::Envelope* env, double phase, double tmax, double lam, double refPoint){
            return std::unique_ptr<ElectricFieldProfileToPotential>(new ElectricFieldProfileToPotential(
                sim->getNumPoints(), fieldProfile, sim->getDX(), phase, tmax, lam, env, sim->findXIdx(refPoint)
            ));
        }),  py::keep_alive<1,3>(),  py::keep_alive<1,4>(), R"V0G0N(
            Pulsed laser potential under dipole approximation. 

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            fieldProfile : FieldProfile
                Electric field profile of laser.
            env : Envelope
                Laser temporal envelope.
            phase : float
                Carrier-envelope phase (CEP), measured with respect to tmax.
            tmax : float
                Center (maximum) of envelope for CEP.
            lam : float
                Laser wavelength.
            refPoint : float
                Potential reference point.

            Returns
            -------
            PulsePotential)V0G0N",
            "sim"_a, "fieldProfile"_a, "env"_a, "phase"_a, "tmax"_a, "lam"_a, "refPoint"_a);
    

    py::class_<CylindricalImageCharge, Potential>(m, "CylindricalImagePotential")
        .def(py::init([](SimulationManager* sim, double ef, double w, double rad, double posMin, double posMax, double surfPos, double refPoint){
            return std::unique_ptr<CylindricalImageCharge>(new CylindricalImageCharge(
                sim->getNumPoints(), sim->getX(), sim->getDX(), ef, w, rad, sim->findXIdx(surfPos), sim->getNElecPtr(), sim->getWeightsPtr(), sim->getRho(), sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Collective image charge potential assuming a cylindrical conductor geometry.
            If the wavefunction has not been initialized upon construction then the potential will be returned as-is (with reference to refPoint).
            If the wavefunction has been initialized, then the potential will be returned as the change in potential with respect to the initial density.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            ef : float
                Fermi energy (for Jellium-like surface mask).
            w : float
                Work function (for Jellium-like surface mask).
            rad : float
                Cylinder radius of curvature.
            posMin : float
                Minimum pos where density is included.
            posMax : float
                Maximum pos where density is included.
            surfPos : float
                Surface position (index) of the cylinder.
            refPoint : float
                Potential reference point.

            Returns
            -------
            CylindricalImagePotential)V0G0N",
            "sim"_a, "ef"_a, "w"_a, "rad"_a, "posMin"_a, "posMax"_a, "surfPos"_a, "refPoint"_a);

    py::class_<PlanarToCylindricalHartree,Potential>(m, "PlanarToCylindricalHartreePotential")
        .def(py::init([](SimulationManager* sim, bool mimickOpenSystem, double rad, double posMin, double posMax, double surfPos, double refPoint){
            return std::unique_ptr<PlanarToCylindricalHartree>(new PlanarToCylindricalHartree(
                mimickOpenSystem, sim->getNumPoints(), sim->getDX(), rad, sim->findXIdx(surfPos), sim->getNElecPtr(), sim->getWeightsPtr(), sim->wavefunctionIsInitialized() ? sim->getRho() : nullptr, sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Nonlocal Hartree potential assuming charge is distributed on a planar geometry for x <= surfPos 
            and on a cylindrical geometry for x > surfPos, with a transition radius of curvature rad (the planar
            charge is assumed to be within the cylinder of radius rad for x > surfPos).
            Charge lost to the left (planar) boundary is re-distributed over the initial density.
            surfPos is set after construction via assemble.
            If the sim's wavefunction has not been initialized upon construction then the potential will be returned as-is (with reference to refPoint).
            If the sim's wavefunction has been initialized, then the potential will be returned as the change in potential with respect to the initial density.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            mimickOpenSystem : bool
                If true, the charge is scaled to conserve the total charge, minus what leaves the left boundary at posMax.
                This can only be used if the density/wavefunction is initialized.
            rad : float
                Cylinder radius of curvature.
            posMin : float
                Minimum pos where density is included.
            posMax : float
                Maximum pos where density is included.
            surfPos : float
                Surface position (index) of the cylinder.
            refPoint : float
                Potential reference point.

            Returns
            -------
            PlanarToCylindricalHartreePotential)V0G0N",
            "sim"_a, "mimickOpenSystem"_a, "rad"_a, "posMin"_a, "posMax"_a, "surfPos"_a, "refPoint"_a);

    py::class_<LDAFunctional, Potential>(m, "LDAFunctional")
        .def(py::init([](SimulationManager* sim, LDAFunctionalType typ, double refPoint){
            return std::unique_ptr<LDAFunctional>(new LDAFunctional(
                typ, sim->getNumPoints(), sim->getDX(), sim->wavefunctionIsInitialized() ? sim->getRho() : nullptr, sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Local density approximation (LDA) functional potential.
            If the wavefunction has not been initialized upon construction then the potential will be returned as-is (with reference to refPoint).
            If the wavefunction has been initialized, then the potential will be returned as the change in potential with respect to the initial density.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            typ : LDAFunctionalType
                Type of LDA functional. Choices include:
                    X_SLATER : Slater exchange.
                    C_PW : Perdew-Wang correlation.
            refPoint : float
                Potential reference point.

            Returns
            -------
            LDAFunctional)V0G0N",
            "sim"_a, "typ"_a, "refPoint"_a);
    
    py::class_<MeasuredPotential, Potential>(m, "MeasuredPotential")
        .def(py::init([](SimulationManager* sim, Potential* pot, Measurers::Measurer* meas, size_t numSteps){
            return std::unique_ptr<MeasuredPotential>(new MeasuredPotential(
                pot, meas, numSteps, numSteps*sim->getDT()
            ));
        }), py::keep_alive<1,3>(), py::keep_alive<1,4>(), R"V0G0N(
            A potential which is also measured when it is called.
            At each evaluation the potential is calculated and then the measurer passed is called using that potential only.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            pot : Potential
                Potential to use and measure.
            meas: Measurer
                Measurer to use.
            numSteps : uint
                Number of time steps to measure.

            Returns
            -------
            MeasuredPotential)V0G0N",
            "sim"_a, "pot"_a, "meas"_a, "numSteps"_a);
            
    py::enum_<LDAFunctionalType>(m, "LDAFunctionalType")
        .value("X_SLATER", LDAFunctionalType::X_SLATER)
        .value("C_PW", LDAFunctionalType::C_PW)
        .export_values();
}