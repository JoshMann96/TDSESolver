#include "PyPotentials.h"
#include "Potentials.h"
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

    py::class_<PyFileFieldProfile, ElectricFieldProfiles::ElectricFieldProfile>(m, "FileFieldProfile")
        .def(py::init<PySimulation*, double, double, double, double, double, std::string>(), R"V0G0N(
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
                Lengscale over which field is decayed.
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

    py::class_<PyFilePotential, Potential>(m, "FilePotential")
        .def(py::init<PySimulation*, double, const std::string, double>(), R"V0G0N(
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


    py::class_<PyJelliumPotential, Potential>(m, "JelliumPotential")
        .def(py::init<PySimulation*, double, double, double, double, double, double>(), R"V0G0N(
            Static Jellium slab potential.

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
            JelliumPotential)V0G0N",
            "sim"_a, "center"_a, "ef"_a, "w"_a, "backStart"_a, "backWidth"_a, "refPoint"_a);

    
    py::class_<PyPulsePotential, Potential>(m, "PulsePotential")
        .def(py::init<PySimulation*, ElectricFieldProfiles::ElectricFieldProfile*, Envelopes::Envelope*, double, double, double, double>(),  py::keep_alive<1,3>(),  py::keep_alive<1,4>(), R"V0G0N(
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
    
    py::class_<PyCylindricalImagePotential, Potential>(m, "CylindricalImagePotential")
        .def(py::init<PySimulation*, double, double, double, double, double, double, double>(), R"V0G0N(
            Collective image charge potential assuming a cylindrical conductor geometry.

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
                Position of surface.
            refPoint : float
                Potential reference point.

            Returns
            -------
            CylindricalImagePotentail)V0G0N",
            "sim"_a, "ef"_a, "w"_a, "rad"_a, "posMin"_a, "posMax"_a, "surfPos"_a, "refPoint"_a)
        .def("negatePotential", &PyCylindricalImagePotential::negatePotential, R"V0G0N(
            Negates potential as evaluated in the Simulation's current state.
            Intended to be used such that the potential only depends on the change in density, leading the initially calculated eigenstates to be the actual eigenstates before perturbation.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.)V0G0N",
            "sim"_a);
}