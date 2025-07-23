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

    py::class_<UniformPotential, Potential>(m, "UniformPotential")
        .def(py::init([](SimulationManager* sim, double v0){
            return std::unique_ptr<UniformPotential>(new UniformPotential(
                sim->getNumPoints(), v0
            ));
        }), R"V0G0N(
            Static uniform potential.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            v0 : float
                Potential value.

            Returns
            -------
            UniformPotential)V0G0N",
            "sim"_a, "v0"_a);
    
    py::class_<OscillatingPotential, Potential>(m, "OscillatingPotential")
        .def(py::init([](SimulationManager* sim, Potential* basePot, double omega, double phase){
            return std::unique_ptr<OscillatingPotential>(new OscillatingPotential(
                sim->getNumPoints(), basePot, omega, phase
            ));
        }), py::keep_alive<1,3>(), R"V0G0N(
            Static potential with sinusoidal oscillation applied to it.
            The scalar product is $\sin(\omega t + \phi)$.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            basePot : Potential
                Base potential to apply the oscillation to.
            omega : float
                Oscillation frequency.
            phase : float
                Oscillation phase.

            Returns
            -------
            OscillatingPotential)V0G0N",
            "sim"_a, "basePot"_a, "omega"_a, "phase"_a);

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

    py::class_<CustomPotential, Potential>(m, "CustomPotential")
        .def(py::init([](SimulationManager* sim, const std::vector<double>& x, const std::vector<double>& v, double refPoint){
            return std::unique_ptr<CustomPotential>(new CustomPotential(
                sim->getNumPoints(), sim->getX(), x.size(), x.data(), v.data(), sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Static potential defined by a custom set of positions and values.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            x : list of float
                Positions at which the potential is defined.
            v : list of float
                Values of the potential at the specified positions.
            refPoint : float
                Potential reference point.

            Returns
            -------
            CustomPotential)V0G0N",
            "sim"_a, "x"_a, "v"_a, "refPoint"_a);

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
        .def(py::init([](SimulationManager* sim, int mimickOpenSystem, int ghostCharge, double rad, double posMin, double posMax, double surfPos, double refPoint){
            return std::unique_ptr<PlanarToCylindricalHartree>(new PlanarToCylindricalHartree(
                mimickOpenSystem, ghostCharge, sim->getNumPoints(), sim->getDX(), rad, sim->findXIdx(surfPos), sim->getNElecPtr(), sim->getWeightsPtr(), sim->wavefunctionIsInitialized() ? sim->getRho() : nullptr, sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(refPoint)
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
            mimickOpenSystem : int
                If nonzero, the charge is scaled to conserve the total charge, minus what leaves the selected boundary (-1 for left, 1 for right).
            ghostCharge : int
                If nonzero, charge lost on either boundary is replaced by a "ghost" of itself to preserve total charge even in truly open systems. -1 for left, 1 for right, 0 for no ghost charge.
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
            "sim"_a, "mimickOpenSystem"_a, "ghostCharge"_a, "rad"_a, "posMin"_a, "posMax"_a, "surfPos"_a, "refPoint"_a);
    
    py::class_<PlanarHartree, Potential>(m, "PlanarHartreePotential")
        .def(py::init([](SimulationManager* sim, double refPoint){
            return std::unique_ptr<PlanarHartree>(new PlanarHartree(
                sim->getNumPoints(), sim->getDX(), sim->wavefunctionIsInitialized() ? sim->getRho() : nullptr, sim->findXIdx(refPoint)
            ));
        }), R"V0G0N(
            Nonlocal Hartree potential assuming charge is distributed on a planar geometry.
            If the wavefunction has not been initialized upon construction then the potential will be returned as-is (with reference to refPoint).
            If the wavefunction has been initialized, then the potential will be returned as the change in potential with respect to the initial density.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            refPoint : float
                Potential reference point.

            Returns
            -------
            PlanarHartreePotential)V0G0N",
            "sim"_a, "refPoint"_a);
    
    py::class_<MixedGeometryHartree, Potential>(m, "MixedGeometryHartreePotential")
        .def(py::init([](SimulationManager* sim, double posMin, double posMax, double surfPos, double maskLength, double shieldLength, int neumannSide, double mRTheta, Densities::Density* dens, double refPoint, bool includeVectorPotential){
            
            // get geometry profile from passed density calculator
            std::unique_ptr<double*> hRad = std::make_unique<double*>(new double[sim->getNumPoints()]);
            std::fill_n(*hRad, sim->getNumPoints(), 1.0);
            dens->applyProfile(sim->getNumPoints(), sim->getNElec(), sim->getDX(), *hRad);
            for(size_t i = 0; i < sim->getNumPoints(); i++)
                (*hRad)[i] = 1.0 / (*hRad)[i]; // invert to get the geometry profile

            return std::unique_ptr<MixedGeometryHartree>(new MixedGeometryHartree(
                sim->getNumPoints(), sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(surfPos), maskLength, shieldLength, neumannSide, sim->getDX(), mRTheta, *hRad, 
                sim->wavefunctionIsInitialized() ? sim->getRho() : nullptr,
                (sim->wavefunctionIsInitialized() && includeVectorPotential) ? sim->getCur() : nullptr, 
                sim->findXIdx(refPoint), includeVectorPotential
            ));
        }), R"V0G0N(
            Hartree potential assuming charge is distributed on a mixed geometry defined by the passed Density calculator.
            Includes a shielded region where the change in charge is diminished in space, akin to, e.g., the Debye length or Thomas-Fermi screening length.
            If the wavefunction has not yet been initialized, the calculated potential will be returned as-is (with reference to refPoint).
            If the wavefunction has been initialized, then the potential will be returned as the change in potential with respect to the initial density.
            The transverse lengthscale may be included via mRTheta. e.g., for radius of curvature R and angular mode 1, mRTheta = 1.0/R.
            Includes a calculation of the vector potential followed by a gauge transformation to include it in the electrostatic potential. Applicable for high electron energies or quickly changing distributions.

            Parameters
            ----------
            sim : Simulation
                Associated simulation.
            posMin : float
                Minimum position where density is included.
            posMax : float
                Maximum position where density is included.
            surfPos : float
                Surface position.
            maskLength : float
                Lengthscale over which the source is masked. Good for truncating the system smoothly. If zero, Heaviside step functions are used.
            shieldLength : float
                Lengthscale over which the potential is shielded.
                Negative values indicate shielding to the left, positive values indicate shielding to the right.
            neumannSide : int
                Side of the Neumann boundary condition. -1 for left, 1 for right. The other boundary is Dirichlet.
            mRTheta : float
                Transverse lengthscale for the geometry profile. e.g., for radius of curvature R and angular mode 1, mRTheta = 1.0/R.
            dens : Density
                Density calculator to use for the geometry profile.
            refPoint : float
                Potential reference point.
            includeVectorPotential : bool
                If true, the vector potential is calculated and included in the electrostatic potential via a gauge transformation which preserves the electric field.
                Defaults to true.
            
            Returns
            -------
            MixedGeometryHartreePotential)V0G0N",
            "sim"_a, "posMin"_a, "posMax"_a, "surfPos"_a, "maskLength"_a, "shieldLength"_a, "neumannSide"_a, "mRTheta"_a, "dens"_a, "refPoint"_a, "includeVectorPotential"_a = true);

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
        .def(py::init([](SimulationManager* sim, Potential* pot, Measurers::Measurer* meas, size_t numSteps, bool measureVirtual){
            return std::unique_ptr<MeasuredPotential>(new MeasuredPotential(
                pot, meas, numSteps, numSteps*sim->getDT(), measureVirtual
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
            measureVirtual : bool
                If true, the potential is also measured during virtual steps (i.e., when getVVirtual is called).
                Defaults to false.

            Returns
            -------
            MeasuredPotential)V0G0N",
            "sim"_a, "pot"_a, "meas"_a, "numSteps"_a, "measureVirtual"_a = false);
            
    py::enum_<LDAFunctionalType>(m, "LDAFunctionalType")
        .value("X_SLATER", LDAFunctionalType::X_SLATER)
        .value("C_PW", LDAFunctionalType::C_PW)
        .export_values();
}