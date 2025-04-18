#pragma once
#include "MathTools.h"
#include "PyCommon.h"
#include "PySimulation.h"
#include "Potentials.h"
#include "Measurers.h"
#include "WfcRhoTools.h"

using namespace Potentials;

class PyFileFieldProfile
    : public ElectricFieldProfiles::FileFieldProfile{
        public:
        PyFileFieldProfile(PySimulation * sim, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, std::string fil)
            : ElectricFieldProfiles::FileFieldProfile(sim->getNumPoints(), sim->getX(), offset, rightDecayPos, leftDecayPos, decayLength, emax, fil){}
    };

// simplify instantiation of potentials
class PyFilePotential
    : public FilePotential{
        public:
        PyFilePotential(PySimulation * sim, double offset, std::string fil, double refPoint)
            : FilePotential(sim->getNumPoints(), sim->getX(), offset, fil, sim->findXIdx(refPoint)){}
    };

class PyJelliumPotential
    : public JelliumPotentialBacked{
        public:
        PyJelliumPotential(PySimulation * sim, double center, double ef, double w, double backStart, double backWidth, double refPoint)
            : JelliumPotentialBacked(sim->getNumPoints(), sim->getX(), center, ef, w, backStart, backWidth, sim->findXIdx(refPoint)){}
    };

class PyPulsePotential
    : public ElectricFieldProfileToPotential{
        public:
        PyPulsePotential(PySimulation* sim, ElectricFieldProfiles::ElectricFieldProfile* fieldProfile, Envelopes::Envelope * env, double phase, double tmax, double lam, double refPoint)
            : ElectricFieldProfileToPotential(sim->getNumPoints(), fieldProfile, sim->getDX(), phase, tmax, lam, env, sim->findXIdx(refPoint)){}
    };

class PyCylindricalImagePotential
    : public CylindricalImageCharge{
        public:
        PyCylindricalImagePotential(PySimulation* sim, double ef, double w, double rad, double posMin, double posMax, double surfPos, double refPoint)
            : CylindricalImageCharge(sim->getNumPoints(), sim->getX(), sim->getDX(), ef, w, rad, sim->findXIdx(surfPos), sim->getNElecPtr(), sim->getWeightsPtr(), sim->getRho(), sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(refPoint)){}
    };

class PyPlanarToCylindricalHartreePotential
    : public PlanarToCylindricalHartree{
        public:
        PyPlanarToCylindricalHartreePotential(PySimulation* sim, double rad, double posMin, double posMax, double surfPos, double refPoint)
            : PlanarToCylindricalHartree(sim->getNumPoints(), sim->getDX(), rad, sim->findXIdx(surfPos), sim->getNElecPtr(), sim->getWeightsPtr(), sim->getRho(), sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(refPoint)){}
    };

class PyLDAFunctional
    : public LDAFunctional{
        public:
        PyLDAFunctional(PySimulation* sim, LDAFunctionalType typ, double refPoint)
            : LDAFunctional(typ, sim->getNumPoints(), sim->getDX(), sim->getRho(), sim->findXIdx(refPoint)){}
    };

class PyMeasuredPotential
    : public MeasuredPotential{
        public:
        PyMeasuredPotential(PySimulation* sim, Potential* pot, Measurers::Measurer* meas, size_t numSteps)
            : MeasuredPotential(pot, meas, numSteps, numSteps*sim->getDT()){};
    };

void init_Potentials(py::module &m);