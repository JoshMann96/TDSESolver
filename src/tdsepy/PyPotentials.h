#pragma once
#include "MathTools.h"
#include "PyCommon.h"
#include "PySimulation.h"
#include "Potentials.h"
#include "WfcRhoTools.h"

class PyFileFieldProfile
    : public ElectricFieldProfiles::FileFieldProfile{
        public:
        PyFileFieldProfile(PySimulation * sim, double offset, double rightDecayPos, double leftDecayPos, double decayLength, double emax, std::string fil)
            : ElectricFieldProfiles::FileFieldProfile(sim->getNumPoints(), sim->getXPtr(), offset, rightDecayPos, leftDecayPos, decayLength, emax, fil.c_str()){}
    };

// simplify instantiation of potentials
class PyFilePotential
    : public Potentials::FilePotential{
        public:
        PyFilePotential(PySimulation * sim, double offset, std::string fil, double refPoint)
            : Potentials::FilePotential(sim->getNumPoints(), sim->getXPtr(), offset, fil.c_str(), sim->findXIdx(refPoint)){}
    };

class PyJelliumPotential
    : public Potentials::JelliumPotentialBacked{
        public:
        PyJelliumPotential(PySimulation * sim, double center, double ef, double w, double backStart, double backWidth, double refPoint)
            : Potentials::JelliumPotentialBacked(sim->getNumPoints(), sim->getXPtr(), center, ef, w, backStart, backWidth, sim->findXIdx(refPoint)){}
    };

class PyPulsePotential
    : public Potentials::ElectricFieldProfileToPotential{
        public:
        PyPulsePotential(PySimulation* sim, ElectricFieldProfiles::ElectricFieldProfile* fieldProfile, Envelopes::Envelope * env, double phase, double tmax, double lam, double refPoint)
            : Potentials::ElectricFieldProfileToPotential(sim->getNumPoints(), fieldProfile, sim->getDX(), phase, tmax, lam, env, sim->findXIdx(refPoint)){}
    };

class PyCylindricalImagePotential
    : public Potentials::CylindricalImageCharge{
        public:
        PyCylindricalImagePotential(PySimulation* sim, double ef, double w, double rad, double posMin, double posMax, double surfPos, double refPoint)
            : Potentials::CylindricalImageCharge(sim->getNumPoints(), sim->getXPtr(), sim->getDX(), sim->getDT(), ef, w, rad, sim->getNElecPtr(), sim->getPotPointer(), sim->getWght(), sim->getDens(), sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(surfPos), sim->findXIdx(refPoint)){}
    };

void init_Potentials(py::module &m);