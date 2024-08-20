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
            : ElectricFieldProfiles::FileFieldProfile(sim->getNumPoints(), sim->getXPtr(), offset, rightDecayPos, leftDecayPos, decayLength, emax, fil.c_str()){}
    };

// simplify instantiation of potentials
class PyFilePotential
    : public FilePotential{
        public:
        PyFilePotential(PySimulation * sim, double offset, std::string fil, double refPoint)
            : FilePotential(sim->getNumPoints(), sim->getXPtr(), offset, fil.c_str(), sim->findXIdx(refPoint)){}
    };

class PyJelliumPotential
    : public JelliumPotentialBacked{
        public:
        PyJelliumPotential(PySimulation * sim, double center, double ef, double w, double backStart, double backWidth, double refPoint)
            : JelliumPotentialBacked(sim->getNumPoints(), sim->getXPtr(), center, ef, w, backStart, backWidth, sim->findXIdx(refPoint)){}
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
            : CylindricalImageCharge(sim->getNumPoints(), sim->getXPtr(), sim->getDX(), ef, w, rad, sim->getNElecPtr(), sim->getWeightsPtr(), sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(surfPos), sim->findXIdx(refPoint)){}
        
        void negatePotential(PySimulation* sim){
            CylindricalImageCharge::negateGroundEffects(sim->getRho(), sim->getPsi());
        }
    };

class PyPlanarToCylindricalHartreePotential
    : public PlanarToCylindricalHartree{
        public:
        PyPlanarToCylindricalHartreePotential(PySimulation* sim, double rad, double posMin, double posMax, double surfPos, double refPoint)
            : PlanarToCylindricalHartree(sim->getNumPoints(), sim->getXPtr(), sim->getDX(), rad, sim->getNElecPtr(), sim->getWeightsPtr(), sim->findXIdx(posMin), sim->findXIdx(posMax), sim->findXIdx(surfPos), sim->findXIdx(refPoint)){}
        
        void negatePotential(PySimulation* sim){
            PlanarToCylindricalHartree::negateGroundEffects(sim->getRho(), sim->getPsi());
        }
    };

class PyLDAFunctional
    : public LDAFunctional{
        public:
        PyLDAFunctional(PySimulation* sim, LDAFunctionalType typ, double refPoint)
            : LDAFunctional(typ, sim->getNumPoints(), sim->getDX(), sim->findXIdx(refPoint)){}
        
        void negatePotential(PySimulation* sim){
            LDAFunctional::negateGroundEffects(sim->getRho(), sim->getPsi());
        }
    };

class PyMeasuredPotential
    : public MeasuredPotential{
        public:
        PyMeasuredPotential(PySimulation* sim, Potential* pot, Measurers::Measurer* meas)
            : MeasuredPotential(pot, meas, sim->getNumSteps(), sim->getMaxT()){}
    };

void init_Potentials(py::module &m);