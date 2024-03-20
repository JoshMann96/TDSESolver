#pragma once
#include "PyCommon.h"
#include "Measurers.h"
#include "PySimulation.h"

void init_Measurers(py::module &m);

using namespace Measurers;

class PyConstant
    : public DoubleConst {
        public:
        PyConstant(double c, std::string filName, std::string fol)
            : DoubleConst(c, filName.c_str(), fol.c_str()){}
    };

class PyBasic
    : public BasicMeasurers{
        public:
        PyBasic(PySimulation* sim, std::string fol)
            : BasicMeasurers(sim->getNumPoints(), sim->getDX(), sim->getDT(), sim->getXPtr(), fol.c_str()){}
    };

class PyXS
    : public XS{
        public:
        PyXS(PySimulation* sim, std::string fol)
            : XS(sim->getNumPoints(), sim->getXPtr(), fol.c_str()){}
    };

class PyTS
    : public TS{
        public:
        PyTS(std::string fol)
            : TS(fol.c_str()){}
    };

class PyOrigPot
    : public OrigPot{
        public:
        PyOrigPot(PySimulation* sim, std::string fol)
            : OrigPot(sim->getNumPoints(), fol.c_str()){}
    };

class PyPsi2t
    : public Psi2t{
        public:
        PyPsi2t(PySimulation* sim, int nx, int nt, std::string fol)
            : Psi2t(sim->getNumPoints(), nx, nt, sim->getMaxT(), sim->getDT(), sim->getXPtr(), fol.c_str()){}
    };

class PyVfunct
    : public Vfunct{
        public:
        PyVfunct(PySimulation* sim, int nx, int nt, std::string fol)
            : Vfunct(sim->getNumPoints(), nx, nt, sim->getMaxT(), sim->getXPtr(), fol.c_str()){}
    };

/*
____________________
ADD WARNING
USES 3-POINT STENCIL DERIVATIVE FOR KINETIC OPERATOR
SHOULD BE UPDATED IF THIS BECOMES IMPORTANT
APPLIES FOR
PyExpectE, PyExpectP
____________________
*/
class PyExpectE
    : public ExpectE{
        public:
        PyExpectE(PySimulation* sim, std::string fol)
            : ExpectE(sim->getNumPoints(), sim->getDX(), fol.c_str()){}
    };

class PyExpectE0
    : public ExpectE0{
        public:
        PyExpectE0(PySimulation* sim, std::string fol)
            : ExpectE0(sim->getNumPoints(), sim->getDX(), fol.c_str()){}
    };

class PyExpectX
    : public ExpectX{
        public:
        PyExpectX(PySimulation* sim, std::string fol)
            : ExpectX(sim->getNumPoints(), sim->getXPtr(), sim->getDX(), fol.c_str()){}
    };

//Computationally expensive
class PyExpectP
    : public ExpectP{
        public:
        PyExpectP(PySimulation* sim, std::string fol)
            : ExpectP(sim->getNumPoints(), sim->getDX(), fol.c_str()){}
    };

class PyExpectA
    : public ExpectA{
        public:
        PyExpectA(PySimulation* sim, std::string fol)
            : ExpectA(sim->getNumPoints(), sim->getDX(), fol.c_str()){}
    };

class PyTotProb
    : public TotProb{
        public:
        PyTotProb(PySimulation* sim, std::string fol)
            : TotProb(sim->getNumPoints(), sim->getDX(), fol.c_str()){}
    };

// name is 4 characters
class PyVDProbCurrent
    : public VDProbCurrent{
        public:
        PyVDProbCurrent(PySimulation* sim, double vdPos, int vdNum, std::string name, std::string fol)
            : VDProbCurrent(sim->getNumPoints(), sim->getDX(), sim->findXIdx(vdPos), vdNum, name.c_str(),fol.c_str()){}
    };

class PyVDPsi
    : public VDPsi{
        public:
        PyVDPsi(PySimulation* sim, double vdPos, int vdNum, std::string name, std::string fol)
            : VDPsi(sim->findXIdx(vdPos), vdNum, name.c_str(),fol.c_str()){}
    };

class PyVDPot
    : public VDPot{
        public:
        PyVDPot(PySimulation* sim, double vdPos, int vdNum, std::string name, std::string fol)
            : VDPot(sim->findXIdx(vdPos), vdNum, name.c_str(),fol.c_str()){}
    };

class PyVDFluxSpec
    : public VDFluxSpec{
        public:
        PyVDFluxSpec(PySimulation* sim, double vdPos, int vdNum, int nSamp, double emax, std::string name, std::string fol)
            : VDFluxSpec(sim->findXIdx(vdPos), vdNum, sim->getNElecPtr(), nSamp, emax, sim->getMaxT(), name.c_str(),fol.c_str()){}
    };

class PyPsiT
    : public PsiT{
        public:
        PyPsiT(PySimulation* sim, double meaT, int vdNum, std::string name, std::string fol)
            : PsiT(sim->getNumPoints(), meaT, vdNum, name.c_str(),fol.c_str()){}
    };

class PyPotT
    : public PotT{
        public:
        PyPotT(PySimulation* sim, double meaT, int vdNum, std::string name, std::string fol)
            : PotT(sim->getNumPoints(), meaT, vdNum, name.c_str(),fol.c_str()){}
    };

class PyNElec
    : public ElectronNumber{
        public:
        PyNElec(PySimulation* sim, std::string fol)
            : ElectronNumber(sim->getNElecPtr(), fol.c_str()){}
    };

class PyWeights
    : public WfcRhoWeights{
        public:
        PyWeights(PySimulation* sim, std::string fol)
            : WfcRhoWeights(sim->getNElecPtr(), sim->getNumPoints(), sim->getDX(), sim->getWght(), fol.c_str()){}
    };