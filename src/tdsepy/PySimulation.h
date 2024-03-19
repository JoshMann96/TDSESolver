#pragma once
#include "MathTools.h"
#include "PyCommon.h"
#include "SimulationManager.h"
#include "WfcRhoTools.h"

// couple particular grid with SimulationManager for ease of use
class PySimulation 
    : public SimulationManager {
    private:
        double* x;
        int nPts;
        WfcToRho::Weight* wght = nullptr;
        WfcToRho::Density* dens = nullptr;
    public:
        PySimulation(double xmin, double xmax, double dx, double dt, double maxT)
            : SimulationManager((int)((xmax-xmin)/dx) + 1, dx, dt, maxT, NULL){
                nPts = getNumPoints();
                x = new double[nPts];
                for(int i = 0; i < nPts; i++)
                    x[i] = i*dx + xmin;
        }
        ~PySimulation(){SimulationManager::~SimulationManager(); delete[] x;}

        void addPotential(Potentials::Potential * pot){
            SimulationManager::addPotential(pot);
        }

        double* getXPtr(){return x;}

        WfcToRho::Weight* getWght(){return wght;}
        void setWght(WfcToRho::Weight* wght){this->wght = wght;}
        WfcToRho::Density* getDens(){return dens;}
        void setDens(WfcToRho::Density* dens){this->dens = dens;}

        std::vector<double> getX(){return std::vector<double>(x, x + nPts);}
        int findXIdx(double xp){return vtls::findValue(nPts, x, xp);}
    };

void init_Simulation(py::module &m);