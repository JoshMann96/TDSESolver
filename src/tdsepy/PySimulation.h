#pragma once
#include "PyCommon.h"
#include "SimulationManager.h"

// couple particular grid with SimulationManager for ease of use
class PySimulation 
    : public SimulationManager {
    private:
        double* x;
        int nPts;
    public:
        PySimulation(double xmin, double xmax, double dx, double dt, double maxT)
            : nPts((int)((xmax-xmin)/dx)), 
              SimulationManager(nPts, dx, dt, maxT, NULL){
                x = new double[nPts];
                for(int i = 0; i < nPts; i++)
                    x[i] = i*dx + xmin;
        }
        ~PySimulation(){SimulationManager::~SimulationManager(); delete[] x;}

        void addPotential(Potentials::Potential * pot){
            SimulationManager::addPotential(pot);
        }

        double* getXPtr(){return x;}
        std::vector<double> getX(){return std::vector<double>(x, x + nPts);}
        int findXIdx(double xp){return vtls::findValue(nPts, x, xp);}
    };

void init_Simulation(py::module &m);