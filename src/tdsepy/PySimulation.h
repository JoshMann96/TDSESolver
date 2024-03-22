#pragma once
#include "KineticOperator.h"
#include "MathTools.h"
#include "PyCommon.h"
#include "pybind11/attr.h"
#include "pybind11/functional.h"
#include "SimulationManager.h"
#include "WfcRhoTools.h"
#include "AbsorptiveRegions.h"

// couple particular grid with SimulationManager for ease of use
class PySimulation 
    : public SimulationManager {
    private:
        double* x;
        int nPts;
        WfcToRho::Weight* wght = nullptr;
        WfcToRho::Density* dens = nullptr;
        char* wisdomFile = new char[50];
    public:
        PySimulation(double xmin, double xmax, double dx, double dt, double maxT, const std::function<void(int)> &callback)
            : SimulationManager((int)((xmax-xmin)/dx) + 1, dx, dt, maxT, callback){
                nPts = getNumPoints();
                x = new double[nPts];
                for(int i = 0; i < nPts; i++)
                    x[i] = i*dx + xmin;
                
			    std::snprintf(wisdomFile, 50, "fftw_nt_%04d.wisdom", omp_get_max_threads());
                fftw_init_threads();
				fftw_import_wisdom_from_filename(wisdomFile);
        }

        ~PySimulation(){fftw_export_wisdom_to_filename(wisdomFile); SimulationManager::~SimulationManager(); delete[] x;}

        void addPotential(Potentials::Potential * pot){
            SimulationManager::addPotential(pot);
        }

        double* getXPtr(){return x;}

        WfcToRho::Weight* getWght(){return wght;} //LOOK HERE: keep_alive
        void setWght(WfcToRho::Weight* wght){this->wght = wght;}
        WfcToRho::Density* getDens(){return dens;}
        void setDens(WfcToRho::Density* dens){this->dens = dens;}

        void addLeftAbsBdy(double rate, double width){addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(nPts, findXIdx(x[0]+width), 0, rate*getDT()));}
        void addRightAbsBdy(double rate, double width){addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(nPts, findXIdx(x[nPts-1]-width), nPts-1, rate*getDT()));}

        void findEigenStates(double minE, double maxE){
            SimulationManager::findEigenStates(minE, maxE, 0.0, 0.0);
        }

        std::vector<double> getX(){return std::vector<double>(x, x + nPts);}
        int findXIdx(double xp){return vtls::findValue(nPts, x, xp);}
    };

void init_Simulation(py::module &m);