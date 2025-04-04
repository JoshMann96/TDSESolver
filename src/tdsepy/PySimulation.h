#pragma once
#include "KineticOperator.h"
#include "MathTools.h"
#include "PyCommon.h"
#include "pybind11/attr.h"
#include "pybind11/functional.h"
#include "SimulationManager.h"
#include "WfcRhoTools.h"
#include "AbsorptiveRegions.h"

class PySimulation 
    : public SimulationManager {
    private:
        char* wisdomFile = new char[50];
    public:
        PySimulation(double xmin, double xmax, double dx, double dt, const std::function<void(int)> &callback)
            : SimulationManager(xmin, xmax, dx, dt, callback){
			    std::snprintf(wisdomFile, 50, "fftw_nt_%04d.wisdom", omp_get_max_threads());
                fftw_init_threads();
				fftw_import_wisdom_from_filename(wisdomFile);
        }

        ~PySimulation(){
            fftw_export_wisdom_to_filename(wisdomFile);
            fftw_cleanup_threads();
            delete[] wisdomFile;
        }

        void addPotential(Potentials::Potential * pot){
            SimulationManager::addPotential(pot);
        }

        void addLeftAbsBdy(double rate, double width){addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(getNumPoints(), findXIdx(getX()[0]+width), 0, rate*getDT()).get());}
        void addRightAbsBdy(double rate, double width){addSpatialDamp(AbsorptiveRegions::getSmoothedSpatialDampDecay(getNumPoints(), findXIdx(getX()[getNumPoints()-1]-width), getNumPoints()-1, rate*getDT()).get());}

        void findEigenStates(double minE, double maxE){
            SimulationManager::findEigenStates(minE, maxE);
        }

        void findInhomogeneousEigenStates(int nElec, const std::vector<double> &energies){
            SimulationManager::findInhomogeneousEigenStates(nElec, energies.data());
        }

        std::vector<double> getXVec(){return std::vector<double>(SimulationManager::getX(), SimulationManager::getX() + getNumPoints());}

        int findElectricalSurfaceCentroidRule(double minPos, double maxPos){
            return SimulationManager::findElectricalSurfaceCentroidRule(findXIdx(minPos), findXIdx(maxPos));
        }
    };

void init_Simulation(py::module &m);