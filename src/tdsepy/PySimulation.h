#pragma once
#include "KineticOperator.h"
#include "MathTools.h"
#include "PyCommon.h"
#include "pybind11/attr.h"
#include "pybind11/functional.h"
#include "SimulationManager.h"
#include "Densities.h"

class PySimulation 
    : public SimulationManager {
    private:
        char* wisdomFile = new char[50];
    public:
        PySimulation(double xmin, double xmax, double dx, double dt, const std::optional<std::function<void(double)>> &callback, std::optional<size_t> numCallbackCalls)
            : SimulationManager(xmin, xmax, dx, dt, callback.has_value() ? callback.value() : nullptr, numCallbackCalls.has_value() ? numCallbackCalls.value() : 101){
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

        void addLeftAbsBdy(double rate, double width){addSpatialDamp(vtls::getPolynomialSmoothBoundary(getNumPoints(), findXIdx(getX()[0]+width), 0, rate*getDT()).get());}
        void addRightAbsBdy(double rate, double width){addSpatialDamp(vtls::getPolynomialSmoothBoundary(getNumPoints(), findXIdx(getX()[getNumPoints()-1]-width), getNumPoints()-1, rate*getDT()).get());}

        void findEigenStates(double minE, double maxE){
            SimulationManager::findEigenStates(minE, maxE);
        }

        void findInhomogeneousEigenStates(size_t nElec, const py::array_t<double, py::array::c_style | py::array::forcecast> energies){
            SimulationManager::findInhomogeneousEigenStates(nElec, energies.data());
        }

        std::vector<double> getX(){return std::vector<double>(SimulationManager::getX(), SimulationManager::getX() + getNumPoints());}

        std::vector<double> getRho(){
            try{
                return std::vector<double>(SimulationManager::getRho(), SimulationManager::getRho() + getNumPoints());
            } catch (const std::runtime_error &e) {
                throw py::value_error(e.what());
            }
        }

        py::array_t<std::complex<double>> getPsi(){
            try{
                return py::array_t<std::complex<double>>(getNumPoints() * getNElec(), SimulationManager::getPsi());
            } catch (const std::runtime_error &e) {
                throw py::value_error(e.what());
            }
        }

        std::vector<double> getV(){
            try{
                return std::vector<double>(SimulationManager::getV(), SimulationManager::getV() + getNumPoints());
            } catch (const std::runtime_error &e) {
                throw py::value_error(e.what());
            }
        }

        std::vector<double> getWeightValues(){
            try{
                return std::vector<double>(SimulationManager::getWeightValues(), SimulationManager::getWeightValues() + getNElec());
            } catch (const std::runtime_error &e) {
                throw py::value_error(e.what());
            }
        }

        size_t findElectricalSurfaceCentroidRule(double minPos, double maxPos){
            return SimulationManager::findElectricalSurfaceCentroidRule(findXIdx(minPos), findXIdx(maxPos));
        }
    };

void init_Simulation(py::module &m);