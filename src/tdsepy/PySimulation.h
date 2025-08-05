#pragma once
#include "KineticOperator.h"
#include "MathTools.h"
#include "PyCommon.h"
#include "pybind11/attr.h"
#include "pybind11/functional.h"
#include "SimulationManager.h"
#include "Densities.h"
#include <boost/asio/ip/host_name.hpp>

class PySimulation 
    : public SimulationManager {
    private:
        char* wisdomFile = nullptr;
        void initFFTW(const std::string &fftwWisdomPrefix) {
            const size_t maxLen = 64 + boost::asio::ip::host_name().length() + fftwWisdomPrefix.length();

            if (wisdomFile){
                delete[] wisdomFile;
                wisdomFile = nullptr;
            }
            wisdomFile = new char[maxLen];

            std::snprintf(wisdomFile, maxLen, "%sfftw_nt_%04d_%s.wisdom", 
                fftwWisdomPrefix.c_str(), 
                omp_get_max_threads(), 
                boost::asio::ip::host_name().c_str());

            fftw_init_threads();
			fftw_import_wisdom_from_filename(wisdomFile);
        }
    public:
        PySimulation(double xmin, double xmax, double dx, double dt, const std::optional<std::function<void(double)>> &callback, std::optional<size_t> numCallbackCalls, std::optional<std::string> fftwWisdomPrefix)
            : SimulationManager(xmin, xmax, dx, dt, callback.has_value() ? callback.value() : nullptr, numCallbackCalls.has_value() ? numCallbackCalls.value() : 101)
            {initFFTW(fftwWisdomPrefix.has_value() ? fftwWisdomPrefix.value() : "");}

        PySimulation(size_t nPts, double xmin, double dx, double dt, const std::optional<std::function<void(double)>> &callback, std::optional<size_t> numCallbackCalls, std::optional<std::string> fftwWisdomPrefix)
            : SimulationManager(nPts, xmin, dx, dt, callback.has_value() ? callback.value() : nullptr, numCallbackCalls.has_value() ? numCallbackCalls.value() : 101)
            {initFFTW(fftwWisdomPrefix.has_value() ? fftwWisdomPrefix.value()  : "");}
        
        PySimulation(double xmin, double xmax, size_t nPts, double dt, const std::optional<std::function<void(double)>> &callback, std::optional<size_t> numCallbackCalls, std::optional<std::string> fftwWisdomPrefix)
            : SimulationManager(xmin, xmax, nPts, dt, callback.has_value() ? callback.value() : nullptr, numCallbackCalls.has_value() ? numCallbackCalls.value() : 101)
            {initFFTW(fftwWisdomPrefix.has_value() ? fftwWisdomPrefix.value()  : "");}

        ~PySimulation(){
            fftw_export_wisdom_to_filename(wisdomFile);
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

        std::vector<double> getCur(){
            try{
                return std::vector<double>(SimulationManager::getCur(), SimulationManager::getCur() + getNumPoints());
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

        void setPsi(const py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> &psi, size_t nElec = 1, Densities::NormalizationScheme norm = Densities::UNNORMALIZED) {
            if (psi.size() != getNumPoints() * nElec) {
                throw py::value_error("Wavefunction size does not match the number of grid points and states.");
            }
            SimulationManager::setPsi(psi.data(), nElec, norm);
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