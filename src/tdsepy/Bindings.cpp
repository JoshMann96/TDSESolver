#include "SimulationManager.h"
#include "KineticOperator.h"
#include "Potentials.h"
#include "Measurers.h"
#include "PhysCon.h"
#include <memory>
#include <new>

using namespace std;

extern "C" {

/*
    """
    Creates a TDSE Simulation object.

    Parameters
    ----------
    xmin : float
        Left boundary position.
    xmax : float
        Right boundary position.
    dx : float
        Spatial step size.
    dt : float
        Temporal step size.
    maxT : float
        Maximum time in simulation (end time).

    Returns
    -------
    ctypes.c_void_p
        Pointer to simulation instance.
    """
*/
    void* createSimulation( double xmin, double xmax, double dx, double dt, double maxT){
        int nPts = (int)((xmax-xmin)/dx);
        //Simulation is oblivious to particular coordinates, set here and pass as extra arg to store and retrieve later.
        double* xs = new double[nPts];
        for (int i = 0; i < nPts; i++)
            xs[i] = xmin + dx*i;

        SimulationManager* sim = new(std::nothrow) SimulationManager(nPts, dx, dt, maxT, NULL, xs);

        return reinterpret_cast<void*>(sim);
    }

/*
    """
    Deletes a TDSE simulation instance.

    Parameters
    ----------
    ptr : c_void_p
        Pointer to simulation instance.

    Returns
    -------
    int
        0 if successful.
    """
*/
    int deleteSimulation( void* ptr ){
        void* args = reinterpret_cast<SimulationManager*>(ptr)->getargs();
        delete[] reinterpret_cast<double*>(args);
        delete reinterpret_cast<SimulationManager*>(ptr);
        return 0;
    }

/*
    """
    Adds potential stored in file.

    Parameters
    ----------
    ptr : c_void_p
        Pointer to simulation instance.
    offset : float
        Positional offset.
    fil : str
        File containing potential.
    refPoint : int
        Maximum time in simulation (end time).
    """
*/
    void addPot_FilePotential( void* ptr, double offset, string fil, int refPoint ){
        SimulationManager* sim = reinterpret_cast<SimulationManager*>(ptr);

        Potentials::Potential* nPot = new Potentials::FilePotential(
            sim->getNumPoints(), reinterpret_cast<double*>(sim->getargs()), offset, fil.c_str(), refPoint);

        sim->addPotential(nPot);
    }

/*
    """
    Adds Jellium potential.

    Parameters
    ----------
    ptr : c_void_p
        Pointer to simulation instance.
    center : float
        Jellium center (sigmoid center).
    ef : float
        Fermi energy.
    w : float
        Work function.
    backStart : float
        Internal start position of polynomial smooth spline backing.
    backWidth : float
        Width of polynomial smooth spline backing.
    refPoint : int
        Maximum time in simulation (end time).
    """
*/
    void addPot_JelliumPotentialBacked( void* ptr, double center, double ef, double w, double backStart, double backWidth, int refPoint ){
        SimulationManager* sim = reinterpret_cast<SimulationManager*>(ptr);

        Potentials::Potential* nPot = new Potentials::JelliumPotentialBacked(
            sim->getNumPoints(), reinterpret_cast<double*>(sim->getargs()), center, ef, w, backStart, backWidth, refPoint);

        std::cout << "0" << std::endl;

        sim->addPotential(nPot);

        std::cout << "1" << std::endl;
    }

    int testPtr( void* ptr ){
        SimulationManager* sim = reinterpret_cast<SimulationManager*>(ptr);

        return sim->getNumPoints();
    }

}