#include "SimulationManager.h"
#include "KineticOperator.h"
#include "Potentials.h"
#include "Measurers.h"
#include "PhysCon.h"
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
    int
        Pointer to simulation instance.
    """
*/
    void* createSimulation( double xmin, double xmax, double dx, double dt, double maxT){
        int nPts = (int)((xmax-xmin)/dx);
        //Simulation is oblivious to particular coordinates, set here and pass as extra arg to store and retrieve later.
        double* xs = new double[nPts];
        for (int i = 0; i < nPts; i++)
            xs[i] = xmin + dx*i;

        return new(std::nothrow) SimulationManager(nPts, dx, dt, maxT, NULL, xs);
    }

/*
    """
    Deletes a TDSE simulation instance.

    Parameters
    ----------
    ptr : int
        Pointer to simulation instance.

    Returns
    -------
    int
        0 if successful.
    """
*/
    int deleteSimulation( void* ptr ){
        delete reinterpret_cast<SimulationManager*>(ptr);
        return 0;
    }

/*
    """
    Adds potential stored in file.

    Parameters
    ----------
    ptr : int
        Pointer to simulation instance.
    offset : float
        Positional offset.
    fil : str
        File containing potential.
    refPoint : int
        Maximum time in simulation (end time).

    Returns
    -------
    int
        Pointer to simulation instance.
    """
*/
    void addPot_FilePotential( void* ptr, double offset, string fil, int refPoint ){
        SimulationManager* sim = reinterpret_cast<SimulationManager*>(ptr);

        Potentials::Potential* nPot = new Potentials::FilePotential(
            sim->getNumPoints(), reinterpret_cast<double*>(sim->getargs()), offset, fil.c_str(), refPoint);

        sim->addPotential(nPot);
    }
}