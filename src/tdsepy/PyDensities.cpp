#include "PyDensities.h"
#include "Densities.h"
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pytypes.h>

void init_Densities(py::module &m) {
    py::class_<Densities::Weight>(m, "Weight");

    py::class_<Densities::UniformWeight, Densities::Weight>(m, "UniformWeight")
        .def(py::init<double>(), R"V0G0N(
            Sets the weights to be a constant value for all states.

            Parameters
            ----------
            weight : float
                Weight to be assigned to each state.

            Returns
            -------
            UniformWeight)V0G0N",
            "weight"_a);

    py::class_<Densities::BoundFermiGas, Densities::Weight>(m, "BoundFermiGas")
        .def(py::init<double>(), R"V0G0N(
            Uses 3-D Fermi gas distribution at zero temperature to convert 1-D wavefunctions to an effective 3-D density.
            Mapping uses wavefunction initial eigenstates assuming they are bound and normalized.

            Parameters
            ----------
            ef : float
                Fermi energy.

            Returns
            -------
            BoundFermiGas)V0G0N",
            "ef"_a);
    
    py::class_<Densities::SemiInfiniteFermiGas, Densities::Weight>(m, "SemiInfiniteFermiGas")
        .def(py::init<double>(), R"V0G0N(
            Uses 3-D Fermi gas distribution at zero temperature to convert 1-D wavefunctions to an effective 3-D density.
            Mapping uses wavefunction initial eigenstates assuming they are eigenstates of the open system, in contact with a zero-temperature free electron gas (FEG).
            The boundary conditions on the FEG side must be inhomogeneous and have an incoming component of magnitude 1.

            Parameters
            ----------
            ef : float
                Fermi energy.

            Returns
            -------
            SemiInfiniteFermiGas)V0G0N",
            "ef"_a);

// DENSITY PROCESSING

    py::class_<Densities::Density>(m, "Density");

    py::class_<Densities::DirectDensity, Densities::Density>(m, "DirectDensity")
        .def(py::init<>(), R"V0G0N(
            Uses no preprocessing in calculating the final density.
            Density = sum over states (weight x psi*psi)

            Returns
            -------
            DirectDensity)V0G0N");

    py::class_<Densities::GaussianSmoothedDensity, Densities::Density>(m, "GaussianSmoothedDensity")
        .def(py::init<double>(), R"V0G0N(
            Uses a Gaussian smoothing function to calculate the final density.
            Density = sum over states (weight x psi*psi) * Gaussian
                (* = convolution)

            Parameters
            ----------
            sigma : float
                Standard deviation of Gaussian. Typically the inverse of the Thomas-Fermi wavenumber.

            Returns
            -------
            GaussianSmoothedDensity)V0G0N",
            "sigma"_a)
        .def(py::init<double, Densities::Density*>(), py::keep_alive<1,3>(), R"V0G0N(
            Uses a Gaussian smoothing function to calculate the final density.
            Density = sum over states (weight x psi*psi) * Gaussian
                (* = convolution)

            Parameters
            ----------
            sigma : float
                Standard deviation of Gaussian. Typically the inverse of the Thomas-Fermi wavenumber.
            baseDens : Density
                Base density calculator to use for the raw density. This is applied first, then the smoothing is applied.

            Returns
            -------
            GaussianSmoothedDensity)V0G0N",
            "sigma"_a, "baseDens"_a);
    
    py::class_<Densities::CylindricalDensity, Densities::Density>(m, "CylindricalDensity")
        .def(py::init<double, double, double>(), R"V0G0N(
            Applies geometric dispersion to density assuming a cylindrical geometry (azimuthal symmetry) with some definite radius.


            Parameters
            ----------
            center : float
                Center of azimuthal symmetry.
            radius : float 
                Radius of azimuthal symmetry.
            minX : float
                Left boundary position.

            Returns
            -------
            CylindricalDensity)V0G0N",
            "center"_a, "radius"_a, "minX"_a);
}