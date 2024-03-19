#include "PyMeasurers.h"

using namespace Measurers;

void init_Measurers(py::module &m) {
    py::class_<Measurer>(m, "Measurer");

    py::class_<PyConstant Measurer>(m, "Constant")
        .def(py::init<double, std::string, std::string>, R"V0G0N(
            Records a constant.

            Parameters
            ----------
            c : float
                Value to be recorded.
            fileName : str
                Name of output file (no extension).
            fol : str
                Directory to contain file.

            Returns
            -------
            Constant)V0G0N",
            "ef"_a);
}