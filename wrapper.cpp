#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
//#include <pybind11/numpy.h>
#include "GreenAnisotropic2D.hpp"

namespace py = pybind11;

using namespace pybind11::literals;

PYBIND11_PLUGIN(GreenAnisotropic2D) {
    py::module m("GreenAnisotropic2D", "A module to compute the Green's function of 2D anisotropic elastic media and its gradients.");
    py::class_<medium>(m, "medium")
        .def(py::init<>())
        .def_readonly("isym",&medium::isym)
        .def_readonly("P0",&medium::P0)
        .def_readonly("P1",&medium::P1)
        .def_readonly("R0",&medium::R0)
        .def_readonly("R1",&medium::R1)
        .def_readonly("T0",&medium::T0)
        .def_readonly("T1",&medium::T1)
        .def_readonly("K",&medium::K)
        .def_readonly("Th",&medium::Th)
        .def_readonly("S",&medium::S)
        .def_readonly("H",&medium::H)
        .def_readonly("L1111",&medium::L1111)
        .def_readonly("L1122",&medium::L1122)
        .def_readonly("L1112",&medium::L1112)
        .def_readonly("L2211",&medium::L2211)
        .def_readonly("L2222",&medium::L2222)
        .def_readonly("L2212",&medium::L2212)
        .def_readonly("L1211",&medium::L1211)
        .def_readonly("L1222",&medium::L1222)
        .def_readonly("L1212",&medium::L1212)
        .def_readonly("L1121",&medium::L1121)
        .def_readonly("L2221",&medium::L2221)
        .def_readonly("L2111",&medium::L2111)
        .def_readonly("L2122",&medium::L2122)
        .def_readonly("L2112",&medium::L2112)
        .def_readonly("L1221",&medium::L1221)
        .def_readonly("L2121",&medium::L2121)
        .def("set_sym",&medium::set_sym)
		.def("get_S",&medium::get_S)
		.def("get_H",&medium::get_H)
		.def("dnGi",&medium::dnGi);
		// Access to the following functions is not granted through
		// the Python module. Otherwise, uncomment. 
		//.def("dkN1_anisotropic", &medium::dkN1_anisotropic)
		//.def("dkN2_anisotropic", &medium::dkN2_anisotropic)
		//.def("dkN1_orthotropic", &medium::dkN1_orthotropic)
		//.def("dkN2_orthotropic", &medium::dkN2_orthotropic)
		//.def("dkN1_r0_orthotropic", &medium::dkN1_r0_orthotropic)
		//.def("dkN2_r0_orthotropic", &medium::dkN2_r0_orthotropic)
		//.def("dkN1_square_symmetric", &medium::dkN1_square_symmetric)
		//.def("dkN2_square_symmetric", &medium::dkN2_square_symmetric)
		//.def("dkN1_isotropic", &medium::dkN1_isotropic)
		//.def("dkN2_isotropic", &medium::dkN2_isotropic)
		//.def("dh",&medium::dh)
		//.def("h",&medium::h)
		//.def("dkh1i",&medium::dkh1i)
		//.def("dknvec",&medium::dknvec)
		//.def("dkmvec",&medium::dkmvec);	
    //m.def("Binom", &Binom, "Returns a binomial coefficient",py::arg("n"),py::arg("k"));	
return m.ptr();
}



