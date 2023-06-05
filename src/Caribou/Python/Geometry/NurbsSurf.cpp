#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/NurbsSurf.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename NurbsSurfType>
void declare_nurbs_surf(py::module & m, const std::string & name) {
    using BaseNurbsSurfType = typename NurbsSurfType::Base;

    // BaseNurbsSurf
    std::string base_name = "Base" + name;
    declare_element<NurbsSurfType>(m, base_name);
    py::class_<BaseNurbsSurfType, Element<NurbsSurfType>> (m, base_name.c_str());

    // NurbsSurf
    py::class_<NurbsSurfType, BaseNurbsSurfType> c (m, name.c_str());
    c.def("__str__", [&name](const NurbsSurfType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if constexpr(caribou::geometry::traits<NurbsSurfType>::Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << name << " : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_nurbs_surf(pybind11::module & m) {

    declare_nurbs_surf<NurbsSurf<_2D>>(m, "NurbsSurf_2D");

    declare_nurbs_surf<NurbsSurf<_3D>>(m, "NurbsSurf_3D");

    m.def("NurbsSurf", [](){
        return py::cast(NurbsSurf<_2D>());
    });


    m.def("NurbsSurf", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(NurbsSurf<_2D>());
        } else {
            return py::cast(NurbsSurf<_3D>());
        }
    }, py::arg("dimension"));

    using dyn_vec = NurbsSurf<_2D>::Dyn_Vector;

    // 2D
    m.def("NurbsSurf", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 9, 2> & nodes,dyn_vec knot_1,
          const dyn_vec knot_2, const dyn_vec weights, const dyn_vec knot_span)
    {
        return py::cast(NurbsSurf<_2D>(nodes, knot_1, knot_2, weights, knot_span));
    }, py::arg("nodes"), py::arg("knot_1"), py::arg("knot_2"), py::arg("weights"), py::arg("knot_span"));

    m.def("jacobian_papa", &NurbsSurf<_2D>::jacobian_papa);

    // 3D

    m.def("NurbsSurf", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 9, 3> & nodes, dyn_vec knot_1,
          const dyn_vec knot_2, const dyn_vec weights, const dyn_vec knot_span)
    {
        return py::cast(NurbsSurf<_3D>(nodes, knot_1, knot_2, weights, knot_span));
    }, py::arg("nodes"), py::arg("knot_1"), py::arg("knot_2"), py::arg("weights"), py::arg("knot_span"));

    m.def("jacobian_papa", &NurbsSurf<_3D>::jacobian_papa);
}

}
