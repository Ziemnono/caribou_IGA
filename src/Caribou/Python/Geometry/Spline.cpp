#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Spline.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename SplineType>
void declare_spline(py::module & m, const std::string & name) {
    using BaseSplineType = typename SplineType::Base;

    // BaseSpline
    std::string base_name = "Base" + name;
    declare_element<SplineType>(m, base_name);
    py::class_<BaseSplineType, Element<SplineType>> (m, base_name.c_str());

    // Spline
    py::class_<SplineType, BaseSplineType> c (m, name.c_str());
    c.def("__str__", [&name](const SplineType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if constexpr(caribou::geometry::traits<SplineType>::Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << name << " : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_spline(pybind11::module & m) {
    declare_spline<Spline<_1D>>(m, "Spline_1D");

    declare_spline<Spline<_2D>>(m, "Spline_2D");

    declare_spline<Spline<_3D>>(m, "Spline_3D");

    m.def("Spline", [](){
        return py::cast(Spline<_1D>());
    });


    m.def("Spline", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_1D) {
            return py::cast(Spline<_1D>());
        } else if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(Spline<_2D>());
        } else {
            return py::cast(Spline<_3D>());
        }
    }, py::arg("dimension"));


    // Linear creation
    // 1D
    m.def("Spline", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & nodes) {
        return py::cast(Spline<_1D>(nodes));
    }, py::arg("nodes"));

    m.def("Spline", [](const FLOATING_POINT_TYPE & n0, const FLOATING_POINT_TYPE & n1, const FLOATING_POINT_TYPE & n2) {
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> nodes (n0, n1, n2);
        return py::cast(Spline<_1D>(nodes));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 2D
    m.def("Spline", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 2> & nodes) {
        return py::cast(Spline<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("Spline", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2) {
        return py::cast(Spline<_2D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 3D
    m.def("Spline", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> & nodes) {
        return py::cast(Spline<_3D>(nodes));
    }, py::arg("nodes"));

    m.def("Spline", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2) {
        return py::cast(Spline<_3D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));


}

}
