#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/BezierCrv.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename BezierCrvType>
void declare_bezier_crv(py::module & m, const std::string & name) {
    using BaseBezierCrvType = typename BezierCrvType::Base;

    // BaseBezierCrv
    std::string base_name = "Base" + name;
    declare_element<BezierCrvType>(m, base_name);
    py::class_<BaseBezierCrvType, Element<BezierCrvType>> (m, base_name.c_str());

    // BezierCrv
    py::class_<BezierCrvType, BaseBezierCrvType> c (m, name.c_str());
    c.def("__str__", [&name](const BezierCrvType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if constexpr(caribou::geometry::traits<BezierCrvType>::Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << name << " : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_bezier_crv(pybind11::module & m) {
    declare_bezier_crv<BezierCrv<_1D>>(m, "BezierCrv_1D");

    declare_bezier_crv<BezierCrv<_2D>>(m, "BezierCrv_2D");

    declare_bezier_crv<BezierCrv<_3D>>(m, "BezierCrv_3D");

    m.def("BezierCrv", [](){
        return py::cast(BezierCrv<_1D>());
    });


    m.def("BezierCrv", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_1D) {
            return py::cast(BezierCrv<_1D>());
        } else if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(BezierCrv<_2D>());
        } else {
            return py::cast(BezierCrv<_3D>());
        }
    }, py::arg("dimension"));


    // Linear creation
    // 1D
    m.def("BezierCrv", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & nodes) {
        return py::cast(BezierCrv<_1D>(nodes));
    }, py::arg("nodes"));

    m.def("BezierCrv", [](const FLOATING_POINT_TYPE & n0, const FLOATING_POINT_TYPE & n1, const FLOATING_POINT_TYPE & n2) {
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> nodes (n0, n1, n2);
        return py::cast(BezierCrv<_1D>(nodes));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 2D
    m.def("BezierCrv", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 2> & nodes) {
        return py::cast(BezierCrv<_2D>(nodes));
    }, py::arg("nodes"));

    m.def("BezierCrv", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> & n2) {
        return py::cast(BezierCrv<_2D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

    // 3D
    m.def("BezierCrv", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 3> & nodes) {
        return py::cast(BezierCrv<_3D>(nodes));
    }, py::arg("nodes"));

    m.def("BezierCrv", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n0, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n1, const Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> & n2) {
        return py::cast(BezierCrv<_3D>(n0, n1, n2));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"));

}

}
