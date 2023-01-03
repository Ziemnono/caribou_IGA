#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Eigen/Core>
#include <Caribou/constants.h>
#include <Caribou/Geometry/BezierSurf.h>
#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Geometry/Element.h>

namespace py = pybind11;

namespace caribou::geometry::bindings {

template<typename BezierSurfType>
void declare_bezier_surf(py::module & m, const std::string & name) {
    using BaseBezierSurfType = typename BezierSurfType::Base;

    // BaseBezierSurf
    std::string base_name = "Base" + name;
    declare_element<BezierSurfType>(m, base_name);
    py::class_<BaseBezierSurfType, Element<BezierSurfType>> (m, base_name.c_str());

    // BezierSurf
    py::class_<BezierSurfType, BaseBezierSurfType> c (m, name.c_str());
    c.def("__str__", [&name](const BezierSurfType & s) {
        Eigen::IOFormat f(std::numeric_limits<double>::digits10 + 2, 0, ", ", ", ", "[", "]", "[", "]");
        if constexpr(caribou::geometry::traits<BezierSurfType>::Dimension == 1) {
            f.rowPrefix = "";
            f.rowSuffix = "";
        }
        std::stringstream ss;
        ss << name << " : ";
        ss << s.nodes().format(f);
        return ss.str();
    });
}

void create_bezier_surf(pybind11::module & m) {

    declare_bezier_surf<BezierSurf<_2D>>(m, "BezierSurf_2D");

    declare_bezier_surf<BezierSurf<_3D>>(m, "BezierSurf_3D");

    m.def("BezierSurf", [](){
        return py::cast(BezierSurf<_2D>());
    });


    m.def("BezierSurf", [](const caribou::bindings::Dimension & dim) {
        if (dim == caribou::bindings::Dimension::_2D) {
            return py::cast(BezierSurf<_2D>());
        } else {
            return py::cast(BezierSurf<_3D>());
        }
    }, py::arg("dimension"));


    // 2D
    m.def("BezierSurf", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 9, 2> & nodes) {
        return py::cast(BezierSurf<_2D>(nodes));
    }, py::arg("nodes"));

    using bezier_surf_2d_node = Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1>;

    m.def("BezierSurf", [](const bezier_surf_2d_node & n0, const bezier_surf_2d_node & n1, const bezier_surf_2d_node & n2,
                           const bezier_surf_2d_node & n3, const bezier_surf_2d_node & n4, const bezier_surf_2d_node & n5,
                           const bezier_surf_2d_node & n6, const bezier_surf_2d_node & n7, const bezier_surf_2d_node & n8) {
        return py::cast(BezierSurf<_2D>(n0, n1, n2, n3, n4, n5, n6, n7, n8));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"),
                                                    py::arg("n6"), py::arg("n7"), py::arg("n8"));

    // 3D
    m.def("BezierSurf", [](const Eigen::Matrix<FLOATING_POINT_TYPE, 9, 3> & nodes) {
        return py::cast(BezierSurf<_3D>(nodes));
    }, py::arg("nodes"));

    using bezier_surf_3d_node = Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>;

    m.def("BezierSurf", [](const bezier_surf_3d_node & n0, const bezier_surf_3d_node & n1, const bezier_surf_3d_node & n2,
                           const bezier_surf_3d_node & n3, const bezier_surf_3d_node & n4, const bezier_surf_3d_node & n5,
                           const bezier_surf_3d_node & n6, const bezier_surf_3d_node & n7, const bezier_surf_3d_node & n8) {
        return py::cast(BezierSurf<_3D>(n0, n1, n2, n3, n4, n5, n6, n7, n8));
    }, py::arg("n0"), py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("n4"), py::arg("n5"),
                                                   py::arg("n6"), py::arg("n7"), py::arg("n8"));

}

}
