#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Caribou/constants.h>
#include <Caribou/Topology/SplinePatch.h>

namespace py = pybind11;

namespace caribou::topology::bindings {

template < UNSIGNED_INTEGER_TYPE Dim, typename NodeIndex, typename MatrixType>
void declare_splinepatch(py::module & m){
    using M = SplinePatch<Dim, NodeIndex, EigenSplineNodesHolder<MatrixType>>;

    std::string name = typeid(M).name();
    py::class_<M, BaseSplinePatch> c(m, name.c_str());

    c.def("dimension", &M::dimension);
    c.def("number_of_nodes", &M::number_of_nodes);
    c.def("canonical_dimension", &M::canonical_dimension);
    c.def("number_of_nodes_per_elements", &M::number_of_nodes_per_elements);
    c.def("number_of_elements", &M::number_of_elements);
    c.def("print_spline_patch", &M::print_spline_patch);

    // Positions

    c.def("positions", [](const M & m, const std::vector<INTEGER_TYPE> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());
    c.def("positions", [](const M & m, const std::vector<int> & indices) {
        return m.positions(indices);
    }, py::arg("indices").noconvert());

    // Spline functions

    // Weights
//    c.def("weights", [](const S & s, const std::vector<UNSIGNED_INTEGER_TYPE> & indices) {
//        return s.weights(indices);
//    }. py::arg("indices").noconvert());

    c.def("weights", [](const M & m, const std::vector<INTEGER_TYPE> & indices) {
        return m.weights(indices);
    }, py::arg("indices").noconvert());
    c.def("weights", [](const M & m, const std::vector<int> & indices) {
        return m.weights(indices);
    }, py::arg("indices").noconvert());

    // knot_1 related functions
    c.def("size_knot_1", &M::size_knot_1);
    c.def("knot_1", &M::knot_1);

    // knot_2 related functions
    c.def("size_knot_2", &M::size_knot_2);
    c.def("knot_2", &M::knot_2);

    // Element knot ranges
    c.def("element_knotranges", [](const M & m, const UNSIGNED_INTEGER_TYPE & index){
        return m.element_knotranges(index);
    }, py::arg("index"));

    // Constructions
    using DynVector = typename M::DynVector;
    using IndicesMatrix = Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using Double_Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using ElementsIndices = typename M::ElementsIndices;
    using ElementsKnotrange = typename M::ElementsKnotrange;

    m.def("SplinePatch", [](const MatrixType & positions, const DynVector & weights, const ElementsIndices & indices,
                            const DynVector & knot_1, const DynVector & knot_2, const ElementsKnotrange & knotrange) {
        return py::cast(M(positions, weights, indices, knot_1, knot_2, knotrange));
    }, py::arg("positions").noconvert(), py::arg("weights").noconvert(), py::arg("indices").noconvert(),
       py::arg("knot_1").noconvert(), py::arg("knot_2").noconvert(), py::arg("knotrange").noconvert());

    m.def("SplinePatch", [](const Double_Matrix & positions, const DynVector & weights, const IndicesMatrix & indices,
                            const DynVector & knot_1, const DynVector & knot_2, const Double_Matrix & knotrange) {
        return py::cast(M(positions, weights, indices, knot_1, knot_2, knotrange));
    }, py::arg("positions").noconvert(), py::arg("weights").noconvert(), py::arg("indices").noconvert(),
       py::arg("knot_1").noconvert(), py::arg("knot_2").noconvert(), py::arg("knotrange").noconvert());

}

template <UNSIGNED_INTEGER_TYPE Dim>
void declare_splinepatch(py::module &m) {
    declare_splinepatch<Dim, UNSIGNED_INTEGER_TYPE, Eigen::Matrix<double, Eigen::Dynamic, Dim, (Dim>1?Eigen::RowMajor:Eigen::ColMajor)>>(m);
}

void create_splinepatch(py::module & m) {
    py::class_<BaseSplinePatch> a (m, "BaseSplinePatch");

    declare_splinepatch<2>(m);
//    declare_splinepatch<3>(m);
}

}
