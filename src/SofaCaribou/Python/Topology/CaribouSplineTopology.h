#pragma once
#include <pybind11/pybind11.h>
#include <SofaCaribou/Topology/CaribouSplineTopology.inl>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

namespace SofaCaribou::topology::python {

template<typename Element>
void bind_caribou_spline_topology(pybind11::module &m) {
    using sofa::core::objectmodel::BaseObject;
    const std::string name = "CaribouSplineTopology<" + CaribouSplineTopology<Element>::templateName() + ">";
    pybind11::class_<CaribouSplineTopology<Element>, BaseObject, sofapython3::py_shared_ptr<CaribouSplineTopology<Element>>> c(m, name.c_str());

    c.def("splinepatch", &CaribouSplineTopology<Element>::splinepatch);
}

void addCaribouSplineTopology(pybind11::module &m);

}
