#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <SofaCaribou/Forcefield/TractionSplineForcefield.h>
#include <SofaCaribou/Forcefield/TractionSplineForcefield.inl>

namespace SofaCaribou::forcefield::python {

template<typename Element>
void bind_traction_spline_forcefield(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "TractionSplineForcefield<" + template_name + ">";

    pybind11::class_<TractionSplineForcefield<Element>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<TractionSplineForcefield<Element>>> c(m, name.c_str());

    c.def("apply_load", &TractionSplineForcefield<Element>::apply_load );

    sofapython3::PythonFactory::registerType<TractionSplineForcefield<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<TractionSplineForcefield<Element>*>(o));
    });
}

template<typename Element>
void bind_traction_spline_forcefield(pybind11::module &m) {
    bind_traction_spline_forcefield<Element>(m, TractionSplineForcefield<Element>::templateName());
}

void addTractionSplineForcefield(pybind11::module &m);
}
