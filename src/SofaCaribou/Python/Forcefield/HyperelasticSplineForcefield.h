#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <SofaCaribou/Forcefield/HyperelasticSplineForcefield.h>
#include <SofaCaribou/Forcefield/HyperelasticSplineForcefield.inl>

namespace SofaCaribou::forcefield::python {

template<typename Element>
void bind_hyper_elastic_spline_forcefield(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "HyperelasticSplineForcefield<" + template_name + ">";

    pybind11::class_<HyperelasticSplineForcefield<Element>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<HyperelasticSplineForcefield<Element>>> c(m, name.c_str());

    c.def("K", &HyperelasticSplineForcefield<Element>::K);
    c.def("cond", &HyperelasticSplineForcefield<Element>::cond);
    c.def("eigenvalues", &HyperelasticSplineForcefield<Element>::eigenvalues);
    c.def("assemble_stiffness", [](HyperelasticSplineForcefield<Element> & self, const Eigen::Matrix<double, Eigen::Dynamic, HyperelasticSplineForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));
    c.def("assemble_stiffness", [](HyperelasticSplineForcefield<Element> & self, const Eigen::Matrix<float, Eigen::Dynamic, HyperelasticSplineForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));

    sofapython3::PythonFactory::registerType<HyperelasticSplineForcefield<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<HyperelasticSplineForcefield<Element>*>(o));
    });
}

template<typename Element>
void bind_hyper_elastic_spline_forcefield(pybind11::module &m) {
    bind_hyper_elastic_spline_forcefield<Element>(m, HyperelasticSplineForcefield<Element>::templateName());
}

void addHyperelasticSplineForcefield(pybind11::module &m);
}
