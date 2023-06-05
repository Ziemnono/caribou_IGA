#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <SofaCaribou/Forcefield/ElasticSplineForcefield.h>
#include <SofaCaribou/Forcefield/ElasticSplineForcefield.inl>

namespace SofaCaribou::forcefield::python {

template<typename Element>
void bind_elastic_spline_forcefield(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "ElasticSplineForcefield<" + template_name + ">";

    pybind11::class_<ElasticSplineForcefield<Element>, sofa::core::objectmodel::BaseObject, sofapython3::py_shared_ptr<ElasticSplineForcefield<Element>>> c(m, name.c_str());

    c.def("K", &ElasticSplineForcefield<Element>::K);
    c.def("cond", &ElasticSplineForcefield<Element>::cond);
    c.def("eigenvalues", &ElasticSplineForcefield<Element>::eigenvalues);
    c.def("assemble_stiffness", [](ElasticSplineForcefield<Element> & self, const Eigen::Matrix<double, Eigen::Dynamic, ElasticSplineForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));
    c.def("assemble_stiffness", [](ElasticSplineForcefield<Element> & self, const Eigen::Matrix<float, Eigen::Dynamic, ElasticSplineForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));

    sofapython3::PythonFactory::registerType<ElasticSplineForcefield<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<ElasticSplineForcefield<Element>*>(o));
    });
}

template<typename Element>
void bind_elastic_spline_forcefield(pybind11::module &m) {
    bind_elastic_spline_forcefield<Element>(m, ElasticSplineForcefield<Element>::templateName());
}

void addElasticSplineForcefield(pybind11::module &m);
}
