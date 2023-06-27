#include "HyperelasticSplineForcefield.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/NurbsSurf.h>

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addHyperelasticSplineForcefield(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_hyper_elastic_spline_forcefield<NurbsSurf<_2D>>(m);
//    bind_elastic_spline_forcefield<NurbsSurf<_3D>>(m);

}

}
