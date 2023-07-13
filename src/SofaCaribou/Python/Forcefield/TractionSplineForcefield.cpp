#include "TractionSplineForcefield.h"

#include <Caribou/constants.h>
#include <Caribou/Geometry/NurbsSurf.h>

#include <pybind11/eigen.h>

namespace SofaCaribou::forcefield::python {

void addTractionSplineForcefield(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;
    bind_traction_spline_forcefield<NurbsSurf<_2D>>(m);

}

}
