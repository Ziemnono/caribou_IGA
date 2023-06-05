#include <SofaCaribou/Python/Topology/CaribouSplineTopology.h>

#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>

#include <pybind11/pybind11.h>

namespace SofaCaribou::topology::python {

void addCaribouSplineTopology(pybind11::module &m) {
    using namespace caribou;
    using namespace caribou::geometry;

    bind_caribou_spline_topology<NurbsSurf<_2D>>(m);
    bind_caribou_spline_topology<NurbsSurf<_3D>>(m);

}

} // namespace SofaCaribou::topology::python
