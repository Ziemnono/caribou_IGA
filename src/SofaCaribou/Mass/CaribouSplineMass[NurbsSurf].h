#pragma once

#include <SofaCaribou/Mass/CaribouSplineMass.h>
#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>
#include <Caribou/Geometry/NurbsSurf.h>
using namespace caribou;

namespace SofaCaribou::mass {

// NurbsSurf specialization
extern template class CaribouSplineMass<caribou::geometry::NurbsSurf<_2D>>;

} // namespace SofaCaribou::mass
