#pragma once
#include <SofaCaribou/Forcefield/HyperelasticSplineForcefield.h>
#include <SofaCaribou/Forcefield/CaribouSplineForcefield[NurbsSurf].h>
#include <Caribou/Geometry/NurbsSurf.h>

namespace SofaCaribou::forcefield {

// Hexahedron quadratic specialization
extern template class HyperelasticSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>>;

}
