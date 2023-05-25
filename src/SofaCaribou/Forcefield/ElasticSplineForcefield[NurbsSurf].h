#pragma once
#include <SofaCaribou/Forcefield/ElasticSplineForcefield.h>
#include <SofaCaribou/Forcefield/CaribouSplineForcefield[NurbsSurf].h>
#include <Caribou/Geometry/NurbsSurf.h>

namespace SofaCaribou::forcefield {

// Hexahedron quadratic specialization
extern template class ElasticSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>>;

}
