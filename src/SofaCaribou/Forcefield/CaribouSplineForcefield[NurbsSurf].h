#pragma once
#include <SofaCaribou/Forcefield/CaribouSplineForcefield.h>
#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>
#include <Caribou/Geometry/NurbsSurf.h>

namespace SofaCaribou::forcefield {

// Quad linear specialization
// 2D
template <> auto CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>>::templateName(const CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>> *) -> std::string;
extern template class CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>>;

// 3D
template <> auto CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_3D>>::templateName(const CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_3D>> *) -> std::string;
extern template class CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_3D>>;

}
