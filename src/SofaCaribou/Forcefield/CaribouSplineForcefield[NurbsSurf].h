#pragma once
#include <SofaCaribou/Forcefield/CaribouSplineForcefield.h>
#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>
#include <Caribou/Geometry/NurbsSurf.h>

namespace SofaCaribou::forcefield {

// Quad linear specialization
// 2D
template <> auto CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>>::templateName(const CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>> *) -> std::string;
template <> void CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>>::triangulate_face(const caribou::geometry::NurbsSurf<caribou::_2D> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>>;

// 3D
template <> auto CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_3D>>::templateName(const CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_3D>> *) -> std::string;
template <> void CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_3D>>::triangulate_face(const caribou::geometry::NurbsSurf<caribou::_3D> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouSplineForcefield<caribou::geometry::NurbsSurf<caribou::_3D>>;

}
