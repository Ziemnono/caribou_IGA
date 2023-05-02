#pragma once

#include <SofaCaribou/Topology/CaribouSplineTopology.h>
#include <Caribou/Geometry/BezierSurf.h>

namespace SofaCaribou::topology {

// BezierSurf 2D linear specialization
template<>
auto CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_2D>>::templateName(
        const CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_2D>> *) -> std::string;

//template <> auto CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_2D>>::mesh_is_compatible(
//        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

//template <> auto CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_2D>>::get_indices_data_from(
//        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_2D>>;

// BezierSurf 3D linear specialization
template<>
auto CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_3D>>::templateName(
        const CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_3D>> *) -> std::string;

//template <> auto CaribouSplineTopology<caribou::geometry::BezierSurf <caribou::_3D>>::mesh_is_compatible(
//        const sofa::core::topology::BaseMeshTopology * topology) -> bool;

//template <> auto CaribouSplineTopology<caribou::geometry::BezierSurf <caribou::_3D>>::get_indices_data_from(
//        const sofa::core::topology::BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData *;

extern template
class CaribouSplineTopology<caribou::geometry::BezierSurf<caribou::_3D>>;

} // namespace SofaCaribou::topology
