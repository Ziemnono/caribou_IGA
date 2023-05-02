#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouSplineTopology.inl>
#include <SofaCaribou/Topology/CaribouSplineTopology[BezierSurf].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/QuadSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace sofa::component::topology;
using namespace sofa::core::topology;
using namespace caribou;

namespace SofaCaribou::topology {

// BezierSurf 2D linear specialization
template<>
auto CaribouSplineTopology<BezierSurf<_2D>>::templateName(const CaribouSplineTopology<BezierSurf<_2D>> *) -> std::string {
    return "BezierSurf_2D";
}

//template <>
//auto CaribouSplineTopology<BezierSurf<_2D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
//    return (
//    dynamic_cast<const BezierSurfSetTopologyContainer*>(topology) != nullptr
//    );
//}

//template <>
//auto CaribouSplineTopology<BezierSurf<_2D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
//    return topology->findData("BezierSurfs");
//}

// BezierSurf 3D linear specialization
template<>
auto
CaribouSplineTopology<BezierSurf<_3D>>::templateName(const CaribouSplineTopology<BezierSurf<_3D>> *) -> std::string {
    return "BezierSurf";
}

//template <>
//auto CaribouSplineTopology<BezierSurf<_3D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
//    return (
//            dynamic_cast<const BezierSurfSetTopologyContainer*>(topology) != nullptr
//    );
//}

//template <>
//auto CaribouSplineTopology<BezierSurf<_3D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
//    return topology->findData("BezierSurfs");
//}

// This will force the compiler to compile the class with some template type
template class CaribouSplineTopology<BezierSurf<_2D>>;
template class CaribouSplineTopology<BezierSurf<_3D>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou topology")
.add<CaribouSplineTopology<BezierSurf<_2D>>>()
.add<CaribouSplineTopology<BezierSurf<_3D>>>()
;
}
