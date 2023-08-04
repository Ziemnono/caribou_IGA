#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouSplineTopology.inl>
#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>

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

// NurbsSurf 2D linear specialization
template<>
auto CaribouSplineTopology<NurbsSurf<_2D>>::templateName(const CaribouSplineTopology<NurbsSurf<_2D>> *) -> std::string {
    return "NurbsSurf_2D";
}

//template <>
//auto CaribouSplineTopology<NurbsSurf<_2D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
//    return (
//    dynamic_cast<const NurbsSurfSetTopologyContainer*>(topology) != nullptr
//    );
//}

//template <>
//auto CaribouSplineTopology<NurbsSurf<_2D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
//    return topology->findData("NurbsSurfs");
//}

// NurbsSurf 3D linear specialization
template<>
auto
CaribouSplineTopology<NurbsSurf<_3D>>::templateName(const CaribouSplineTopology<NurbsSurf<_3D>> *) -> std::string {
    return "NurbsSurf";
}

//template <>
//auto CaribouSplineTopology<NurbsSurf<_3D>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
//    return (
//            dynamic_cast<const NurbsSurfSetTopologyContainer*>(topology) != nullptr
//    );
//}

//template <>
//auto CaribouSplineTopology<NurbsSurf<_3D>>::get_indices_data_from(const BaseMeshTopology * topology) -> sofa::core::objectmodel::BaseData * {
//    return topology->findData("NurbsSurfs");
//}

// This will force the compiler to compile the class with some template type
template class CaribouSplineTopology<NurbsSurf<_2D>>;
template class CaribouSplineTopology<NurbsSurf<_3D>>;

} // namespace SofaCaribou::topology

namespace sofa::core::objectmodel {
using namespace SofaCaribou::topology;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou spline topology")
.add<CaribouSplineTopology<NurbsSurf<_2D>>>()
.add<CaribouSplineTopology<NurbsSurf<_3D>>>()
;
}
