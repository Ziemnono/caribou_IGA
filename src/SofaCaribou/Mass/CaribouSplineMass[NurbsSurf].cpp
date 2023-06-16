#include <SofaCaribou/Mass/CaribouSplineMass[NurbsSurf].h>
#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>
#include <SofaCaribou/Mass/CaribouSplineMass.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::objectmodel;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::mass {

// This will force the compiler to compile the following templated class
template class CaribouSplineMass<NurbsSurf<_2D>>;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou mass")
        .add<CaribouSplineMass<NurbsSurf<_2D>>>();

} // namespace SofaCaribou::mass
