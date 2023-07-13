#include <SofaCaribou/Forcefield/TractionSplineForcefield[NurbsSurf].h>
#include <SofaCaribou/Forcefield/TractionSplineForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// -----------------------------------
// Hexahedron quadratic specialization
// -----------------------------------

// This will force the compiler to compile the following templated class
template class TractionSplineForcefield<NurbsSurf<_2D>>;

} // namespace SofaCaribou::forcefield

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou elastic spline force field")
    .add<TractionSplineForcefield<NurbsSurf<_2D>>>();
}
