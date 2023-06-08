#include <SofaCaribou/config.h>
#include <SofaCaribou/Material/PlaneStressMaterial.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using namespace sofa::core;

namespace SofaCaribou::material {

static int PlaneStressClass = RegisterObject("Caribou PlaneStressMaterial elastic material")
//    .add< NeoHookeanMaterial<sofa::defaulttype::Vec1Types> >()
    .add< PlaneStressMaterial<sofa::defaulttype::Vec2Types> >();
//    .add< NeoHookeanMaterial<sofa::defaulttype::Vec3Types> >();

} // namespace SofaCaribou::material
