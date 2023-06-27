#include <SofaCaribou/config.h>
#include <SofaCaribou/Material/SaintVenantKirchhoff2DMaterial.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using namespace sofa::core;

namespace SofaCaribou::material {

static int SaintVenantKirchhoff2DClass = RegisterObject("Caribou Saint-Venant-Kirchhoff hyperelastic material")
                                                 .add< SaintVenantKirchhoff2DMaterial<sofa::defaulttype::Vec2Types> >();

} // namespace SofaCaribou::material
