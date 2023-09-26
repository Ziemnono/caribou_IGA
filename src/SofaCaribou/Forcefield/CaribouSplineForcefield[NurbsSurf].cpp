#include <SofaCaribou/Forcefield/CaribouSplineForcefield[NurbsSurf].h>
#include <SofaCaribou/Forcefield/CaribouSplineForcefield.inl>
#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::objectmodel;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// NurbsSurf linear specialization
// --------------------------------

// 2D
template <>
auto CaribouSplineForcefield<NurbsSurf<_2D>>::templateName(const CaribouSplineForcefield<NurbsSurf<_2D>> *) -> std::string {
    return SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_2D>>::templateName();
}

template <>
void CaribouSplineForcefield<NurbsSurf<_2D>>::triangulate_face(const NurbsSurf<_2D> & /*e*/, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & /*triangles_nodes*/)
{
    std::cout << "Trinagulate the faces 2D NURBS Surf \n";
//    int divs = 4;

//    int ncells = divs * divs;

//    using LocalCoordinates = NurbsSurf<_2D>::LocalCoordinates;


}

// This will force the compiler to compile the following templated class
template class CaribouSplineForcefield<NurbsSurf<_2D>>;

// 3D
template <>
auto CaribouSplineForcefield<NurbsSurf<_3D>>::templateName(const CaribouSplineForcefield<NurbsSurf<_3D>> *) -> std::string {
    return SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_3D>>::templateName();
}

template <>
void CaribouSplineForcefield<NurbsSurf<_3D>>::triangulate_face(const NurbsSurf<_3D> & /*e*/, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & /*triangles_nodes*/)
{
    std::cout << "Trinagulate the faces 3D NURBS Surf \n";
//    int divs = 4;

//    int ncells = divs * divs;

//    using LocalCoordinates = NurbsSurf<_2D>::LocalCoordinates;


}

// This will force the compiler to compile the following templated class
template class CaribouSplineForcefield<NurbsSurf<_3D>>;

} // namespace SofaCaribou::forcefield
