#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
#include <sofa/helper/testing/BaseTest.h>
#else
#include <sofa/testing/BaseTest.h>
#endif
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaSimulationGraph/DAGSimulation.h>
#include <SofaSimulationGraph/SimpleApi.h>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/simulation/Node.h>
DISABLE_ALL_WARNINGS_BEGIN

#include <SofaCaribou/Topology/CaribouSplineTopology[BezierSurf].h>

#include <Caribou/Topology/SplinePatch.h>
#include <Caribou/Topology/IO/NURBSReader.h>
#include <Caribou/Topology/IO/NURBSReader.cpp>

#include "../sofacaribou_test.h"

using sofa::helper::system::PluginManager ;
using namespace sofa::simulation;
using namespace sofa::simpleapi;
using namespace sofa::helper::logging;

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
using namespace sofa::helper::testing;
#else
using namespace sofa::testing;
#endif

TEST(CaribouSplineTopology, BezierSurf2DAttachPatch) {

    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using SplinePatch = io::NURBSReader<_2D, PointID>::PatchType;
    std::string path = executable_directory_path + "/meshes/splines/xy_rectangle.txt";
    auto reader = io::NURBSReader<_2D, PointID>::Read(path);

    const auto * splinepatch = dynamic_cast<const SplinePatch *>(reader.patch_ptr());
    EXPECT_EQ(splinepatch->number_of_nodes(), 9);
    EXPECT_EQ(splinepatch->number_of_nodes_per_elements(), 9);
    EXPECT_EQ(splinepatch->number_of_elements(), 1);


//    EXPECT_EQ(splinepatch->number_of_domains(), 2);

//    // First domain is the segment contour, second is the quad surface domain
//    // Get the second one
//    const auto * domain = dynamic_cast<const Domain * >(splinepatch->domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
//    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
//    auto topo = dynamic_cast<SofaCaribou::topology::CaribouSplineTopology<BezierSurf<_2D>> *> (
//            createObject(root, "CaribouSplineTopology", {{"template", "BezierSurf_2D"}}).get() );
//    EXPECT_NE(topo, nullptr);

    // Attach the spline patch
//    topo->attachSplinePatch(splinepatch);

//    // Make sure the data parameter `indices` has been filled-up correctly
//    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 9>>>;
//    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
////    EXPECT_EQ(indices.size(), 16); // Number of elements

}


