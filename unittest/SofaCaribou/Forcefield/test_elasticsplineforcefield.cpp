#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
#include <sofa/helper/testing/BaseTest.h>
#else
#include <sofa/testing/BaseTest.h>
#endif
#include <sofa/simulation/Node.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaSimulationGraph/DAGSimulation.h>
#include <SofaSimulationGraph/SimpleApi.h>
#include <sofa/helper/system/PluginManager.h>
DISABLE_ALL_WARNINGS_END

#include <SofaCaribou/Forcefield/ElasticSplineForcefield.h>
#include <SofaCaribou/Forcefield/ElasticSplineForcefield[NurbsSurf].h>
//#include <Caribou/Topology/IO/NURBSReader.cpp>

#include "../sofacaribou_test.h"

using sofa::helper::system::PluginManager ;
using namespace sofa::simulation;
using namespace sofa::simpleapi;
using namespace sofa::helper::logging;

using namespace caribou;
using namespace caribou::topology;

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
using namespace sofa::helper::testing;
#else
using namespace sofa::testing;
#endif

TEST(ElasticSplineForcefield, NurbsSurf_from_SOFA) {

    std::string path = executable_directory_path + "/meshes/splines/knot_test_geo.txt";
//    auto reader = io::NURBSReader<_2D>::Read(path);
//    auto patch = reader.patch();

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 201200)
    createObject(root, "RequiredPlugin", {{"pluginName", "SofaBoundaryCondition SofaEngine"}});
#else
    createObject(root, "RequiredPlugin", {{"pluginName", "SofaComponentAll"}});
#endif
#if (defined(SOFA_VERSION) && SOFA_VERSION > 201299)
    createObject(root, "RequiredPlugin", {{"pluginName", "SofaTopologyMapping"}});
#endif
    auto meca = createChild(root, "meca");
    // Create the ODE system
    createObject(meca, "StaticODESolver", {{"newton_iterations", "10"}, {"correction_tolerance_threshold", "1e-5"}, {"residual_tolerance_threshold", "1e-5"}});
    createObject(meca, "LDLTSolver");

    // Complete Nurbs topology container
//    createObject(meca, "CaribouSplineTopology", name =  "topology", indices = patch.indices());



}
