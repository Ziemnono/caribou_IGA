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

TEST(ElasticSplineForcefield, NurbsSurf_from_SOFA) {
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
    createObject(root, "NurbsSurf", {{"name", "nurbs"}});

    auto meca = createChild(root, "meca");
    // Create the ODE system
    createObject(meca, "StaticODESolver", {{"newton_iterations", "10"}, {"correction_tolerance_threshold", "1e-5"}, {"residual_tolerance_threshold", "1e-5"}});
    createObject(meca, "LDLTSolver");
    createObject(meca, "MechanicalObject", {{"name", "mo"}, {"src", "@" + executable_directory_path + "/meshes/splines/xy_rectangle.txt"}});

//    // Complete Nurbs topology container
    createObject(meca, "CaribouSplineTopology", {{"name", "mechanical_topology"}, {"src", "@" + executable_directory_path + "/meshes/splines/xy_rectangle.txt"}});

//    // Mechanics
    createObject(meca, "SaintVenantKirchhoffMaterial", {{"young_modulus", "3000"}, {"poisson_ratio", "0.499"}});
    auto ff = dynamic_cast<SofaCaribou::forcefield::ElasticSplineForcefield<caribou::geometry::NurbsSurf<caribou::_2D>> *> (
        createObject(meca, "ElasticSplineForcefield").get()
    );

//    // Fix the left side of the beam
//    createObject(meca, "BoxROI", {{"name", "fixed_roi"}, {"quad", "@surface_topology.quad"}, {"box", "-7.5 -7.5 -0.9 7.5 7.5 0.1"}});
//    createObject(meca, "FixedConstraint", {{"indices", "@fixed_roi.indices"}});

//    // Apply traction on the right side of the beam
//    createObject(meca, "BoxROI", {{"name", "top_roi"}, {"quad", "@surface_topology.quad"}, {"box", "-7.5 -7.5 79.9 7.5 7.5 80.1"}});
//    createObject(meca, "QuadSetTopologyContainer", {{"name", "traction_container"}, {"quads", "@top_roi.quadInROI"}});
//    createObject(meca, "TractionForcefield", {{"traction", "0 -30 0"}, {"slope", "0.2"}, {"topology", "@traction_container"}});

    getSimulation()->init(root.get());

    EXPECT_EQ(ff->number_of_elements(), 2);
}
