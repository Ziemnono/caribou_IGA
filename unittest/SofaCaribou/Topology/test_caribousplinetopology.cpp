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

#include <SofaCaribou/Topology/CaribouSplineTopology[NurbsSurf].h>

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

TEST(CaribouSplineTopology, NS2D_AttachPatch_Rectangle) {

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
    std::cout << "\nINITIAL TEST ARE EXECUTED \n";



    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_2D>> *> (
            createObject(root, "CaribouSplineTopology", {{"template", "NurbsSurf_2D"}}).get() );
    EXPECT_NE(topo, nullptr);
    std::cout << "TOPO IS CREATED\n";
    // Attach the spline patch
    topo->attachSplinePatch(splinepatch);
    std::cout << "SPLINE PATCH ATTACHED \n";
    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 9>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
//    EXPECT_EQ(indices.size(), 16); // Number of elements
    std::cout << "Inidices starts : " << indices.size() << "\n";

    using real_val = SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_2D>>::Real;
    using Dataknots = Data<sofa::type::vector<sofa::type::fixed_array<real_val, 4>>>;
    auto knots = ReadAccessor<Dataknots> (dynamic_cast<Dataknots*>(topo->findData("knot_spans")));
    std::cout << "Knots starts : " << knots[0] << "\n";
    using Datawgts = Data<sofa::type::vector<real_val>>;
    auto weights = ReadAccessor<Datawgts> (dynamic_cast<Datawgts*>(topo->findData("weights")));
    std::cout << "Weights Start : " << weights[0] << " \n";

//    using Dataextrs = Data<sofa::type::vector<sofa::type::fixed_array<sofa::type::fixed_array<real_val, 9>, 9>>>;
//    auto extras = ReadAccessor<Dataextrs> (dynamic_cast<Dataextrs*>(topo->findData("extractions")));
//    std::cout << "\n Extraaction matrix Start : \n" << extras[0][0] << "\n"; // First row

    using vcord = SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_2D>>::VecCoord;
    using Datapositions = Data<vcord>;
    auto positions = ReadAccessor<Datapositions> (dynamic_cast<Datapositions*>(topo->findData("position")));
    std::cout << "Positions 1 : " << positions[0] << "\n"; // First row
    std::cout << "Positions 2 : " << positions[1] << "\n"; // First row
    std::cout << "Positions 3 : " << positions[2] << "\n"; // First row
    std::cout << "Positions 4 : " << positions[4] << "\n"; // First row
}

TEST(CaribouSplineTopology, NS2D_AttachPatch_PlateHole) {

    std::cout << "\n  ------------------- Plate Hole ----------------- \n";
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using SplinePatch = io::NURBSReader<_2D, PointID>::PatchType;
    std::string path = executable_directory_path + "/meshes/splines/plate_hole_geo.txt";
    auto reader = io::NURBSReader<_2D, PointID>::Read(path);

    const auto * splinepatch = dynamic_cast<const SplinePatch *>(reader.patch_ptr());
    EXPECT_EQ(splinepatch->number_of_nodes(), 12);
    EXPECT_EQ(splinepatch->number_of_nodes_per_elements(), 9);
    EXPECT_EQ(splinepatch->number_of_elements(), 2);
    std::cout << "\n  ---------------------- INITIAL TESTS ARE EXECUTED \n";



    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_2D>> *> (
            createObject(root, "CaribouSplineTopology", {{"template", "NurbsSurf_2D"}}).get() );
    EXPECT_NE(topo, nullptr);
    std::cout << "\n ----------------------  TOPO IS CREATED\n";
    // Attach the spline patch
    topo->attachSplinePatch(splinepatch);
    std::cout << "\n ---------------------- SPLINE PATCH ATTACHED \n";
    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 9>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
//    EXPECT_EQ(indices.size(), 16); // Number of elements
    std::cout << "Inidices starts : " << indices.size() << "\n";

    using real_val = SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_2D>>::Real;
    using Dataknots = Data<sofa::type::vector<sofa::type::fixed_array<real_val, 4>>>;
    auto knots = ReadAccessor<Dataknots> (dynamic_cast<Dataknots*>(topo->findData("knot_spans")));
    std::cout << "1st elem span : " << knots[0] << "\n";
    std::cout << "2nd elem span : " << knots[1] << "\n";

    using Datawgts = Data<sofa::type::vector<real_val>>;
    auto weights = ReadAccessor<Datawgts> (dynamic_cast<Datawgts*>(topo->findData("weights")));
    std::cout << "Weights size : " << weights.size() << "\n";
    std::cout << "Weights : \n";
    for (int i = 0; i < static_cast<int>(weights.size()); i++){
        std::cout << weights[i] << " ";
    }
    std::cout << "\n";

    auto knot_1 = ReadAccessor<Datawgts> (dynamic_cast<Datawgts*>(topo->findData("knot_1")));
    auto knot_2 = ReadAccessor<Datawgts> (dynamic_cast<Datawgts*>(topo->findData("knot_2")));

    std::cout << "Knots 1 : ";
    for (int i = 0; i < static_cast<int>(knot_1.size()); i++){
        std::cout << knot_1[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Knots 2 : ";
    for (int i = 0; i < static_cast<int>(knot_2.size()); i++){
        std::cout << knot_2[i] << " ";
    }
    std::cout << "\n";


    using vcord = SofaCaribou::topology::CaribouSplineTopology<NurbsSurf<_2D>>::VecCoord;
    using Datapositions = Data<vcord>;
    auto positions = ReadAccessor<Datapositions> (dynamic_cast<Datapositions*>(topo->findData("position")));
    std::cout << "Positions 1 : " << positions[0] << "\n";
    std::cout << "Positions 2 : " << positions[1] << "\n";
    std::cout << "Positions 3 : " << positions[2] << "\n";
    std::cout << "Positions 4 : " << positions[3] << "\n";
}


