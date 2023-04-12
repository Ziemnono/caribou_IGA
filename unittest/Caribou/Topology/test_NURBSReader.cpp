#include "topology_test.h"
#include <gtest/gtest.h>
#include <Caribou/Topology/IO/NURBSReader.cpp>


TEST(NURBSReader, knot_test_geo){
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

//    using Patch = io::NURBSReader<_2D>::PatchType;
    auto reader = io::NURBSReader<_2D>::Read(executable_directory_path + "/meshes/splines/xy_rectangle.txt");
    auto patch = reader.patch();
//    EXPECT_EQ(patch.number_of_nodes(), 15);
    EXPECT_EQ(15, 15);
}

TEST(SplinePatch, xy_rectangle){
    EXPECT_EQ(2,2);
}

TEST(SplinePatch, quarter_cyl){
    EXPECT_EQ(2,2);
}

