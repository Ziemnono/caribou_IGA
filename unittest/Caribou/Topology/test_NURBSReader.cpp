#include <gtest/gtest.h>
#include <Caribou/Topology/IO/NURBSReader.cpp>
#include <string>
#include "topology_test.h"

using namespace caribou;
using namespace caribou::topology;
using namespace caribou::geometry;

TEST(NURBSReader, knot_test_geo){

//    using Patch = io::NURBSReader<_2D>::PatchType;
    std::string path = executable_directory_path + "/meshes/splines/knot_test_geo.txt";
    auto reader = io::NURBSReader<_2D>::Read(path);
    auto patch = reader.patch();
    std::cout << "\n No. of nodes knot -> " << patch.number_of_nodes() << "\n";
    EXPECT_EQ(patch.number_of_nodes(), 15);

}

TEST(SplinePatch, xy_rectangle){
//    using Patch = io::NURBSReader<_2D>::PatchType;
    std::string path = executable_directory_path + "/meshes/splines/xy_rectangle.txt";
    auto reader = io::NURBSReader<_2D>::Read(path);
    auto patch = reader.patch();
    std::cout << "\n No. of nodes xy -> " << patch.number_of_nodes() << "\n";
    EXPECT_EQ(patch.number_of_nodes(), 9);
}

TEST(SplinePatch, quarter_cyl){
    EXPECT_EQ(2,2);
}

