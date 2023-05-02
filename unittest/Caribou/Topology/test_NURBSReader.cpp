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
//    std::cout << "\n ############  START -- XY-Rectangle  ########## \n";
//    std::cout << "\n No. of nodes xy -> " << patch.number_of_nodes() << "\n";
    EXPECT_EQ(patch.number_of_nodes(), 9);
//    std::cout << "No. of elements \n" << patch.number_of_elements() << "\n this is end \n";
//    std::cout << "Element Knot range 1 \n" << patch.element_knotranges(0) << "\n this is end \n";
//    std::cout << "\n ############  END -- XY-Rectangle  ########## \n";
}

TEST(SplinePatch, quarter_cyl){
    std::string path = executable_directory_path + "/meshes/splines/quarter_cylinder.txt";
    auto reader = io::NURBSReader<_2D>::Read(path);
    auto patch = reader.patch();

    EXPECT_EQ(patch.number_of_nodes(), 12);


}


TEST(SplinePatch, quarter_cyl_ptr){
    using PatchType = io::NURBSReader<_2D>::PatchType;
    std::string path = executable_directory_path + "/meshes/splines/quarter_cylinder.txt";
    auto reader = io::NURBSReader<_2D>::Read(path);
    const PatchType * patch = reader.patch_ptr();
    std::cout << "\n ############  START -- UNIQUE_PTR  ########## \n";
    std::cout << "\n No. of nodes xy -> " << patch->number_of_nodes() << "\n";
//    EXPECT_EQ(patch.number_of_nodes(), 9);
    std::cout << "No. of elements \n" << patch->number_of_elements() << "\n this is end \n";
    std::cout << "No. og nodes \n" << patch->number_of_nodes() << "\n this is end \n";
    std::cout << "Element Knot range 1 \n" << patch->element_knotranges(0) << "\n this is end \n";
    std::cout << "Element Knot range 2 \n" << patch->element_knotranges(1) << "\n this is end \n";
    std::cout << "Patch weights  \n" << patch->weights({0,1,2,3,4,5,6,7,8,9,10,11}) << "\n this is end \n";
    std::cout << "\n ############  END -- UNIQUE_PTR  ########## \n";
//    free (patch);
}

