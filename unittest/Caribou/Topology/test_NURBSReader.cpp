#include <gtest/gtest.h>
#include <Caribou/Topology/IO/NURBSReader.cpp>
#include <string>
#include <cmath>
#include "topology_test.h"

using namespace caribou;
using namespace caribou::topology;
using namespace caribou::geometry;


TEST(NURBSReader, knot_test_geo){

//    using Patch = io::NURBSReader<_2D>::PatchType;
    std::string path = executable_directory_path + "/meshes/splines/knot_test_geo.txt";
    auto reader = io::NURBSReader<_2D>::Read(path);
    auto patch = reader.patch();
//    std::cout << "\n No. of nodes knot -> " << patch.number_of_nodes() << "\n";
    EXPECT_EQ(patch.number_of_nodes(), 15);

}

TEST(NURBSReader, xy_rectangle){
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

TEST(NURBSReader, quarter_cyl){
    std::string path = executable_directory_path + "/meshes/splines/quarter_cylinder.txt";
    auto reader = io::NURBSReader<_2D>::Read(path);
    auto patch = reader.patch();

    EXPECT_EQ(patch.number_of_nodes(), 12);

//    patch.element(0).print();


}

TEST(NURBSReader, plate_hole){
    std::string path = executable_directory_path + "/meshes/splines/plate_hole_geo.txt";
    using DynVector = io::NURBSReader<_2D>::DynVector;

    auto reader = io::NURBSReader<_2D>::Read(path);
    auto patch = reader.patch();

    EXPECT_EQ(patch.number_of_nodes(), 12);
    EXPECT_EQ(patch.number_of_elements(), 2);
    EXPECT_EQ(patch.number_of_nodes_per_elements(), 9);


    DynVector wgts(12);
    wgts << 1.000000000000000,   0.853553390593274,   0.853553390593274,  1.000000000000000,   1.000000000000000,
            1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.000000000000000,
            1.000000000000000,   1.000000000000000 ;

    EXPECT_MATRIX_EQUAL(patch.weights(), wgts);

    DynVector span1(4);
    span1 << 0, 0, 0.5, 1;
    DynVector span2(4);
    span2 << 0.5, 0, 1, 1;

    EXPECT_MATRIX_EQUAL(patch.element_knotranges(0), span1.transpose());
    EXPECT_MATRIX_EQUAL(patch.element_knotranges(1), span2.transpose());


}

TEST(NURBSReader, plate_hole_integration){
    std::string path = executable_directory_path + "/meshes/splines/plate_hole_geo.txt";
//    using DynVector = io::NURBSReader<_2D>::DynVector;

    auto reader = io::NURBSReader<_2D>::Read(path);
    auto patch = reader.patch();

    EXPECT_EQ(patch.number_of_nodes(), 12);
    EXPECT_EQ(patch.number_of_elements(), 2);
    EXPECT_EQ(patch.number_of_nodes_per_elements(), 9);

    Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> tval;
    tval << 0, 0;

    for (UNSIGNED_INTEGER_TYPE i = 0; i < patch.number_of_elements(); ++i) {
        auto nurbs_elem = patch.element(i);
        for (const auto & g : nurbs_elem.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto detJ = nurbs_elem.jacobian(x).determinant();
            auto jp2p  = nurbs_elem.jacobian_parent2para(x);
            const auto N = nurbs_elem.world_coordinates(x);
            tval = tval + N * jp2p * detJ * w;
        }
    }

    Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> matlab_val;

    EXPECT_EQ(std::round(tval[0] * 100.0) /100.0, -31.67);
    EXPECT_EQ(std::round(tval[1]* 100.0) /100.0, 31.67);

//    std::cout << "Integration values : " << tval.transpose() << "\n";

}



