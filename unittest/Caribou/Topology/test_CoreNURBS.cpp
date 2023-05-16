#include "topology_test.h"
#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/IO/CoreNURBS.h>
//#include <Caribou/Topology/IO/CoreNURBS_1.h>
//#include <Caribou/Topology/SplinePatch.h>
//#include <Caribou/Topology/IO/NURBSReader.h>
//#include <Caribou/Topology/IO/trial_io.h>

using NodeIndex = UNSIGNED_INTEGER_TYPE;

TEST(CoreNURBS, para_topo) {
    // Testing para_topo data structure.
    /*
     * p = 2; q = 2;
     * knot_u = [0,0,0,0.5,1,1,1]
     * knot_v = [0,0,0,1,1,1]
     * +----------+----------+ v = 1
       |          |          |
       |          |          |
       |          |          |
       |          |          |
       |          |          |
       |          |          |
       +----------+----------+ v = 0
       0        u = 0.5      1

    element knot ranges will be [ [0, 0, 0.5, 1], [0.5, 0, 1, 1] ]
    element wise node indices / connectivity information
    [ [0, 1, 2, 4, 5, 6, 8, 9,  10]
      [1, 2, 3, 5, 6, 7, 9, 10, 11] ]
*/
    using namespace caribou;
    using namespace caribou::topology::io;
    Double_Matrix element_ranges(2,4);
    element_ranges << 0, 0, 0.5, 1,
                      0.5, 0, 1, 1;
    Matrix<NodeIndex> element_connectivity(2, 9);
    element_connectivity << 0, 1, 2, 4, 5, 6, 8, 9, 10,
                            1, 2, 3, 5, 6, 7, 9, 10, 11;

    utils::para_topo<NodeIndex> topo_info(element_ranges, element_connectivity);
    EXPECT_MATRIX_EQUAL(topo_info.get_elrange(), element_ranges);
    EXPECT_MATRIX_EQUAL(topo_info.get_elconn(), element_connectivity);
}

TEST(CoreNURBS, reader_knot) {
    using namespace caribou;
    using namespace caribou::topology::io;
    coreNurbs<NodeIndex> nurbs_patch;
    nurbs_patch.SetFileName(executable_directory_path + "/meshes/splines/knot_test_geo.txt");
    nurbs_patch.Update();

    EXPECT_EQ(nurbs_patch.GetP(), 2);
    EXPECT_EQ(nurbs_patch.GetQ(), 2);
    EXPECT_EQ(nurbs_patch.get_no_pnts_u(), 5);
    EXPECT_EQ(nurbs_patch.get_no_pnts_v(), 3);
    EXPECT_EQ(nurbs_patch.GetNumberOfPoints(), 15);
    EXPECT_EQ(nurbs_patch.get_no_elems_u(), 3);
    EXPECT_EQ(nurbs_patch.get_no_elems_v(), 1);
    EXPECT_EQ(nurbs_patch.GetNumberOfElements(), 3);

    // Knot vector
    Double_Vector knot_u(8), knot_v(6);
    knot_u << 0, 0, 0, 0.33, 0.66, 1, 1, 1;
    knot_v << 0, 0, 0, 1, 1, 1;
    EXPECT_MATRIX_EQUAL(nurbs_patch.get_knot_u(), knot_u);
    EXPECT_MATRIX_EQUAL(nurbs_patch.get_knot_v(), knot_v);
}


TEST(CoreNURBS, reader_rectangle) {
    using namespace caribou;
    using namespace caribou::topology::io;
    coreNurbs<NodeIndex> nurbs_patch;
    nurbs_patch.SetFileName(executable_directory_path + "/meshes/splines/xy_rectangle.txt");
    nurbs_patch.Update();

    EXPECT_EQ(nurbs_patch.GetP(), 2);
    EXPECT_EQ(nurbs_patch.GetQ(), 2);

    // Points
    Double_Matrix pnts(9,2);
    pnts << 0, 0,
            1, 0,
            2, 0,
            0, 1,
            1, 1,
            2, 1,
            0, 2,
            1, 2,
            2, 2;
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetPoints(), pnts);
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetPoint(4), pnts.row(4));
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetPoint(6), pnts.row(6));
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetPoint(8), pnts.row(8));

}

TEST(CoreNURBS, quarter_cylinder) {
    using namespace caribou;
    using namespace caribou::topology::io;
    coreNurbs<NodeIndex> nurbs_patch;
    nurbs_patch.SetFileName(executable_directory_path + "/meshes/splines/quarter_cylinder.txt");
    nurbs_patch.Update();

    EXPECT_EQ(nurbs_patch.GetP(), 2);
    EXPECT_EQ(nurbs_patch.GetQ(), 2);

    // Weights
    Double_Vector wgts(12);
    wgts << 1.0000,0.853553390593274,0.853553390593274,1.0000, 1.0000,0.853553390593274,0.853553390593274,1.0000, 1.0000,0.853553390593274,0.853553390593274,1.0000 ;

    EXPECT_MATRIX_EQUAL(nurbs_patch.GetWeights(), wgts);
    EXPECT_EQ(nurbs_patch.GetWeight(0), wgts(0));
    EXPECT_EQ(nurbs_patch.GetWeight(1), wgts(1));

    // Indices - Unsigned integer type
    USInt_Matrix indices(nurbs_patch.GetNumberOfElements(), nurbs_patch.GetNumberOfElementPoints());
    indices << 0, 1, 2, 4, 5, 6, 8, 9, 10,
               1, 2, 3, 5, 6, 7, 9, 10, 11;
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetIndices(), indices);

    // KnotRanges
    Double_Matrix  knotranges(nurbs_patch.GetNumberOfElements(), 4);
    knotranges << 0, 0, 0.5, 1,
                  0.5, 0, 1, 1;
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetKnotRanges(), knotranges);

}

TEST(CoreNURBS, plate_hole) {
    using namespace caribou;
    using namespace caribou::topology::io;
    coreNurbs<NodeIndex> nurbs_patch;
    nurbs_patch.SetFileName(executable_directory_path + "/meshes/splines/plate_hole_geo.txt");
    nurbs_patch.Update();

    EXPECT_EQ(nurbs_patch.GetP(), 2);
    EXPECT_EQ(nurbs_patch.GetQ(), 2);

    // Weights
    Double_Vector wgts(12);
    wgts << 1.000000000000000,   0.853553390593274,   0.853553390593274,  1.000000000000000,   1.000000000000000,
            1.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000,  1.000000000000000,
            1.000000000000000,   1.000000000000000 ;

//    std::cout << "Weights : \n" << nurbs_patch.GetWeights().transpose() << "\n";
//    std::cout << "Control points : \n" << nurbs_patch.GetPoints().transpose() << "\n";
//    std::cout << "Knot U : \n" << nurbs_patch.get_knot_u().transpose() << "\n";
//    std::cout << "Knot V : \n" << nurbs_patch.get_knot_v().transpose() << "\n";
//    std::cout << "Spans : \n" << nurbs_patch.GetKnotRanges() << "\n";
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetWeights(), wgts);
    EXPECT_EQ(nurbs_patch.GetWeight(0), wgts(0));
    EXPECT_EQ(nurbs_patch.GetWeight(1), wgts(1));

    // Indices - Unsigned integer type
    USInt_Matrix indices(nurbs_patch.GetNumberOfElements(), nurbs_patch.GetNumberOfElementPoints());
    indices << 0, 1, 2, 4, 5, 6, 8, 9, 10,
               1, 2, 3, 5, 6, 7, 9, 10, 11;
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetIndices(), indices);

    // KnotRanges
    Double_Matrix  knotranges(nurbs_patch.GetNumberOfElements(), 4);
    knotranges << 0, 0, 0.5, 1,
                  0.5, 0, 1, 1;
    EXPECT_MATRIX_EQUAL(nurbs_patch.GetKnotRanges(), knotranges);

}

