#include "topology_test.h"
#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/IO/CoreNURBS.h>
#include <Caribou/Topology/IO/CoreNURBS_1.h>
#include <Caribou/Topology/SplinePatch.h>
#include <Caribou/Topology/IO/NURBSReader.h>
#include <Caribou/Topology/IO/trial_io.h>

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
    USInt_Matrix element_connectivity(2, 9);
    element_connectivity << 0, 1, 2, 4, 5, 6, 8, 9, 10,
                            1, 2, 3, 5, 6, 7, 9, 10, 11;

    para_topo topo_info(element_ranges, element_connectivity);
    EXPECT_MATRIX_EQUAL(topo_info.get_elrange(), element_ranges);
    EXPECT_MATRIX_EQUAL(topo_info.get_elconn(), element_connectivity);
}

TEST(CoreNURBS, reader) {
    using namespace caribou;
    using namespace caribou::topology::io;
    coreNurbs nurbs_patch;
    nurbs_patch.SetFileName(executable_directory_path + "/meshes/splines/3D_spline_surface.txt");
    nurbs_patch.Update();

    EXPECT_EQ(nurbs_patch.GetP(), 2);
    EXPECT_EQ(nurbs_patch.GetQ(), 2);
    EXPECT_EQ(nurbs_patch.get_no_pnts_u(), 4);
    EXPECT_EQ(nurbs_patch.get_no_pnts_v(), 3);
    EXPECT_EQ(nurbs_patch.GetNumberOfPoints(), 12);
    EXPECT_EQ(nurbs_patch.get_no_elems_u(), 2);
    EXPECT_EQ(nurbs_patch.get_no_elems_v(), 1);

    Double_Vector knot_u(7), knot_v(6);
    knot_u << 0, 0, 0, 0.5, 1, 1, 1;
    knot_v << 0, 0, 0, 1, 1, 1;

    EXPECT_MATRIX_EQUAL(nurbs_patch.get_knot_u(), knot_u);
    EXPECT_MATRIX_EQUAL(nurbs_patch.get_knot_v(), knot_v);


}

//TEST(CoreNURBS, tt1) {
//    using namespace caribou;
//    using namespace caribou::topology::io;
//    trial t;
//    t.SetFileName(executable_directory_path + "/meshes/splines/test_io.txt");
//    t.Update();
//    EXPECT_EQ(2, 2);
//    EXPECT_EQ(t.GetPdim(), 10);
//    EXPECT_EQ(t.GetRealDim(), 11);
//    Int_Vector b(2);
//    b << 10, 11;
//    EXPECT_MATRIX_EQUAL(t.get_vec(), b);
//    std::cout << "Vector " << b << " \n";
//}


//TEST(CoreNURBS, reader_1) {
//    using namespace caribou;
//    using namespace caribou::topology::io;
//    coreNurbs_1 nurbs_patch;
//    nurbs_patch.SetFileName(executable_directory_path + "/meshes/splines/3D_spline_surface.txt");
//    nurbs_patch.Update();

//    EXPECT_EQ(nurbs_patch.GetP(), 2);
//    EXPECT_EQ(nurbs_patch.GetQ(), 2);
//    EXPECT_EQ(nurbs_patch.get_no_pnts_u(), 4);
//    EXPECT_EQ(nurbs_patch.get_no_pnts_v(), 3);
//    EXPECT_EQ(nurbs_patch.GetNumberOfPoints(), 12);
//    Double_Vector knot_u(7), knot_v(6);
//    knot_u << 0, 0, 0, 0.5, 1, 1, 1;
//    knot_v << 0, 0, 0, 1, 1, 1;

//    EXPECT_MATRIX_EQUAL(nurbs_patch.get_knot_u(), knot_u);
//    EXPECT_MATRIX_EQUAL(nurbs_patch.get_knot_v(), knot_v);
//}
