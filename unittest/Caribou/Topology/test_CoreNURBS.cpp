#include "topology_test.h"
#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/IO/CoreNURBS.h>
#include <Caribou/Topology/SplinePatch.h>
#include <Caribou/Topology/IO/NURBSReader.h>

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


