#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/SplinePatch.h>
#include <Caribou/Geometry/BezierSurf.h>

using namespace caribou;
using namespace caribou::geometry;

template<typename dtype>
using Vector = Eigen::Matrix<dtype, 1, Eigen::Dynamic, Eigen::RowMajor>;

template<typename dtype>
using Matrix = Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Double_Vector = Vector<double>;
using Int_Vector = Vector<int>;
using Double_Matrix = Matrix<double>;
using Int_Matrix = Matrix<int>;
using USInt_Matrix = Matrix<UNSIGNED_INTEGER_TYPE>;


TEST (Splinepatch, XY_RECTANGLE) {
    Double_Matrix initial_positions(9,2);
    initial_positions << 0,0,
            1,0,
            2,0,
            0,1,
            1,1,
            2,1,
            0,2,
            1,2,
            2,2;

    Double_Vector weights(9);
    weights << 1, 1, 1, 1, 0.75, 0.7171, 1, 1, 1;

    USInt_Matrix indices(1,9);
    indices << 0, 1, 2, 3, 4, 5, 6, 7, 8;

    Double_Matrix knot_ranges(1,4);
    knot_ranges << 0, 0, 1, 1;

    Double_Matrix mat(9,9);
    mat <<  1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1;
    std::vector<Double_Matrix> ext_mat  ={ mat };


    caribou::topology::SplinePatch<_2D> patch(initial_positions, weights, indices, knot_ranges, ext_mat);

//    std::cout << " wgts 1, 2 L = \n" << patch.weights({4,5}) << " ## \n";

    EXPECT_EQ(patch.number_of_nodes(), 9);
    EXPECT_EQ(patch.canonical_dimension(), 2);
    EXPECT_EQ(patch.number_of_nodes_per_elements(), 9);
    EXPECT_EQ(patch.number_of_elements(), 1);

    EXPECT_EQ(patch.element(0).center()[0], initial_positions(4,0));

    EXPECT_MATRIX_EQUAL(patch.element_indices(0), indices);

    EXPECT_MATRIX_EQUAL(patch.element_knotranges(0), knot_ranges.row(0));

    EXPECT_MATRIX_EQUAL(patch.element_extraction(0), ext_mat[0]);

    EXPECT_MATRIX_EQUAL(patch.position(0), initial_positions.row(0));

    Double_Matrix posi(2,2);
    posi << initial_positions.block<2,2>(1,0);
    EXPECT_MATRIX_EQUAL(patch.positions({1,2}), posi );


    EXPECT_EQ(patch.weight(4), 0.75);

    Double_Vector wgts(2);
    wgts << 0.75, 0.7171;
    EXPECT_MATRIX_EQUAL(patch.weights({4,5}), wgts );

}

