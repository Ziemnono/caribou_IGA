#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/SplinePatch.h>
#include <Caribou/Geometry/NurbsCrv.h>
#include <Caribou/Geometry/NurbsSurf.h>

using namespace caribou;
using namespace caribou::geometry;

template<typename dtype>
using Vector = Eigen::Matrix<dtype, Eigen::Dynamic, 1>;

template<typename dtype>
using Matrix = Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Double_Vector = Vector<FLOATING_POINT_TYPE>;
using Int_Vector = Vector<int>;
using Double_Matrix = Matrix<FLOATING_POINT_TYPE>;
using Int_Matrix = Matrix<int>;
using USInt_Matrix = Matrix<UNSIGNED_INTEGER_TYPE>;


TEST (Splinepatch, XY_RECTANGLE) {
    Double_Matrix initial_positions(9,2);
    initial_positions << 0, 0,
            1, 0,
            2, 0,
            0, 1,
            1, 1,
            2, 1,
            0, 2,
            1, 2,
            2, 2;

    Double_Vector weights(9);
    weights << 1,1,1,1,1,1,1,1,1;

    USInt_Matrix indices(1,9);
    indices << 0, 1, 2, 3, 4, 5, 6, 7, 8;

    Double_Matrix knot_ranges(1,4);
    knot_ranges << 0, 0.5, 1, 1;

    Double_Vector knot1(6);
    knot1 << 0, 0, 0, 1, 1, 1;

    Double_Vector knot2(6);
    knot2 <<  0, 0, 0, 1, 1, 1;

    EXPECT_EQ(2,2);

    caribou::topology::SplinePatch<_2D> patch(initial_positions, weights, indices, knot1, knot2, knot_ranges);

    std::cout << " wgts 1, 2 L = \n" << patch.weights({4,5}).transpose() << " ## \n";

    patch.print_spline_patch();

    EXPECT_EQ(patch.number_of_nodes(), 9);
    EXPECT_EQ(patch.canonical_dimension(), 2);
    EXPECT_EQ(patch.number_of_nodes_per_elements(), 9);
    EXPECT_EQ(patch.number_of_elements(), 1);

    EXPECT_EQ(patch.element(0).center()[0], initial_positions(4,0));

    EXPECT_MATRIX_EQUAL(patch.element_indices(0), indices);
    std::cout << "Indices \n " << patch.element_indices(0) << "\n";

    EXPECT_MATRIX_EQUAL(patch.element_knotranges(0), knot_ranges.row(0));
    std::cout << "Knot range \n " << patch.element_knotranges(0) << "\n";

    EXPECT_MATRIX_EQUAL(patch.position(0), initial_positions.row(0));

    Double_Matrix posi(2,2);
    posi << initial_positions.block<2,2>(1,0);
    EXPECT_MATRIX_EQUAL(patch.positions({1,2}), posi );

    EXPECT_EQ(patch.weight(4), 1);

    Double_Vector wgts;
    wgts.resize(2);
    wgts << 1, 1;
    EXPECT_MATRIX_EQUAL(patch.weights({4,5}), wgts );
    std::cout << "Weights \n " << patch.weights({0,1,2,3,4,5,6,7,8}) << "\n";

    std::cout << "all weights " << patch.weights() << "\n";
    std::cout << "all indices " << patch.indices() << "\n";
    std::cout << "all positions " << patch.positions().node(0) << "\n";
    std::cout << "Boundary 1 \n" << patch.element_boundary_nodes(1,0) << "\n";
    std::cout << "Boundary 2 \n" << patch.element_boundary_nodes(2,0) << "\n";
    std::cout << "Boundary 3 \n" << patch.element_boundary_nodes(3,0) << "\n";
    std::cout << "Boundary 4 \n" << patch.element_boundary_nodes(4,0) << "\n";

}

TEST (Splinepatch, bounday_indices) {
    Double_Matrix initial_positions(20,2);
    initial_positions <<
            0, 0,
            1, 0,
            2, 0,
            3, 0,
            4, 0,
            0, 1,
            1, 1,
            2, 1,
            3, 1,
            4, 1,
            0, 2,
            1, 2,
            2, 2,
            3, 2,
            4, 2,
            0, 3,
            1, 3,
            2, 3,
            3, 3,
            4, 3;

    Double_Vector weights(20);
    weights << 1,1,1,1,1,
               1,1,1,1,1,
               1,1,1,1,1,
               1,1,1,1,1;

    USInt_Matrix indices(6,9);
    indices << 0,1,2,5,6,7,10,11,12,
               1,2,3,6,7,8,11,12,13,
               2,3,4,7,8,9,12,13,14,
               5,6,7,10,11,12,15,16,17,
               6,7,8,11,12,13,16,17,18,
               7,8,9,12,13,14,17,18,19;

    Double_Matrix knot_ranges(6,4);
    knot_ranges << 0,      0, 0.33, 0.5,
                   0.33,   0, 0.66, 0.5,
                   0.66,   0,    1, 0.5,
                      0, 0.5, 0.33, 1,
                   0.33, 0.5, 0.66, 1,
                   0.66, 0.5,    1, 1;

    Double_Vector knot1(8);
    knot1 << 0, 0, 0, 0.33, 0.66, 1, 1, 1;

    Double_Vector knot2(7);
    knot2 <<  0, 0, 0, 0.5, 1, 1, 1;

    caribou::topology::SplinePatch<_2D> patch(initial_positions, weights, indices, knot1, knot2, knot_ranges);


    patch.print_spline_patch();

    EXPECT_EQ(patch.number_of_nodes(), 20);
    EXPECT_EQ(patch.canonical_dimension(), 2);
    EXPECT_EQ(patch.number_of_nodes_per_elements(), 9);
    EXPECT_EQ(patch.number_of_elements(), 6);


    std::cout << "all weights " << patch.weights() << "\n";
    std::cout << "all indices " << patch.indices() << "\n";
    std::cout << "all positions " << patch.positions().node(0) << "\n";
    std::cout << "Boundary 1 \n" << patch.element_boundary_nodes(1,0) << "\n";
    std::cout << "Boundary 2 \n" << patch.element_boundary_nodes(2,0) << "\n";
    std::cout << "Boundary 3 \n" << patch.element_boundary_nodes(3,0) << "\n";
    std::cout << "Boundary 3 : 2 \n" << patch.element_boundary_nodes(3,2) << "\n";
    std::cout << "Boundary 4 \n" << patch.element_boundary_nodes(4,0) << "\n";


    auto crv = patch.boundary_element(1,0);
    std::cout << "\nprinting the bottom crv\n";
    crv.print();

}
