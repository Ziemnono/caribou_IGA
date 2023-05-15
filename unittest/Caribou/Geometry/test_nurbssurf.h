#pragma once


#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/NurbsSurf.h>
#include <iostream>

TEST(NurbsSurf, Basis) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsSurf<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;

        Eigen::Matrix<double, 9, 1> wgts;
        wgts << 1, 1, 1, 1, 1, 1, 1, 1, 1;

        Eigen::Matrix<double, 9, 2> pnts;
        pnts << 0, 0,
                1, 0,
                2, 0,
                0, 1,
                1, 1,
                2, 1,
                0, 2,
                1, 2,
                2, 2;

        Eigen::Matrix<double, 7, 1> knot1;
        knot1 << 0, 0, 0, 0.5, 1, 1, 1;

        Eigen::Matrix<double, 7, 1> knot2;
        knot2 << 0, 0, 0, 0.5, 1, 1, 1;

        Eigen::Matrix<double, 4, 1> span;
        span << 0, 0, 0.5, 0.5;

        NURBS elem_1(pnts, knot1, knot2, wgts, span);
        elem_1.print();

        std::cout << "Basis functions at " << elem_1.L(LocalCoordinates(0,1)) << "\n";
        std::cout << "Basis derivative :  \n " << elem_1.dL(LocalCoordinates(0,1)).transpose() << "\n";
    }
}

TEST(NurbsSurf, Center) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsSurf<_2D>;
//        using LocalCoordinates = NURBS::LocalCoordinates;
        using WorldCoordinates  = NURBS::WorldCoordinates ;
    //        double u = 0.25;

        Eigen::Matrix<double, 9, 1> wgts;
        wgts << 1, 1, 1, 1, 1, 1, 1, 1, 1;

        Eigen::Matrix<double, 9, 2> pnts;
        pnts << 0, 0,
                1, 0,
                2, 0,
                0, 1,
                1, 1,
                2, 1,
                0, 2,
                1, 2,
                2, 2;

        Eigen::Matrix<double, 6, 1> knot1;
        knot1 << 0, 0, 0, 1, 1, 1;

        Eigen::Matrix<double, 6, 1> knot2;
        knot2 << 0, 0, 0, 1, 1, 1;

        Eigen::Matrix<double, 4, 1> span;
        span << 0, 0, 1.0, 1.0;


        NURBS n(pnts, knot1, knot2, wgts, span);

        n.print();


        // Center
        std::cout << "Single element nurbs center : " << n.center() << "\n";
        WorldCoordinates center_node = {1,1};
        EXPECT_DOUBLE_EQ(n.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(n.center()[1], center_node[1]);

    }
}
