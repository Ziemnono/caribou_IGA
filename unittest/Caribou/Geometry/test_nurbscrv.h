#pragma once


#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/NurbsCrv.h>
#include <iostream>

TEST(NurbsCrv, Basis) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsCrv<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;

        Eigen::Matrix<double, 3, 1> wgts;
        wgts << 1, 1, 1;

        Eigen::Matrix<double, 3, 2> pnts;
        pnts << 0, 0,
                1, 0,
                2, 0;

        Eigen::Matrix<double, 7, 1> knots;
        knots << 0, 0, 0, 0.5, 1, 1, 1;

        Eigen::Matrix<double, 2, 1> span;
        span << 0, 0.5;

        NURBS elem_1(pnts, knots, wgts, span);

        std::cout << "\n ############# IN NURBS curve ###############\n";
        elem_1.print();
        std::cout << "\n at xi = 0 \n";
        std::cout << "Basis functions at " << elem_1.L(LocalCoordinates(0)) << "\n";
        std::cout << "Basis derivative :  \n " << elem_1.dL(LocalCoordinates(0)).transpose() << "\n";
        std::cout << "\n at xi = 1 \n";
        std::cout << "Basis functions at " << elem_1.L(LocalCoordinates(1)) << "\n";
        std::cout << "Basis derivative :  \n " << elem_1.dL(LocalCoordinates(1)).transpose() << "\n";

        std::cout << "\n ############# OUT NURBS curve ###############\n";
    }
}

TEST(NurbsCrv, Center) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsCrv<_2D>;
//        using LocalCoordinates = NURBS::LocalCoordinates;
        using WorldCoordinates  = NURBS::WorldCoordinates ;
    //        double u = 0.25;

        Eigen::Matrix<double, 3, 1> wgts;
        wgts << 1, 1, 1;

        Eigen::Matrix<double, 3, 2> pnts;
        pnts << 0, 0,
                1, 0,
                2, 0;

        Eigen::Matrix<double, 6, 1> knots;
        knots << 0, 0, 0, 1, 1, 1;


        Eigen::Matrix<double, 2, 1> span;
        span << 0, 1.0;


        NURBS n(pnts, knots, wgts, span);

        // Center
        std::cout << "Single element nurbs center : " << n.center() << "\n";
        WorldCoordinates center_node = {1,0};
        EXPECT_DOUBLE_EQ(n.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(n.center()[1], center_node[1]);

    }
}
