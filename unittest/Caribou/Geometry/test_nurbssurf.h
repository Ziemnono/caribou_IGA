#pragma once


#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/NurbsSurf.h>
#include <iostream>

//TEST(NurbsSurf, Basis) {
//    using namespace caribou;
//    {
//        using NURBS = caribou::geometry::NurbsSurf<_2D>;
//        using LocalCoordinates = NURBS::LocalCoordinates;

//        Eigen::Matrix<UNSIGNED_INTEGER_TYPE, 2,1> degrees;
//        degrees << 2,2;

//        Eigen::Matrix<double, 9, 1> wgts;
//        wgts << 1, 1, 1, 1, 1, 1, 1, 1, 1;

//        Eigen::Matrix<double, 9, 2> pnts;
//        pnts << 0, 0,
//                1, 0,
//                2, 0,
//                0, 1,
//                1, 1,
//                2, 1,
//                0, 2,
//                1, 2,
//                2, 2;

//        Eigen::Matrix<double, 7, 1> knot1;
//        knot1 << 0, 0, 0, 0.5, 1, 1, 1;

//        Eigen::Matrix<double, 7, 1> knot2;
//        knot2 << 0, 0, 0, 0.5, 1, 1, 1;

//        Eigen::Matrix<double, 4, 1> span;
//        span << 0, 0, 0.5, 0.5;

//        NURBS elem_1(degrees, pnts, knot1, knot2, wgts, span);

//    }
//}

TEST(NurbsSurf, Linear) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsSurf<_2D>;
//        using LocalCoordinates = NURBS::LocalCoordinates;

        Eigen::Matrix<UNSIGNED_INTEGER_TYPE, 2,1> degrees;
        degrees << 1,1;

        Eigen::Matrix<double, 4, 1> wgts;
        wgts << 1, 1, 1, 1;

        Eigen::Matrix<double, 4, 2> pnts;
        pnts << 0, 0,
                2, 0,
                0, 2,
                2, 2;

        Eigen::Matrix<double, 4, 1> knot1;
        knot1 << 0, 0, 1, 1;

        Eigen::Matrix<double, 4, 1> knot2;
        knot2 << 0, 0, 1, 1;

        Eigen::Matrix<double, 4, 1> span;
        span << 0, 0, 1, 1;

        NURBS elem_1(degrees, pnts, knot1, knot2, wgts, span);

    }
}

TEST(NurbsSurf, LinearQuadratic) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsSurf<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;

        Eigen::Matrix<UNSIGNED_INTEGER_TYPE, 2,1> degrees;
        degrees << 2, 1;

        Eigen::Matrix<double, 6, 1> wgts;
        wgts << 1, 1, 1, 1, 1, 1;

        Eigen::Matrix<double, 6, 2> pnts;
        pnts << 0, 0,
                1, 0,
                2, 0,
                0, 1,
                1, 1,
                2, 1;

        Eigen::Matrix<double, 6, 1> knot1;
        knot1 << 0, 0, 0, 1, 1, 1;

        Eigen::Matrix<double, 4, 1> knot2;
        knot2 << 0, 0, 1, 1;

        Eigen::Matrix<double, 4, 1> span;
        span << 0, 0, 1, 1;

        NURBS elem_1(degrees, pnts, knot1, knot2, wgts, span);
//        elem_1.print();

        LocalCoordinates l{0,0};
        std::cout << "At 0, 0  --> " << elem_1.world_coordinates(l) << "\n";

        LocalCoordinates l2{0.25,0.25};
        std::cout << "At 0.25, 0.25  --> " << elem_1.world_coordinates(l2) << "\n";

    }
}

TEST(NurbsSurf, Quadratic_1) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsSurf<_2D>;

        Eigen::Matrix<UNSIGNED_INTEGER_TYPE, 2,1> degrees;
        degrees << 2,2;

        Eigen::Matrix<double, 9, 1> wgts;
        wgts << 1.0000, 0.8536, 0.8536, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000;

        Eigen::Matrix<double, 9, 2> pnts;
        pnts << -1.0000, 0,
                -1.0000, 0.4142,
                -0.4142, 1.0000,
                -2.5000, 0,
                -2.5000, 0.7500,
                -0.7500, 2.5000,
                -4.0000, 0,
                -4.0000, 4.0000,
                -4.0000, 4.0000;

        Eigen::Matrix<double, 7, 1> knot1;
        knot1 << 0, 0, 0, 0.5, 1, 1, 1;
        Eigen::Matrix<double, 6, 1> knot2;
        knot2 << 0, 0, 0, 1, 1, 1;
        Eigen::Matrix<double, 4, 1> span;
        span << 0, 0.0, 0.5, 1.0;
        NURBS n(degrees, pnts, knot1, knot2, wgts, span);
        Eigen::Matrix<double, 2,2> J;
        J.setZero();
        for (const auto & gauss_node : n.gauss_nodes()) {
            const auto x = gauss_node.position;
            J = J + n.jacobian(x);
        }

        Eigen::Matrix <double, 1, 2> center;
        Eigen::Matrix <double, 2, 2> Jacobian;

        center << -2.4138, 1.2571;
        Jacobian << 6.0756,  -28.0248,
                    38.8534,  20.1949;
        Eigen::Matrix <double, 1, 2> c = n.center();
        EXPECT_MATRIX_NEAR(c, center, 1e-3);
        EXPECT_MATRIX_NEAR(J, Jacobian, 1e-3);

    }
}

TEST(NurbsSurf, Quadratic_2) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsSurf<_2D>;

        Eigen::Matrix<UNSIGNED_INTEGER_TYPE, 2,1> degrees;
        degrees << 2,2;

        Eigen::Matrix<double, 9, 1> wgts;
        wgts << 0.8536, 0.8536, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000;

        Eigen::Matrix<double, 9, 2> pnts;
        pnts << -1.0000, 0.4142,
                -0.4142, 1.0000,
                0, 1.0000,
                -2.5000, 0.7500,
                -0.7500, 2.5000,
                0, 2.5000,
                -4.0000, 4.0000,
                -4.0000, 4.0000,
                0, 4.0000;

        Eigen::Matrix<double, 7, 1> knot1;
        knot1 << 0, 0, 0, 0.5, 1, 1, 1;
        Eigen::Matrix<double, 6, 1> knot2;
        knot2 << 0, 0, 0, 1, 1, 1;
        Eigen::Matrix<double, 4, 1> span;
        span << 0.5, 0.0, 1., 1.0;
        NURBS n(degrees, pnts, knot1, knot2, wgts, span);

        Eigen::Matrix<double, 2,2> J;
        J.setZero();
        for (const auto & gauss_node : n.gauss_nodes()) {
            const auto x = gauss_node.position;
            J = J + n.jacobian(x);
        }

        Eigen::Matrix <double, 1, 2> center;
        Eigen::Matrix <double, 2, 2> Jacobian;

        center << -1.2571, 2.4138;
        Jacobian << 38.8534,  -20.1949,
                     6.0756,   28.0248;

        Eigen::Matrix <double, 1, 2> c = n.center();
        EXPECT_MATRIX_NEAR(c, center, 1e-3);
        EXPECT_MATRIX_NEAR(J, Jacobian, 1e-3);
    }
}
