#pragma once


#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Quadrature.h>
#include <iostream>

TEST(Quadrature, Order_Three) {
    using namespace caribou;
    {
        using Quadrature = caribou::geometry::Quadrature;
        using Dyn_Vector = geometry::Dyn_Vector;

        int order = 3;
        Quadrature q(order);

        Dyn_Vector qpnts(order);
        Dyn_Vector qwgts(order);

        qpnts(0) = 0.774596669241483;
        qpnts(1) =-0.774596669241483;
        qpnts(2) = 0.000000000000000;

        qwgts(0) = 0.555555555555556;
        qwgts(1) = 0.555555555555556;
        qwgts(2) = 0.888888888888889;

        for (int i = 0; i < order; ++i) {
            EXPECT_DOUBLE_EQ(q.get_point(i), qpnts(i));
            EXPECT_DOUBLE_EQ(q.get_weight(i), qwgts(i));
        }

    }
}

//TEST(NurbsUtils, degree2) {
//    using namespace caribou;
//    {
//        using Nutils = caribou::geometry::NURBS_utils;
//        using Dyn_Vector = geometry::Dyn_Vector;
//        int degree = 2;
//        Nutils N;
//        Dyn_Vector knots(7);
//        knots << 0., 0., 0., 0.5, 1., 1., 1.;

//        Dyn_Vector span(2);
//        span << 0., 0.5;
//    }
//}
