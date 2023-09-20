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

TEST(Quadrature, Two_dimensional) {
    using namespace caribou;
    {
        using Quadrature = caribou::geometry::Quadrature;

        const int gauss_u = 3; // Gauss points in U direction
        const int gauss_v = 3; // Gauss points in V direction

        std::cout << "Gauss nodes test is runnig\n";
        Quadrature quad_u(gauss_u);
        Quadrature quad_v(gauss_v);
        int count = 0;
        for (int j = 0; j < gauss_v; ++j) {
            for (int i = 0; i < gauss_u; ++i) {

                if ((count == 5) || (count == 7)){
                    std::cout << "============== " << count << " ===============\n";
                    std::cout << "q x : " << quad_v.get_point(j) << "\n";
                    std::cout << "q y : " << quad_u.get_point(i) << "\n";
                    std::cout << "w   : " << quad_v.get_weight(i) * quad_u.get_weight(j) << "\n";
                }
                count++;
            }
        }

    }
}

