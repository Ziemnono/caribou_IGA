#pragma once

#include <Eigen/Dense>
#include <Caribou/constants.h>
#include<Caribou/Geometry/BezierSurf.h>

TEST(BezierSurf, Quadratic) {
    using namespace caribou;

    // Shape functions
    {
        using BezierSurf = caribou::geometry::BezierSurf<_3D>;
        using LocalCoordinates = BezierSurf::LocalCoordinates ;
        BezierSurf b;
        EXPECT_EQ(b.L(LocalCoordinates(-1,-1))[0], 1);
        EXPECT_EQ(b.L(LocalCoordinates(-1,-1))[8], 0);

        EXPECT_EQ(b.L(LocalCoordinates(1, 1))[8], 1);

    }
    // 2D
    {
        using BezierSurf = caribou::geometry::BezierSurf<_2D>;
        using WorldCoordinates  = BezierSurf::WorldCoordinates ;

        WorldCoordinates node_0 {0,0};
        WorldCoordinates node_1 {1,0};
        WorldCoordinates node_2 {2,0};
        WorldCoordinates node_3 {0,1};
        WorldCoordinates node_4 {1,1};
        WorldCoordinates node_5 {2,1};
        WorldCoordinates node_6 {0,2};
        WorldCoordinates node_7 {1,2};
        WorldCoordinates node_8 {2,2};
        BezierSurf beziersurf(node_0, node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8);

        // Center
        WorldCoordinates center_node = {1,1};
        EXPECT_DOUBLE_EQ(beziersurf.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(beziersurf.center()[1], center_node[1]);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 9, 1> values;
        values << p1(node_0), p1(node_1), p1(node_2), p1(node_3), p1(node_4),
                p1(node_5), p1(node_6), p1(node_7), p1(node_8);
//        std::cout << "Values\n " << values << "\n" ;
        for (const auto & gauss_node : beziersurf.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(beziersurf.interpolate(x, values), p1(beziersurf.world_coordinates(x)));
        }

        FLOATING_POINT_TYPE numerical_solution = 0;
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 2> JJ;
        for (const auto & gauss_node : beziersurf.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  beziersurf.jacobian(x);
            JJ = beziersurf.nodes().transpose() *beziersurf.dL(x);
            const auto detJ = J.norm();
            numerical_solution += p1(beziersurf.world_coordinates(x)) * w * detJ;
        }

//        // Integration
//        const auto x0 = node_0[0];
//        const auto y0 = node_0[1];
//        const auto w0 = node_0[2];
//        const auto x2 = node_2[0];
//        const auto y2 = node_2[1];
//        const auto w2 = node_2[2];
//        const auto d = (node_2 - node_0).norm();
//        FLOATING_POINT_TYPE analytic_solution = d * (x0 + x2 + 3*y0/2. + 3*y2/2. + 2*w0 + 2*w2 + 5);
//        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);

    }
    // 3D
    {
        using BezierSurf = caribou::geometry::BezierSurf<_3D>;
        using WorldCoordinates  = BezierSurf::WorldCoordinates ;
        WorldCoordinates node_0 {0,0,0};
        WorldCoordinates node_1 {1,0,0};
        WorldCoordinates node_2 {2,0,0};
        WorldCoordinates node_3 {0,1,0};
        WorldCoordinates node_4 {1,1,0};
        WorldCoordinates node_5 {2,1,0};
        WorldCoordinates node_6 {0,2,0};
        WorldCoordinates node_7 {1,2,0};
        WorldCoordinates node_8 {2,2,0};
        BezierSurf beziersurf(node_0, node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8);

        // Center
        WorldCoordinates center_node = {1,1,0};
        EXPECT_DOUBLE_EQ(beziersurf.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(beziersurf.center()[1], center_node[1]);
        EXPECT_DOUBLE_EQ(beziersurf.center()[2], center_node[2]);


        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 9, 1> values;
        values << p1(node_0), p1(node_1), p1(node_2), p1(node_3), p1(node_4),
                p1(node_5), p1(node_6), p1(node_7), p1(node_8);
    //        std::cout << "Values\n " << values << "\n" ;
        for (const auto & gauss_node : beziersurf.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(beziersurf.interpolate(x, values), p1(beziersurf.world_coordinates(x)));
        }


//        // Integration
//        FLOATING_POINT_TYPE numerical_solution = 0;
//        for (const auto & gauss_node : beziersurf.gauss_nodes()) {
//            const auto x = gauss_node.position;
//            const auto w = gauss_node.weight;
//            const auto J =  beziersurf.jacobian(x);
//            const auto detJ = J.norm();
//            numerical_solution += p1(beziersurf.world_coordinates(x)) * w * detJ;
//        }

//        std::cout << "Numerical sol \n " << numerical_solution << "\n";

//        const auto x0 = node_0[0];
//        const auto y0 = node_0[1];
//        const auto z0 = node_0[2];
//        const auto w0 = node_0[3];
//        const auto x2 = node_2[0];
//        const auto y2 = node_2[1];
//        const auto z2 = node_2[2];
//        const auto w2 = node_2[3];
//        const auto d = (node_2 - node_0).norm();
//        FLOATING_POINT_TYPE analytic_solution = d * (x0 + x2 + 3*y0/2. + 3*y2/2. + 2*z0 + 2*z2 + 5*w0/2 + 5*w2/2 + 5);
//        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);

    }
}
