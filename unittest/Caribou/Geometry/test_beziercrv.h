#pragma once

#include <Eigen/Dense>
#include <Caribou/constants.h>
#include<Caribou/Geometry/BezierCrv.h>

TEST(BezierCrv, Quadratic) {
    using namespace caribou;

    // Shape functions
    {
        using BezierCrv = caribou::geometry::BezierCrv<_3D>;
        using LocalCoordinates = BezierCrv::LocalCoordinates ;
        BezierCrv b;
        EXPECT_EQ(b.L(LocalCoordinates(-1))[0], 1);
        EXPECT_EQ(b.L(LocalCoordinates(-1))[1], 0);
        EXPECT_EQ(b.L(LocalCoordinates(-1))[2], 0);
        // At xi = 0, shape functions [0.25,0.5,0.25]
        EXPECT_EQ(b.L(LocalCoordinates(0))[0], 0.25);
        EXPECT_EQ(b.L(LocalCoordinates(0))[1], 0.5);
        EXPECT_EQ(b.L(LocalCoordinates(0))[2], 0.25);

        EXPECT_EQ(b.L(LocalCoordinates(1))[0], 0);
        EXPECT_EQ(b.L(LocalCoordinates(1))[1], 0);
        EXPECT_EQ(b.L(LocalCoordinates(1))[2], 1);

    }
    // 1D
    {
        using BezierCrv = caribou::geometry::BezierCrv<_1D>;
//        using LocalCoordinates  = BezierCrv::LocalCoordinates;
        using WorldCoordinates  = BezierCrv::WorldCoordinates;
        // In 1D curve, WorldCoordinates diemsion 2*1, 0 -> Physical coodinate, 1-> Weight
        BezierCrv beziercrv(WorldCoordinates(-2.2), WorldCoordinates(1.1), WorldCoordinates(4.4));
        EXPECT_MATRIX_NEAR(beziercrv.node(0), WorldCoordinates(-2.2), 1e-5);
        EXPECT_MATRIX_NEAR(beziercrv.node(1), WorldCoordinates(1.1), 1e-5);
        EXPECT_MATRIX_NEAR(beziercrv.node(2), WorldCoordinates(4.4), 1e-5);

        // Center
        EXPECT_DOUBLE_EQ(beziercrv.center()[0], 1.1);

        // Inverse transformation
        for (const auto & gauss_node : beziercrv.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, beziercrv.local_coordinates(beziercrv.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(beziercrv.node(0)), p1(beziercrv.node(1)), p1(beziercrv.node(2)));
        for (const auto & gauss_node : beziercrv.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(beziercrv.interpolate(x, values), p1(beziercrv.world_coordinates(x)));
        }


        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : beziercrv.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = beziercrv.jacobian(x).norm();
            numerical_solution += p1(beziercrv.world_coordinates(x)) * w * detJ;
        }
        const auto x0 = beziercrv.node(0)[0]; // a
        const auto w0 = beziercrv.node(0)[1]; // c
        const auto x2 = beziercrv.node(2)[0]; // b
        const auto w2 = beziercrv.node(2)[1]; // d
        const auto d = (beziercrv.node(2) - beziercrv.node(0)).norm();
        FLOATING_POINT_TYPE analytic_solution = d * (x0 + x2 + 3*w0/2. + 3*w2/2. + 5);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);

    }
    // 2D
    {
        using BezierCrv = caribou::geometry::BezierCrv<_2D>;
        using WorldCoordinates  = BezierCrv::WorldCoordinates ;
        // Third dimension is weight
        WorldCoordinates node_0 {-1.5, -1.5};
        WorldCoordinates node_1 {2, 2};
        WorldCoordinates node_2 {5.5, 5.5};

        BezierCrv beziercrv(node_0, node_1, node_2);

        // Center
        WorldCoordinates center_node = node_0 + (node_2 - node_0).normalized()*(node_2-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(beziercrv.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(beziercrv.center()[1], center_node[1]);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(node_0), p1(node_1), p1(node_2));
        for (const auto & gauss_node : beziercrv.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(beziercrv.interpolate(x, values), p1(beziercrv.world_coordinates(x)));
        }

        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : beziercrv.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  beziercrv.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p1(beziercrv.world_coordinates(x)) * w * detJ;
        }

        // Integration
        const auto x0 = node_0[0]; // a
        const auto y0 = node_0[1]; // c
        const auto x2 = node_2[0]; // b
        const auto y2 = node_2[1]; // d
        const auto d = (node_2 - node_0).norm();
        FLOATING_POINT_TYPE analytic_solution = d * (x0 + x2 + 3*y0/2. + 3*y2/2. + 5);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);

    }
    // 3D
    {
        using BezierCrv = caribou::geometry::BezierCrv<_3D>;
        using WorldCoordinates  = BezierCrv::WorldCoordinates ;

        WorldCoordinates node_0 {-1.5, -1.5, -5.2};
        WorldCoordinates node_1 {2, 2, 24.55};
        WorldCoordinates node_2 {5.5, 5.5, 54.3};

        // Center
        BezierCrv beziercrv(node_0, node_1, node_2);
        WorldCoordinates center_node = node_0 + (node_2 - node_0).normalized()*(node_2-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(beziercrv.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(beziercrv.center()[1], center_node[1]);
        EXPECT_DOUBLE_EQ(beziercrv.center()[2], center_node[2]);

//        // Interpolation
//        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(node_0), p1(node_1), p1(node_2));
//        for (const auto & gauss_node : beziercrv.gauss_nodes()) {
//            const auto x = gauss_node.position;
//            EXPECT_DOUBLE_EQ(beziercrv.interpolate(x, values), p1(beziercrv.world_coordinates(x)));
//        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : beziercrv.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  beziercrv.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p1(beziercrv.world_coordinates(x)) * w * detJ;
        }

        const auto x0 = node_0[0];
        const auto y0 = node_0[1];
        const auto z0 = node_0[2];
        const auto x2 = node_2[0];
        const auto y2 = node_2[1];
        const auto z2 = node_2[2];
        const auto d = (node_2 - node_0).norm();
        FLOATING_POINT_TYPE analytic_solution = d * (x0 + x2 + 3*y0/2. + 3*y2/2. + 2*z0 + 2*z2 + 5);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);

    }
}
