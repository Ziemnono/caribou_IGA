#pragma once

#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Spline.h>

TEST(Spline, Quadratic) {
    using namespace caribou;

    // Shape functions
    {
        using Spline = caribou::geometry::Spline<_3D>;
        using LocalCoordinates = Spline::LocalCoordinates ;
        Spline s;
        // At xi = -1, shape functions [1,0,0]
        EXPECT_EQ(s.L(LocalCoordinates(-1))[0], 1);
        EXPECT_EQ(s.L(LocalCoordinates(-1))[1], 0);
        EXPECT_EQ(s.L(LocalCoordinates(-1))[2], 0);
        // At xi = 0, shape functions [0.25,0.5,0.25]
        EXPECT_EQ(s.L(LocalCoordinates(0))[0], 0.25);
        EXPECT_EQ(s.L(LocalCoordinates(0))[1], 0.5);
        EXPECT_EQ(s.L(LocalCoordinates(0))[2], 0.25);

        EXPECT_EQ(s.L(LocalCoordinates(0.5))[0], 0.0625);
        EXPECT_EQ(s.L(LocalCoordinates(0.5))[1], 0.375);
        EXPECT_EQ(s.L(LocalCoordinates(0.5))[2], 0.5625);

    }
    // 1D
    {
        using Spline = caribou::geometry::Spline<_1D>;
//        using LocalCoordinates  = Spline::LocalCoordinates;
        using WorldCoordinates  = Spline::WorldCoordinates;

        Spline spline(WorldCoordinates(-2.2), WorldCoordinates(1.1), WorldCoordinates(4.4));
        EXPECT_MATRIX_NEAR(spline.node(0), WorldCoordinates(-2.2), 1e-5);
        EXPECT_MATRIX_NEAR(spline.node(1), WorldCoordinates(1.1), 1e-5);
        EXPECT_MATRIX_NEAR(spline.node(2), WorldCoordinates(4.4), 1e-5);

        // Center
        EXPECT_DOUBLE_EQ(spline.center()[0], 1.1);

        // Inverse transformation
        for (const auto & gauss_node : spline.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, spline.local_coordinates(spline.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(spline.node(0)), p1(spline.node(1)), p1(spline.node(2)));
        for (const auto & gauss_node : spline.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(spline.interpolate(x, values), p1(spline.world_coordinates(x)));
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : spline.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = spline.jacobian(x)[0];
            numerical_solution += p1(spline.world_coordinates(x)) * w * detJ;
        }
        FLOATING_POINT_TYPE x0 = spline.node(0)[0];
        FLOATING_POINT_TYPE x2 = spline.node(2)[0];
        auto analytic_solution = static_cast<FLOATING_POINT_TYPE>((5*x2 + x2 * x2) - (5*x0 + x0 * x0)) ;
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);

    }
    // 2D
    {
        using Spline = caribou::geometry::Spline<_2D>;
        using WorldCoordinates  = Spline::WorldCoordinates ;

        WorldCoordinates node_0 {-1.5, -1.5};
        WorldCoordinates node_1 {2, 2};
        WorldCoordinates node_2 {5.5, 5.5};

        Spline spline(node_0, node_1, node_2);

        // Center
        WorldCoordinates center_node = node_0 + (node_2 - node_0).normalized()*(node_2-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(spline.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(spline.center()[1], center_node[1]);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(node_0), p1(node_1), p1(node_2));
        for (const auto & gauss_node : spline.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(spline.interpolate(x, values), p1(spline.world_coordinates(x)));
        }

        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : spline.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  spline.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p1(spline.world_coordinates(x)) * w * detJ;
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
        using Spline = caribou::geometry::Spline<_3D>;
        using WorldCoordinates  = Spline::WorldCoordinates ;

        WorldCoordinates node_0 {-1.5, -1.5, -5.2};
        WorldCoordinates node_1 {2, 2, 24.55};
        WorldCoordinates node_2 {5.5, 5.5, 54.3};

        // Center
        Spline spline(node_0, node_1, node_2);
        WorldCoordinates center_node = node_0 + (node_2 - node_0).normalized()*(node_2-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(spline.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(spline.center()[1], center_node[1]);
        EXPECT_DOUBLE_EQ(spline.center()[2], center_node[2]);

//        // Interpolation
//        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(node_0), p1(node_1), p1(node_2));
//        for (const auto & gauss_node : spline.gauss_nodes()) {
//            const auto x = gauss_node.position;
//            EXPECT_DOUBLE_EQ(spline.interpolate(x, values), p1(spline.world_coordinates(x)));
//        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : spline.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  spline.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p1(spline.world_coordinates(x)) * w * detJ;
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
