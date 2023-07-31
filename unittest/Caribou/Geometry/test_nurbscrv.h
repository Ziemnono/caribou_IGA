#pragma once
#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/NurbsCrv.h>
#include <iostream>

TEST(NurbsCrv, Degree4) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsCrv<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;
        using Scalar = FLOATING_POINT_TYPE;
        using USINT = UNSIGNED_INTEGER_TYPE;

        USINT degree = 4;
        Eigen::Matrix<Scalar, 5, 1> wgts;
        wgts << 1, 1, 1, 1, 1;
        Eigen::Matrix<Scalar, 5, 2> pnts;
        pnts << 0, 0,
                1, 0,
                2, 0,
                3, 0,
                4, 0;
        Eigen::Matrix<Scalar, 10, 1> knots;
        knots << 0, 0, 0, 0, 0, 1, 1, 1, 1, 1;

        Eigen::Matrix<Scalar, 2, 1> span;
        span << 0, 1;

        NURBS elem(degree, pnts, knots, wgts, span);

        std::vector<Scalar> q = {-0.25, 0, 0.5};
        Eigen::Matrix<Scalar, 3, 5> basis;
        basis << 0.1526,    0.3662,   0.3296,    0.1318,    0.0198,
                0.0625,    0.2500,    0.3750,    0.2500,    0.0625,
                0.0039,    0.0469,    0.2109,    0.4219,    0.3164;

        Eigen::Matrix<Scalar, 3, 5> deris;
        deris <<-0.9766,   -0.7813,    0.7031,    0.8438,    0.2109,
                -0.5000,   -1.0000,         0,    1.0000,    0.5000,
                -0.0625,   -0.5000,   -1.1250,         0,    1.6875;

        Eigen::Matrix<Scalar, 2, 1> jacob = {4.,0.};

        for (int i = 0; i < 3; i++) {
            Eigen::Matrix<Scalar, 5, 1> b = basis.row(i); // basis
            Eigen::Matrix<Scalar, 5, 1> cb = elem.L(LocalCoordinates(q[i])); // Caribou basis

            Eigen::Matrix<Scalar, 5, 1> db = deris.row(i); // basis
            Eigen::Matrix<Scalar, 5, 1> cdb = elem.dL(LocalCoordinates(q[i])); // Caribou basis

            EXPECT_MATRIX_NEAR(cb, b, 1e-3);
            EXPECT_MATRIX_NEAR(cdb, db, 1e-3);
            EXPECT_MATRIX_NEAR(elem.jacobian(LocalCoordinates(q[i])), jacob, 1e-3);
        }
    }
}

TEST(NurbsCrv, Degree3Elem1) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsCrv<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;
        using Scalar = FLOATING_POINT_TYPE;
        using USINT = UNSIGNED_INTEGER_TYPE;

        USINT degree = 3;
        Eigen::Matrix<Scalar, 4, 1> wgts;
        wgts << 1, 1, 1, 1;
        Eigen::Matrix<Scalar, 4, 2> pnts;
        pnts << 0, 0,
                1, 1,
                2, 1,
                3, 0;

        Eigen::Matrix<Scalar, 9, 1> knots;
        knots << 0, 0, 0, 0, 0.5, 1, 1, 1, 1;

        Eigen::Matrix<Scalar, 2, 1> span;
        span << 0, 0.5;

        NURBS elem(degree, pnts, knots, wgts, span);

        std::vector<Scalar> q = {-0.25, 0, 0.5};
        Eigen::Matrix<Scalar, 3, 4> basis;
        basis <<0.2441,    0.5845,    0.1582,    0.0132,
                0.1250,    0.5938,    0.2500,    0.0313,
                0.0156,    0.4570,    0.4219,    0.1055;

        Eigen::Matrix<Scalar, 3, 4> deris;
        deris <<-2.3438,    0.7266,    1.4063,    0.2109,
                -1.5000,   -0.3750,    1.5000,    0.3750,
                -0.3750,   -1.5938,    1.1250,    0.8438;

        Eigen::Matrix<Scalar, 3, 2> jacob;
        jacob << 4.1719,    2.1328,
                 3.7500,    1.1250,
                 3.1875,   -0.4688;

        for (int i = 0; i < 3; i++) {
            Eigen::Matrix<Scalar, 4, 1> b = basis.row(i); // basis
            Eigen::Matrix<Scalar, 4, 1> cb = elem.L(LocalCoordinates(q[i])); // Caribou basis

            Eigen::Matrix<Scalar, 4, 1> db = deris.row(i); // basis
            Eigen::Matrix<Scalar, 4, 1> cdb = elem.dL(LocalCoordinates(q[i])); // Caribou basis

            EXPECT_MATRIX_NEAR(cb, b, 1e-3);
            EXPECT_MATRIX_NEAR(cdb, db, 1e-3);
            Eigen::Matrix<Scalar, 2, 1> j = jacob.row(i);
            EXPECT_MATRIX_NEAR(elem.jacobian(LocalCoordinates(q[i])), j, 1e-3);
        }
    }
}

TEST(NurbsCrv, Degree3Elem2) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsCrv<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;
        using Scalar = FLOATING_POINT_TYPE;
        using USINT = UNSIGNED_INTEGER_TYPE;

        USINT degree = 3;
        Eigen::Matrix<Scalar, 4, 1> wgts;
        wgts << 1, 1, 1, 1;
        Eigen::Matrix<Scalar, 4, 2> pnts;
        pnts << 1, 1,
                2, 1,
                3, 0,
                4, 1;

        Eigen::Matrix<Scalar, 9, 1> knots;
        knots << 0, 0, 0, 0, 0.5, 1, 1, 1, 1;

        Eigen::Matrix<Scalar, 2, 1> span;
        span << 0.5, 1.;

        NURBS elem(degree, pnts, knots, wgts, span);

//        elem1.print();
        std::vector<Scalar> q = {-0.25, 0, 0.5};
        Eigen::Matrix<Scalar, 3, 4> basis;
        basis <<0.0610,    0.3418,    0.5444,    0.0527,
                0.0313,    0.2500,    0.5938,    0.1250,
                0.0039,    0.0781,    0.4961,    0.4219;

        Eigen::Matrix<Scalar, 3, 4> deris;
        deris <<-0.5859,   -1.4063,    1.1484,    0.8438,
                -0.3750,   -1.5000,    0.3750,    1.5000,
                -0.0938,   -1.1250,   -2.1563,    3.3750;

        Eigen::Matrix<Scalar, 3, 2> jacob;
        jacob << 3.4219,   -1.1484,
                 3.7500,   -0.3750,
                 4.6875,    2.1563;

        for (int i = 0; i < 3; i++) {
            Eigen::Matrix<Scalar, 4, 1> b = basis.row(i); // basis
            Eigen::Matrix<Scalar, 4, 1> cb = elem.L(LocalCoordinates(q[i])); // Caribou basis

            Eigen::Matrix<Scalar, 4, 1> db = deris.row(i); // basis
            Eigen::Matrix<Scalar, 4, 1> cdb = elem.dL(LocalCoordinates(q[i])); // Caribou basis

            EXPECT_MATRIX_NEAR(cb, b, 1e-3);
            EXPECT_MATRIX_NEAR(cdb, db, 1e-3);
            Eigen::Matrix<Scalar, 2, 1> j = jacob.row(i);
            EXPECT_MATRIX_NEAR(elem.jacobian(LocalCoordinates(q[i])), j, 1e-3);
        }
    }
}

TEST(NurbsCrv, QuardterCircle_1) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsCrv<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;
        using Scalar = FLOATING_POINT_TYPE;
        using USINT = UNSIGNED_INTEGER_TYPE;

        USINT degree = 2;
        Eigen::Matrix<Scalar, 3, 1> wgts;
        wgts << 1, 0.8536, 0.8536;
        Eigen::Matrix<Scalar, 3, 2> pnts;
        pnts << 1., 0,
                1., 0.4142,
                0.4142, 1.;

        Eigen::Matrix<Scalar, 7, 1> knots;
        knots << 0, 0, 0, 0.5, 1, 1, 1;

        Eigen::Matrix<Scalar, 2, 1> span;
        span << 0., 0.5;

        NURBS elem(degree, pnts, knots, wgts, span);

        Eigen::Matrix<Scalar, 3,2> nn = elem.nodes();
        EXPECT_MATRIX_NEAR(nn, pnts, 1e-3);
//        elem.print();
        std::vector<Scalar> q = {-0.25, 0, 0.5};
        Eigen::Matrix<Scalar, 3, 3> basis;
        basis <<   0.428887089469278,   0.505215267007946,   0.065897643522776,
                   0.280835767243316,   0.599303527297237,   0.119860705459447,
                   0.072442770211533,   0.649290060851927,   0.278267168936540;

        Eigen::Matrix<Scalar, 3, 3> deris;
        deris <<  -2.572529075392716,   1.843139906646620,   0.729389168746097,
                  -2.154315083521857,   1.156005473817829,   0.998309609704028,
                  -1.146791513695777,  -0.384519194542194,   1.531310708237971;

        Eigen::Matrix<Scalar, 3, 2> jacob;
        jacob <<   -0.427266282803425,   1.492842715430207,
                   -0.584796229917229,   1.477142755136909,
                   -0.897021044678654,   1.372037642865816;

        for (int i = 0; i < 3; i++) {
            Eigen::Matrix<Scalar, 3, 1> b = basis.row(i); // basis
            Eigen::Matrix<Scalar, 3, 1> cb = elem.L(LocalCoordinates(q[i])); // Caribou basis

            Eigen::Matrix<Scalar, 3, 1> db = deris.row(i); // basis
            Eigen::Matrix<Scalar, 3, 1> cdb = elem.dL(LocalCoordinates(q[i])); // Caribou basis

            EXPECT_MATRIX_NEAR(cb, b, 1e-3);
            EXPECT_MATRIX_NEAR(cdb, db, 1e-3);
            Eigen::Matrix<Scalar, 2, 1> j = jacob.row(i);
            EXPECT_MATRIX_NEAR(elem.jacobian(LocalCoordinates(q[i])), j, 1e-3);

        }
    }
}

TEST(NurbsCrv, QuardterCircle_2) {
    using namespace caribou;
    {
        using NURBS = caribou::geometry::NurbsCrv<_2D>;
        using LocalCoordinates = NURBS::LocalCoordinates;
        using Scalar = FLOATING_POINT_TYPE;
        using USINT = UNSIGNED_INTEGER_TYPE;

        USINT degree = 2;
        Eigen::Matrix<Scalar, 3, 1> wgts;
        wgts << 0.8536, 0.8536, 1;
        Eigen::Matrix<Scalar, 3, 2> pnts;
        pnts << 1., 0.4142,
                0.4142, 1.,
                0., 1;

        Eigen::Matrix<Scalar, 7, 1> knots;
        knots << 0, 0, 0, 0.5, 1, 1, 1;

        Eigen::Matrix<Scalar, 2, 1> span;
        span << 0.5, 1.;

        NURBS elem(degree, pnts, knots, wgts, span);

        Eigen::Matrix<Scalar, 3,2> nn = elem.nodes();
        EXPECT_MATRIX_NEAR(nn, pnts, 1e-3);
//        elem.print();
        std::vector<Scalar> q = {-0.25, 0, 0.5};
        Eigen::Matrix<Scalar, 3, 3> basis;
        basis <<0.190712804747265,   0.648423536140702,   0.160863659112033,
                0.119860705459447,   0.599303527297237,   0.280835767243316,
                0.028500454084086,   0.370505903093114,   0.600993642822800;

        Eigen::Matrix<Scalar, 3, 3> deris;
        deris <<  -1.268469901391292,  -0.406999423506397,   1.675469324897688,
                  -0.998309609704028,  -1.156005473817829,   2.154315083521857,
                  -0.469381269761986,  -2.453898384142845,   2.923279653904830;

        Eigen::Matrix<Scalar, 3, 2> jacob;
        jacob <<   -1.437054582485672,   0.743052464772956,
                   -1.477142755136909,   0.584796229917229,
                   -1.485819261159375,   0.274957181902667;
        for (int i = 0; i < 3; i++) {
            Eigen::Matrix<Scalar, 3, 1> b = basis.row(i); // basis
            Eigen::Matrix<Scalar, 3, 1> cb = elem.L(LocalCoordinates(q[i])); // Caribou basis

            Eigen::Matrix<Scalar, 3, 1> db = deris.row(i); // basis
            Eigen::Matrix<Scalar, 3, 1> cdb = elem.dL(LocalCoordinates(q[i])); // Caribou basis

            EXPECT_MATRIX_NEAR(cb, b, 1e-3);
            EXPECT_MATRIX_NEAR(cdb, db, 1e-3);
            Eigen::Matrix<Scalar, 2, 1> j = jacob.row(i);
            EXPECT_MATRIX_NEAR(elem.jacobian(LocalCoordinates(q[i])), j, 1e-3);
        }
    }
}
