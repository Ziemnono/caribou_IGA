#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseBezierSurf.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct BezierSurf;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<BezierSurf <_Dimension>> {
     // static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; // BezierSurf degree
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2; // Two parametric direction.
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 3*3; // BezierSurf control points (p+1)*(q+1)
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 3*3; // Gauss points per element
 };

 /**
  * Bezier Surface two canonicalDimensions
  *
  * \verbatim
  *        v
  *        ^
  *        |
  *  6-----7-----8
  *  |     |     |
  *  |     |     |
  *  3     4---- 5 --> u
  *  |           |
  *  |           |
  *  0-----1-----2
  *  \endverbatim
  *
  * @tparam _Dimension The world coordinates dimension
  */

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct BezierSurf: public BaseBezierSurf<BezierSurf <_Dimension>> {
    // Types

    using Scalar = FLOATING_POINT_TYPE;

    using Base = BaseBezierSurf<BezierSurf <_Dimension>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    static constexpr UNSIGNED_INTEGER_TYPE Degree_1 = 2; //BezierSurf
    static constexpr UNSIGNED_INTEGER_TYPE Degree_2 = 2; //BezierSurf
    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime; //$$$BezierSurf
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime; //$$$BezierSurf

    static constexpr FLOATING_POINT_TYPE canonical_nodes [NumberOfNodesAtCompileTime][CanonicalDimension] = {
    //    u,  v
        {-1, -1}, // Node 0
        {+0, -1}, // Node 1
        {+1, -1}, // Node 2
        {-1, +0}, // Node 3
        {+0, +0}, // Node 4
        {+1, +0}, // Node 5
        {-1, +1}, // Node 6
        {+0, +1}, // Node 7
        {+1, +1}, // Node 8

    };

    // Constructors
    using Base::Base;
    BezierSurf() : Base() {
        if constexpr (Dimension == 2) {
            // WorldCoordinates(x, y, w)
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1]); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1]); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1]); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1]); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(canonical_nodes[4][0], canonical_nodes[4][1]); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(canonical_nodes[5][0], canonical_nodes[5][1]); // Node 50.5 *
            this->p_nodes.row(6) = WorldCoordinates(canonical_nodes[6][0], canonical_nodes[6][1]); // Node 6
            this->p_nodes.row(7) = WorldCoordinates(canonical_nodes[7][0], canonical_nodes[7][1]); // Node 7
            this->p_nodes.row(8) = WorldCoordinates(canonical_nodes[8][0], canonical_nodes[8][1]); // Node 8
        }
        else {
            // WorldCoordinates(x, y, z, w)
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1], 0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1], 0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1], 0); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1], 0); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(canonical_nodes[4][0], canonical_nodes[4][1], 0); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(canonical_nodes[5][0], canonical_nodes[5][1], 0); // Node 5
            this->p_nodes.row(6) = WorldCoordinates(canonical_nodes[6][0], canonical_nodes[6][1], 0); // Node 6
            this->p_nodes.row(7) = WorldCoordinates(canonical_nodes[7][0], canonical_nodes[7][1], 0); // Node 7
            this->p_nodes.row(8) = WorldCoordinates(canonical_nodes[8][0], canonical_nodes[8][1], 0); // Node 8
        };
    }


private:
    // Implementations
    friend struct Element<BezierSurf <_Dimension>>;
    friend struct BaseBezierSurf<BezierSurf <_Dimension>>;

    inline auto bezier_basis( const int & p, const Scalar & xi, const int & ni) const -> Scalar {
        if (p == 0 && ni == 0) {
            return 1.;
        }
        else if (p == 0 && ni != 0) {
            return 0.;
        }
        else{
            if (ni < 0 || ni > p+1) {
                return 0.;
            }
            else{
                return 0.5*(1-xi)*bezier_basis(p-1, xi, ni) + 0.5 * (1+xi) * bezier_basis(p-1, xi, ni-1);
            }
        }
    }

    inline auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        Vector<NumberOfNodesAtCompileTime> L;
        // U direction basis
        for (size_t j=0; j <= Degree_2; ++j){
            for (size_t i = 0; i <= Degree_1; ++i) {
                L[(Degree_1+1)*j+i] = bezier_basis(Degree_1, xi[0], i) * bezier_basis(Degree_2, xi[1], j);
            }
        }
        return L;
    };

    inline auto get_dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {
        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> dL;
        for (size_t j=0; j <= Degree_2; ++j){
            for (size_t i = 0; i <= Degree_1; ++i) {
                // dL_du = B_v * dB_u
                dL((Degree_1+1)*j+i, 0) = 0.5 * Degree_1 * (bezier_basis(Degree_1-1, xi[0], i-1) - bezier_basis(Degree_1-1, xi[0], i))
                                              * bezier_basis(Degree_2, xi[1], j);
                // dL_dv = B_u * dB_v
                dL((Degree_1+1)*j+i, 1) = bezier_basis(Degree_1, xi[0], i) *
                                          0.5 * Degree_2 * (bezier_basis(Degree_2-1, xi[1], j-1) - bezier_basis(Degree_2-1, xi[1], j));
            }
        }
        return dL;
    };



    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0), -1*sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 0
            GaussNode {LocalCoordinates(0,                -1*sqrt(3.0/5.0)), (8.0/9.0) * (5.0/9.0)}, // Node 1
            GaussNode {LocalCoordinates(sqrt(3.0/5.0),    -1*sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 2

            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0),                0), (5.0/9.0) * (8.0/9.0)}, // Node 3
            GaussNode {LocalCoordinates(0,                               0), (8.0/9.0) * (8.0/9.0)}, // Node 4
            GaussNode {LocalCoordinates(sqrt(3.0/5.0),                   0), (5.0/9.0) * (8.0/9.0)}, // Node 5

            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0),    sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 6
            GaussNode {LocalCoordinates(0,                   sqrt(3.0/5.0)), (8.0/9.0) * (5.0/9.0)}, // Node 7
            GaussNode {LocalCoordinates(sqrt(3.0/5.0),       sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 8

        };
        return gauss_nodes;
    }
};

}
