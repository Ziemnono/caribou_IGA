#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseBezierCrv.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct BezierCrv;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<BezierCrv <_Dimension>> {
     // static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; // BezierCrv degree
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 1; // Single parametric direction.
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 3; // BezierCrv control points
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 3; // Gauss points per element
 };

/**
 * Single element BezierCrv
 *
 * \verbatim
 * P1 : 0-----+-----1 --> u
 * \endverbatim
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct BezierCrv: public BaseBezierCrv<BezierCrv <_Dimension>> {
    // Types

    using Scalar = FLOATING_POINT_TYPE;

    using Base = BaseBezierCrv<BezierCrv <_Dimension>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; //BezierCrv
    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Degree + 1; //$$$BezierCrv
    static constexpr auto NumberOfGaussNodesAtCompileTime = Degree + 1; //$$$BezierCrv

    // Constructors
    using Base::Base;
    BezierCrv() : Base() {
        if constexpr (Dimension == 1) {
            this->p_nodes[0] = -1; //
            this->p_nodes[1] = 0;  //
            this->p_nodes[2] = +1; //
        } else if constexpr (Dimension == 2) {
            this->p_nodes.row(0) = WorldCoordinates(-1, 0);
            this->p_nodes.row(1) = WorldCoordinates(0, 0); //
            this->p_nodes.row(2) = WorldCoordinates(+1, 0); //
        } else {
            this->p_nodes.row(0) = WorldCoordinates(-1, 0, 0);
            this->p_nodes.row(1) = WorldCoordinates(0, 0, 0); //
            this->p_nodes.row(2) = WorldCoordinates(+1, 0, 0); //
        }
    };


private:
    // Implementations
    friend struct Element<BezierCrv <_Dimension>>;
    friend struct BaseBezierCrv<BezierCrv <_Dimension>>;

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
        for (int i=0; i <= 2; ++i){
            L[i] = bezier_basis(Degree, xi[0], i);
        }
        return L;
    };

    inline auto get_dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {

        // Copying the knotvector to U for easy usage.
        Vector<NumberOfNodesAtCompileTime> dL;
        for (int i=0; i <= 2; ++i){
            dL[i] = 0.5 * Degree * (bezier_basis(Degree-1, xi[0], i-1) - bezier_basis(Degree-1, xi[0], i));
        }
        return dL;
    };



    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0)), 5.0/9.0}, // Node 0
            GaussNode {LocalCoordinates(0)             , 8.0/9.0}, // Node 1
            GaussNode {LocalCoordinates(sqrt(3.0/5.0)) , 5.0/9.0}, // Node 2
        };
        return gauss_nodes;
    }
};

}
