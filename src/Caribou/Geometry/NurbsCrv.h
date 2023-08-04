#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseNurbsCrv.h>
#include <Eigen/Core>
#include <Caribou/Geometry/Quadrature.h>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct NurbsCrv;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<NurbsCrv <_Dimension>> {
     // static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; // NurbsCrv degree
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 1; // 1 parametric direction.
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = caribou::Dynamic; // NurbsCrv control points (p+1)*(q+1)
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = caribou::Dynamic; // Gauss points per element
 };

 /**
  * Bezier Surface two canonicalDimensions
  *
  * \verbatim
  *       --> u
  *  0-----1-----2
  *  \endverbatim
  *
  * @tparam _Dimension The world coordinates dimension
  */

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct NurbsCrv: public BaseNurbsCrv<NurbsCrv <_Dimension>>{
    // Types

    using Scalar = FLOATING_POINT_TYPE;

    using Base = BaseNurbsCrv<NurbsCrv<_Dimension>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <INTEGER_TYPE Rows, INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    using Dyn_Matrix = typename Base::Dyn_Matrix;
    using Dyn_Vector = typename Base::Dyn_Vector;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime; //$$$NurbsCrv
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime; //$$$NurbsCrv

    // Constructors
    using Base::Base;
    NurbsCrv() : Base() {
    }

    void print(void){
        std::cout << "Degree : " << this->p_degree_u;
        std::cout << "Control Points : \n" << this->p_nodes.transpose() << "\n";
        std::cout << "Weight Points  : \n" << this->p_weights.transpose() << "\n";
        std::cout << "Knot Vector  : \n" << this->p_knots.transpose() << "\n";
        std::cout << "Knot span      : \n" << this->p_knot_span.transpose() << "\n";

    }

private:
    // Implementations
    friend struct Element<NurbsCrv <_Dimension>>;
    friend struct BaseNurbsCrv<NurbsCrv <_Dimension>>;


    auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        // U direction basis
        Scalar u = 0.5 * (xi[0]*(this->p_knot_span[1] - this->p_knot_span[0]) + (this->p_knot_span[1] + this->p_knot_span[0]));

        NURBS_utils Nutils;

        auto knots = this->p_knots;
        UNSIGNED_INTEGER_TYPE degree = this->p_degree_u;
        UNSIGNED_INTEGER_TYPE num_nodes = degree + 1;
        Vector<NumberOfNodesAtCompileTime> basis;
        basis.resize(num_nodes,1);
        basis = Nutils.bspbasisfun(u, degree, knots);

        basis = basis.array() * this->p_weights.array();
        basis = basis * (1/basis.sum());
        return basis;

    };

    inline auto get_dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {
//        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> dL;
        Scalar u = 0.5 * (xi[0] * (this->p_knot_span[1] - this->p_knot_span[0]) + (this->p_knot_span[1] + this->p_knot_span[0]));

//        std::cout << "\nu - " <<  u << "\n";
//        std::cout << "v - " <<  v << "\n";
        auto knots = this->p_knots;
        UNSIGNED_INTEGER_TYPE degree = this->p_degree_u;
        UNSIGNED_INTEGER_TYPE num_nodes = degree + 1;

        Dyn_Vector N;
        N.resize(num_nodes,1);
        Dyn_Vector Nderis;
        Nderis.resize(num_nodes, 1);

        NURBS_utils Nutils;
        N = Nutils.bspbasisfun(u, degree, knots);
        Nderis = Nutils.basis_deri(u, degree, knots);

        Scalar w = (N.array() * this->p_weights.array()).sum();
        Scalar dwdxi = (Nderis.array() * this->p_weights.array()).sum();

        Dyn_Vector fac = this->p_weights.array()/(w*w);
        Nderis = (Nderis*w - N*dwdxi).array() * fac.array();
        return Nderis;
    };



    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes;
        auto no_gauss_nodes = this->p_degree_u;
        gauss_nodes.resize(no_gauss_nodes);

        Quadrature q(no_gauss_nodes);

        for (int i = 0; i < no_gauss_nodes; ++i) {
            gauss_nodes(i) = GaussNode{LocalCoordinates(q.get_point(i)), q.get_weight(i)};
        }

        return gauss_nodes;
    }
};

}
