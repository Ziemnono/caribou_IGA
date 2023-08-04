#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseNurbsSurf.h>
#include <Eigen/Core>
#include <Caribou/Geometry/Quadrature.h>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct NurbsSurf;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<NurbsSurf <_Dimension>> {
     // static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; // NurbsSurf degree
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2; // Two parametric direction.
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = caribou::Dynamic; // NurbsSurf control points (p+1)*(q+1)
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = caribou::Dynamic; // Gauss points per element
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
struct NurbsSurf: public BaseNurbsSurf<NurbsSurf <_Dimension>> {
    // Types

    using Scalar = FLOATING_POINT_TYPE;

    using Base = BaseNurbsSurf<NurbsSurf<_Dimension>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <INTEGER_TYPE Rows, INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    using Dyn_Matrix = typename Base::Dyn_Matrix;
    using Dyn_Vector = typename Base::Dyn_Vector;

//    static constexpr UNSIGNED_INTEGER_TYPE degree_u = 2; //NurbsSurf
//    static constexpr UNSIGNED_INTEGER_TYPE degree_v = 2; //NurbsSurf
    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime; //$$$NurbsSurf
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime; //$$$NurbsSurf


    // Constructors
    using Base::Base;
    NurbsSurf() : Base() {}

    void print(void){
        std::cout << "Control Points : \n" << this->p_nodes.transpose() << "\n";
        std::cout << "Weight Points  : \n" << this->p_weights.transpose() << "\n";
        std::cout << "Knot Vector 1  : \n" << this->p_knot_1.transpose() << "\n";
        std::cout << "Knot Vector 2  : \n" << this->p_knot_2.transpose() << "\n";
        std::cout << "Knot span      : \n" << this->p_knot_span.transpose() << "\n";

    }

    auto aa(void) const -> int{
        return 555;
    }

private:
    // Implementations
    friend struct Element<NurbsSurf <_Dimension>>;
    friend struct BaseNurbsSurf<NurbsSurf <_Dimension>>;


    auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {

        UNSIGNED_INTEGER_TYPE degree_u = this->p_degrees[0];
        UNSIGNED_INTEGER_TYPE degree_v = this->p_degrees[1];

        // U direction basis
        Scalar u = 0.5 * (xi[0]*(this->p_knot_span[2] - this->p_knot_span[0]) + (this->p_knot_span[2] + this->p_knot_span[0]));
        Scalar v = 0.5 * (xi[1]*(this->p_knot_span[3] - this->p_knot_span[1]) + (this->p_knot_span[3] + this->p_knot_span[1]));

        NURBS_utils Nutils;

        auto knot1 = this->p_knot_1;
        auto basis_u = Nutils.bspbasisfun(u, degree_u, knot1);

        auto knot2 = this->p_knot_2;
        auto basis_v = Nutils.bspbasisfun(v, degree_v, knot2);

        Dyn_Matrix m = basis_v * basis_u.transpose();

        Eigen::Map<Vector<NumberOfNodesAtCompileTime>> basis(m.data(), m.size());
        basis = basis.array() * this->p_weights.array();
        basis = basis * (1/basis.sum());
        return basis;

    };

    inline auto get_dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {
        Scalar u = 0.5 * (xi[0]*(this->p_knot_span[2] - this->p_knot_span[0]) + (this->p_knot_span[2] + this->p_knot_span[0]));
        Scalar v = 0.5 * (xi[1]*(this->p_knot_span[3] - this->p_knot_span[1]) + (this->p_knot_span[3] + this->p_knot_span[1]));

        UNSIGNED_INTEGER_TYPE degree_u = this->p_degrees[0];
        UNSIGNED_INTEGER_TYPE degree_v = this->p_degrees[1];

        UNSIGNED_INTEGER_TYPE num_nodes_u = degree_u+1;
        UNSIGNED_INTEGER_TYPE num_nodes_v = degree_v+1;
        UNSIGNED_INTEGER_TYPE num_nodes = num_nodes_u * num_nodes_v;

//        std::cout << "\nu - " <<  u << "\n";
//        std::cout << "v - " <<  v << "\n";
        auto knot1 = this->p_knot_1;
        auto knot2 = this->p_knot_2;

        Dyn_Vector Nbasis_u(num_nodes_u);
        Dyn_Vector Nbasis_v(num_nodes_v);

        Dyn_Vector Nderis_u(num_nodes_u);
        Dyn_Vector Nderis_v(num_nodes_v);

        NURBS_utils Nutils;

        // U direction basis & derivative
        Nbasis_u = Nutils.bspbasisfun(u, degree_u, knot1);
        Nderis_u = Nutils.basis_deri(u, degree_u, knot1);

        // V direction basis & derivative
        Nbasis_v = Nutils.bspbasisfun(v, degree_v, knot2);
        Nderis_v = Nutils.basis_deri(v, degree_v, knot2);

        Dyn_Matrix global_storage(num_nodes_v,num_nodes_u);
        global_storage =  Nbasis_v * Nbasis_u.transpose();
//        std::cout << "Golbal basis Matrix \n" <<  global_storage << "\n";

        Eigen::Map<Dyn_Vector> basis(global_storage.data(), global_storage.size());
//        std::cout << "Golbal basis \n" <<  basis.transpose() << "\n";

        Dyn_Matrix deris(num_nodes,2);
        // Bv_dBu
        Dyn_Matrix deri_m1(num_nodes_v,num_nodes_u);
        deri_m1 = Nbasis_v * Nderis_u.transpose();
        deris.col(0) = Eigen::Map<Dyn_Vector>(deri_m1.data(), deri_m1.size());

        // dBv_Bu
        Dyn_Matrix deri_m2(num_nodes_v,num_nodes_u);
        deri_m2 = Nderis_v * Nbasis_u.transpose();
        deris.col(1) = Eigen::Map<Dyn_Vector>(deri_m2.data(), deri_m2.size());

//        std::cout << "First deris \n" <<  deris << "\n";

        // Weights * Basis functions component-wise multiplication.
        Dyn_Vector B1 = this->p_weights.array() * basis.array();
        double w1 = B1.sum();

        // Calculating dB_du and dB_dv
        Dyn_Vector dB_du_B2 = this->p_weights.array() * deris.col(0).array();
        double dB_du_w2 = dB_du_B2.sum();
        deris.col(0) = dB_du_B2/w1 - B1 * dB_du_w2/(w1*w1);

        // Calculating  dB_dv
        Dyn_Vector dB_dv_B2 = this->p_weights.array() * deris.col(1).array();
        double dB_dv_w2 = dB_dv_B2.sum();
        deris.col(1) = dB_dv_B2/w1 - B1 * dB_dv_w2/(w1*w1);

        return deris;
    };



    inline auto get_gauss_nodes() const -> const auto & {

        const int gauss_u = this->p_degrees[0]+1; // Gauss points in U direction
        const int gauss_v = this->p_degrees[1]+1; // Gauss points in V direction
        const int no_gauss_nodes = (gauss_u) * (gauss_v);

        static std::vector<GaussNode> gauss_nodes;
        gauss_nodes.resize(no_gauss_nodes);
        std::cout << "Gauss nodes is called \n";
        Quadrature quad_u(gauss_u);
        Quadrature quad_v(gauss_v);
        int count = 0;
        for (int j = 0; j < gauss_u; ++j) {
            for (int i = 0; i < gauss_v; ++i) {
                gauss_nodes[count] = GaussNode{LocalCoordinates(quad_u.get_point(i), quad_v.get_point(j)), quad_u.get_weight(i) * quad_v.get_weight(j)};
                count++;
            }
        }

        return gauss_nodes;
    }
};

}
