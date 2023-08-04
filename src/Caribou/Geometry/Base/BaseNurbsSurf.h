#pragma once

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<typename Derived>
struct BaseNurbsSurf : public Element<Derived> {
    // Types
    using Base = Element<Derived>;

    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <INTEGER_TYPE Rows, INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    using Scalar = typename Base::Scalar;

    using Dyn_Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Dyn_Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime;

    /** Default empty constructor */
    BaseNurbsSurf() = default;

    /** Constructor from an Eigen matrix containing the positions of the segment's nodes */
    template<typename EigenType>
    explicit BaseNurbsSurf(const Eigen::Matrix<UNSIGNED_INTEGER_TYPE, 2, 1> & degrees, const Eigen::EigenBase<EigenType> & nodes,
                           Dyn_Vector knot_1, Dyn_Vector knot_2, Dyn_Vector weights, Dyn_Vector knot_span) :
        p_degrees(degrees), p_nodes(nodes.derived().template cast<typename Base::Scalar>()), p_knot_1(knot_1), p_knot_2(knot_2),
        p_weights(weights), p_knot_span(knot_span) {}

    inline auto jacobian_papa() const -> Scalar {
        auto Jxi = 0.5 * (p_knot_span[2] - p_knot_span[0]);
        auto Jeta = 0.5 * (p_knot_span[3] - p_knot_span[1]);
        return Jxi * Jeta;
    }

private:
    // Implementations
    friend struct Element<Derived>;
    [[nodiscard]]
    inline auto get_number_of_nodes() const {return p_nodes.rows();}
    inline auto get_number_of_gauss_nodes() const {return (p_degrees[0]+1)*(p_degrees[1]+1);}
    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const {return WorldCoordinates(p_nodes.row(index));};
    inline auto get_nodes() const -> const auto & {return p_nodes;};
    inline auto get_center() const {return Base::world_coordinates(LocalCoordinates(0,0));};
    inline auto get_contains_local(const LocalCoordinates & xi, const FLOATING_POINT_TYPE & eps) const -> bool {
        const auto & u = xi[0];
        const auto & v = xi[1];
        return IN_CLOSED_INTERVAL(-1-eps, u, 1+eps) and IN_CLOSED_INTERVAL(-1-eps, v, 1+eps);
    }

protected:
    Eigen::Matrix<UNSIGNED_INTEGER_TYPE, 2, 1> p_degrees;
    Matrix<NumberOfNodesAtCompileTime, Dimension> p_nodes;
    Dyn_Vector p_knot_1;
    Dyn_Vector p_knot_2;
    Dyn_Vector p_weights;
    Dyn_Vector p_knot_span;
};

}
