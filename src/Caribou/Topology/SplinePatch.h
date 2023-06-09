#pragma once

#include <Caribou/config.h>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Topology/BaseSplinePatch.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Geometry/NurbsSurf.h>

#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <array>
#include <algorithm>

namespace caribou::topology {

template<typename dtype>
using Vector = Eigen::Matrix<dtype, Eigen::Dynamic, 1>;

template<typename dtype>
using Matrix = Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Double_Vector = Vector<FLOATING_POINT_TYPE>;
using Int_Vector = Vector<int>;
using Double_Matrix = Matrix<FLOATING_POINT_TYPE>;
using Int_Matrix = Matrix<int>;
using USInt_Matrix = Matrix<UNSIGNED_INTEGER_TYPE>;

template<typename Derived>
struct BaseEigenSplineNodesHolder {
    using Scalar = typename Eigen::MatrixBase<Derived>::Scalar ;
    using MatrixType = Derived;

    /*! Copy constructor */
    BaseEigenSplineNodesHolder(const BaseEigenSplineNodesHolder & other)
            : p_nodes(other.p_nodes)
    {}

    /*! Move constructor */
    BaseEigenSplineNodesHolder(BaseEigenSplineNodesHolder && other) noexcept
    : BaseEigenSplineNodesHolder() {
        this->p_nodes.swap(other.p_nodes);
    }

    /*! copy-and-swap assigment (valid for both copy and move assigment) */
    auto operator=(BaseEigenSplineNodesHolder other) noexcept -> BaseEigenSplineNodesHolder & {
        this->p_nodes.swap(other.p_nodes);
        return *this;
    }

    template<typename Index>
    auto node(Index && index) const -> auto {return this->p_nodes.row(index);}

    template<typename Index>
    auto node(Index && index) -> auto {return this->p_nodes.row(index);}

    template<typename Size1>
    auto resize(Size1 && n) -> auto {return this->p_nodes.resize(std::forward<Size1>(n), this->p_nodes.cols());}

    auto size() const -> auto {return this->p_nodes.rows();}

    BaseEigenSplineNodesHolder()
            : p_nodes()
    {}

protected:
    MatrixType p_nodes;
};

template <typename T>
struct EigenSplineNodesHolder {};

template<typename Scalar_t, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct EigenSplineNodesHolder<Eigen::Matrix<Scalar_t, Rows, Cols, Options, MaxRows, MaxCols>>
: public BaseEigenSplineNodesHolder<Eigen::Matrix<Scalar_t, Rows, Cols, Options, MaxRows, MaxCols>>
{
    using Base = BaseEigenSplineNodesHolder<Eigen::Matrix<Scalar_t, Rows, Cols, Options, MaxRows, MaxCols>>;
    using Scalar = typename Base::Scalar;
    using MatrixType = typename Base::MatrixType;

    using Base::Base;
};

template <
    unsigned int WorldDimension, typename NodeIndex = UNSIGNED_INTEGER_TYPE,
    typename NodeContainerType = EigenSplineNodesHolder<Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, WorldDimension, (WorldDimension>1?Eigen::RowMajor:Eigen::ColMajor)>>
>
class SplinePatch : public BaseSplinePatch {
public:

//    using NodeIndex = UNSIGNED_INTEGER_TYPE;

    using Element = geometry::NurbsSurf<WorldDimension>;
    //
    using ElementIndices = Eigen::Matrix<NodeIndex, geometry::traits<Element>::NumberOfNodesAtCompileTime, 1>;
    using ElementsIndices = Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime, Eigen::RowMajor>;
    using ElementsKnotrange = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, geometry::traits<Element>::CanonicalDimension * 2, Eigen::RowMajor>; // 4 -> [u1, v1, u2, v2]

    using DynVector = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1>;

    /**
     * Holder type that contains the nodes of the SplinePatch
     */
    using NodeContainer_t = std::decay_t<NodeContainerType>;
    /**
     * Eigen vector to store the node weights of the SplinePatch
     */
//    using DynVector = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1>;

    static_assert(WorldDimension == 1 or WorldDimension == 2 or WorldDimension == 3, "The world dimension must be 1, 2 or 3.");

    using Self = SplinePatch<WorldDimension, NodeIndex, NodeContainer_t>;
    static constexpr INTEGER_TYPE Dimension = WorldDimension;
    using Real = typename NodeContainer_t::Scalar;
    using WorldCoordinates = Eigen::Matrix<Real, 1, Dimension>;

    /*!
     * Default constructor for an empty mesh.
     */
    SplinePatch() : p_nodes {}, p_weights {}, p_indices{}, p_knotranges{} {}

    // ================== Our specialization starts =========================

    explicit SplinePatch(const NodeContainer_t & positions, const DynVector & weights,
                         const ElementsIndices & indices, const DynVector & knot_1, const DynVector & knot_2, const ElementsKnotrange & knotrange) :
        p_nodes(positions), p_weights(weights), p_indices(indices), p_knot_1(knot_1), p_knot_2(knot_2), p_knotranges(knotrange)
    {
        std::cout << "\n1st constructor \n";
    }

    template <typename Derived>
    explicit SplinePatch(const Eigen::MatrixBase<Derived> & positions, const DynVector & weights,
                         const ElementsIndices & indices, const DynVector & knot_1, const DynVector & knot_2, const ElementsKnotrange & knotrange)
    {
        static_assert(
            Eigen::MatrixBase<Derived>::ColsAtCompileTime == Dimension or Eigen::MatrixBase<Derived>::ColsAtCompileTime == Eigen::Dynamic,
            "The number of columns at compile time should match the Dimension of the mesh, or by dynamic (known at compile time)."
        );

        // The number of columns must equal the World dimension
        caribou_assert(positions.cols() == Dimension);
        // Do the copy
        const auto n = positions.rows();
        p_nodes.resize(n);
        p_weights.resize(n); // Resizing the p_weight vector
        for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
            auto node = this->p_nodes.node(i);
            p_weights(i) = weights(i); // Weight holder
            for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                node[j] = positions(i,j);
            }
        }
        p_indices = indices;
        p_knot_1 = knot_1;
        p_knot_2 = knot_2;
        p_knotranges = knotrange;


    }

    explicit SplinePatch(const Double_Matrix & positions, const Double_Vector & weights, const Matrix<NodeIndex> & indices,
                         const DynVector & knot_1, const DynVector & knot_2, const Double_Matrix & knotrange)
    {
        // The number of columns must equal the World dnimension
        caribou_assert(positions.cols() == Dimension);

        // Do the copy
        const auto n = positions.rows();
        p_nodes.resize(n);
        p_weights.resize(n);
        for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
            auto node = this->p_nodes.node(i);
            p_weights(i) = weights(i); // Weight holder
            for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                node[j] = positions(i,j);
            }
        }
        p_indices = indices;
        p_knot_1 = knot_1;
        p_knot_2 = knot_2;
        p_knotranges = knotrange;
    }

    explicit SplinePatch(const std::vector<WorldCoordinates> & positions, const Double_Vector & weights, const Matrix<NodeIndex> & indices,
                         const DynVector & knot_1, const DynVector & knot_2, const Double_Matrix & knotrange)
    {
        std::cout << "\n5th constructor \n";
        // The number of columns must equal the World dimension
        caribou_assert(positions[0].size() == Dimension);

        // Do the copy
        const int n = positions.size();
        p_nodes.resize(n);
        p_weights.resize(n); // Resizing the p_weight vector
        for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
            auto node = this->p_nodes.node(i);
            p_weights(i) = weights(i); // Weight holder
            for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                node[j] = positions[i][j];
            }
        }
        p_indices = indices;
        p_knot_1 = knot_1;
        p_knot_2 = knot_2;
        p_knotranges = knotrange;
    }

    // ================== Our specialization ends ===========================


    /*! Copy constructor */
    SplinePatch(const SplinePatch & other)
    : p_nodes (other.p_nodes), p_weights (other.p_weights), p_indices(other.p_indices),
      p_knot_1(other.p_knot_1), p_knot_2(other.p_knot_2), p_knotranges(other.p_knotranges){ }

    /*! Move constructor */
    SplinePatch(SplinePatch && other) noexcept
    : SplinePatch() {
        swap(*this, other);
    }

    /*! copy-and-swap assigment (valid for both copy and move assigment) */
    auto operator=(SplinePatch other) noexcept -> SplinePatch &
    {
        swap(*this, other);
        return *this;
    }

    /*!
     * \copydoc caribou::topology::BaseMesh::dimension
     */
    [[nodiscard]]
    inline auto dimension() const -> UNSIGNED_INTEGER_TYPE final {
        return Dimension;
    };

    /*!
     * Get the number of nodes of the mesh.
     */
    [[nodiscard]]
    inline auto number_of_nodes() const -> UNSIGNED_INTEGER_TYPE final {return p_nodes.size();};

    [[nodiscard]]
    inline auto canonical_dimension() const -> UNSIGNED_INTEGER_TYPE final {
        return geometry::traits<Element>::CanonicalDimension;
    }

    [[nodiscard]]
    inline auto number_of_nodes_per_elements() const -> UNSIGNED_INTEGER_TYPE final {
            return geometry::traits<Element>::NumberOfNodesAtCompileTime;
    }

    [[nodiscard]]
    inline auto number_of_elements() const -> UNSIGNED_INTEGER_TYPE final {
        return p_indices.rows();
    }


    // =============================== Domain helpers Start ===========================

    inline auto indices(void) const -> ElementsIndices {
        return p_indices;
    }
    inline auto element_indices(const UNSIGNED_INTEGER_TYPE & index) const {
        caribou_assert(index < number_of_elements() and "Trying to get an element that does not exists.");

        return Eigen::Map<const Eigen::Matrix<NodeIndex, 1, geometry::traits<Element>::NumberOfNodesAtCompileTime>, Eigen::Unaligned, Eigen::Stride<1, Eigen::Dynamic>>(
            p_indices.row(index).data(), {1, p_indices.innerStride()}
        );
    }

    void print_spline_patch(void) const {
//        std::cout << "Control points : \n" << p_nodes  << "\n";
        std::cout << "Weights \n" << p_weights.transpose() << "\n";
        std::cout << "Indices \n" << p_indices << "\n";
        std::cout << "Knot spans \n" << p_knotranges << "\n";
        std::cout << "Knot U  : " << p_knot_1.transpose() << "\n";
        std::cout << "Knot V  : " << p_knot_2.transpose() << "\n";

    }

    inline auto element(const UNSIGNED_INTEGER_TYPE & element_id) const -> Element {
        caribou_assert(element_id < number_of_elements());

        using NodeMatrix = typename geometry::Element<Element>::template Matrix<geometry::traits<Element>::NumberOfNodesAtCompileTime, Dimension>;

        DynVector knot1;
        DynVector knot2;
        DynVector weights(geometry::traits<Element>::NumberOfNodesAtCompileTime);
        DynVector knot_span(geometry::traits<Element>::CanonicalDimension * 2);

        NodeMatrix node_positions;
        if constexpr (geometry::traits<Element>::NumberOfNodesAtCompileTime == caribou::Dynamic) {
            node_positions.resize(number_of_nodes_per_elements(), Dimension);
        }

        const auto node_indices = element_indices(element_id);
        for (std::size_t i = 0; i < static_cast<std::size_t>(node_indices.size()); ++i) {
            auto p1 = node_positions.row(i);
            auto p2 = this->node(node_indices[i]);
            p1[0] = p2[0];
            if constexpr (Dimension > 1) {
                p1[1] = p2[1];
            }
            if constexpr(Dimension > 2) {
                p1[2] = p2[2];
            }
            weights[i] = p_weights[node_indices[i]];
        }
        knot1 = this->p_knot_1;
        knot2 = this->p_knot_2;
        knot_span = this->p_knotranges.row(element_id);
        return Element(node_positions, knot1, knot2, weights, knot_span);
    }

    inline auto element_knotranges(const UNSIGNED_INTEGER_TYPE & index) const {
        caribou_assert(index < number_of_elements() and "Trying to get an element that does not exists.");

        return Eigen::Map<const Eigen::Matrix<FLOATING_POINT_TYPE, 1, geometry::traits<Element>::CanonicalDimension * 2>, Eigen::Unaligned, Eigen::Stride<1, Eigen::Dynamic>>(
            p_knotranges.row(index).data(), {1, p_knotranges.innerStride()}
        );
    }


    // =============================== Domain helpers end ===========================

    /*!
     * \copydoc caribou::topology::BaseMesh::node
     */
    [[nodiscard]]
    inline auto node(const UNSIGNED_INTEGER_TYPE & node_id) const -> Eigen::Vector3d override {
        const auto p = position(node_id);
        Eigen::Vector3d n;
        n[0] = p[0];
        if constexpr (Dimension > 1) {
            n[1] = p[1];
        }
        if constexpr (Dimension > 2) {
            n[2] = p[2];
        }
        return n;
    }

    /*!
     * Get weight correspond to the node_id
     */
    [[nodiscard]]
    inline auto weight(const UNSIGNED_INTEGER_TYPE & node_id) const -> FLOATING_POINT_TYPE {
        caribou_assert(static_cast<Eigen::Index>(node_id) < p_weights.size()
                       && "Trying to access the weight at a node index greater "
                          "than the number of nodes in the mesh.");
        return p_weights(node_id);
    }

    inline  auto all_positions(void) const ->Double_Matrix {
        Double_Matrix positions;
        positions.resize(number_of_nodes(), Dimension);
        for (int i = 0; i < static_cast<int>(number_of_nodes()); ++i) {
            auto node = p_nodes.node(i);
            for (int j = 0; j < static_cast<int>(Dimension); ++j) {
                positions(i,j) = static_cast<FLOATING_POINT_TYPE>(node(j));
            }

        }
        return positions;
    }

    inline  auto positions(void) const ->NodeContainer_t {
        return p_nodes;
    }
    /*!
    * Get the position coordinates of a node from its index.
    */
    [[nodiscard]]
    inline auto position(UNSIGNED_INTEGER_TYPE index) const {
        caribou_assert(static_cast<Eigen::Index>(index) < p_nodes.size()
                       && "Trying to access the position vector at a node index greater "
                          "than the number of nodes in the mesh.");
        return p_nodes.node(index);
    }

    /*!
     * Get an array of positions from the given indices.
     * @tparam IntegerType Integer type of the indices passed
     * @param indices [IN] The indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     *
     * Example:
     * \code{.cpp}
     * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, 2, Dimension, Eigen::RowMajor>;
     * Nodes positions = mesh.positions({55, 62});
     * auto p55 = positions.row(0); // Position vector of the node id 55
     * \endcode
     *
     */
    template <typename IntegerType, std::size_t N>
    inline auto positions(const IntegerType(&indices)[N]) const {
        Eigen::Matrix<Real, N, Dimension, (Dimension>1?Eigen::RowMajor:Eigen::ColMajor)> positions;
        for (std::size_t i = 0; i < N; ++i) {
            caribou_assert(indices[i] < static_cast<IntegerType>(number_of_nodes()));
            positions.row(i) = p_nodes.node(indices[i]);
        }
        return positions;
    }

    /*!
     * Get an array of positions from the given indices.
     * @tparam IntegerType Integer type of the indices passed
     * @param indices [IN] The indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     *
     * Example:
     * \code{.cpp}
     * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, 2, Dimension, Eigen::RowMajor>;
     * std::array<int, 2> indices = {55, 62};
     * Nodes positions = mesh.positions(indices);
     * auto p55 = positions.row(0); // Position vector of the node id 55
     * \endcode
     *
     */
    template <typename IntegerType, std::size_t N>
    inline auto positions(const std::array<IntegerType, N> & indices) const {
        return positions<IntegerType, N>(indices.data());
    }

    /*!
     * Get an array of positions from the given indices.
     * @tparam IntegerType Integer type of the indices passed
     * @param indices [IN] he indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     *
     * Example. :
     * \code{.cpp}
     * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Dimension, Eigen::RowMajor>;
     * std::vector<int> indices = {55, 62};
     * Nodes positions = mesh.positions(indices);
     * auto p55 = positions.row(0); // Position vector of the node id 55
     * \endcode
     *
     */
    template <typename IntegerType>
    inline auto positions(const std::vector<IntegerType> & indices) const {
        const auto number_of_indices = indices.size();
        Eigen::Matrix<Real, Eigen::Dynamic, Dimension, (Dimension == 1) ? Eigen::ColMajor : Eigen::RowMajor> positions;
        positions.resize(number_of_indices, Dimension);
        for (std::size_t i = 0; i < number_of_indices; ++i) {
            positions.row(i) = p_nodes.node(indices[i]);
        }
        return positions;
    }

    /*!
     * Get an array of positions from the given indices.
     * @tparam EigenDerived The type of Eigen matrix/array containing the indices.
     * @param indices [IN] he indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     *
     * Example:
     * \code{.cpp}
     * using Nodes = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Dimension, Eigen::RowMajor>;
     * Eigen::VectorXi indices(2);
     * indices << 55, 62;
     * Nodes positions = mesh.positions(indices);
     * auto p55 = positions.row(0); // Position vector of the node id 55
     * \endcode
     *
     */
    template <typename EigenDerived>
    inline auto positions(const Eigen::EigenBase<EigenDerived> & indices) const {
        static_assert(EigenDerived::RowsAtCompileTime == 1 or EigenDerived::ColsAtCompileTime == 1,
            "Only vector type matrix can be used as the indices container.");
        const auto number_of_indices = indices.size();
        Eigen::Matrix<Real,
                      EigenDerived::RowsAtCompileTime == 1 ? EigenDerived::ColsAtCompileTime : EigenDerived::RowsAtCompileTime,
                      Dimension,
                      (Dimension>1 ? Eigen::RowMajor : Eigen::ColMajor)> positions;
        positions.resize(number_of_indices, Dimension);
        for (Eigen::Index i = 0; i < number_of_indices; ++i) {
            positions.row(i) = p_nodes.node(indices.derived()[i]);
        }
        return positions;
    }


    // ===================== Weights ===============================



    /*!
     * Get an array of weights from the given indices.
     * @tparam IntegerType Integer type of the indices passed
     * @param indices [IN] The indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     */
    template <typename IntegerType, std::size_t N>
    inline auto weights(const IntegerType(&indices)[N]) const {
        Eigen::Matrix<FLOATING_POINT_TYPE, N, 1> weights;
        for (std::size_t i = 0; i < N; ++i) {
            caribou_assert(indices[i] < static_cast<IntegerType>(number_of_nodes()));
            weights(i) = p_weights(indices[i]);
        }
        return weights;
    }

    /*!
     * Get an array of positions from the given indices.
     * @tparam IntegerType Integer type of the indices passed
     * @param indices [IN] The indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     */
    template <typename IntegerType, std::size_t N>
    inline auto weights(const std::array<IntegerType, N> & indices) const {
        return weights<IntegerType, N>(indices.data());
    }

    /*!
     * Get an array of positions from the given indices.
     * @tparam IntegerType Integer type of the indices passed
     * @param indices [IN] he indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     */
    template <typename IntegerType>
    inline auto weights(const std::vector<IntegerType> & indices) const {
        const auto number_of_indices = indices.size();
        Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1> weights;
        weights.resize(number_of_indices);
        for (std::size_t i = 0; i < number_of_indices; ++i) {
            weights(i) = p_weights(indices[i]);
        }
        return weights;
    }

    /*!
     * Get an array of positions from the given indices.
     * @tparam EigenDerived The type of Eigen matrix/array containing the indices.
     * @param indices [IN] he indices from which the positions array will be filled
     * @return The position coordinates of every nodes queried as an Eigen dense matrix.
     */
    template <typename EigenDerived>
    inline auto weights(const Eigen::EigenBase<EigenDerived> & indices) const {
        static_assert(EigenDerived::RowsAtCompileTime == 1 or EigenDerived::ColsAtCompileTime == 1,
            "Only vector type matrix can be used as the indices container.");
        const auto number_of_indices = indices.size();
        Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, 1> weights;
        weights.resize(number_of_indices);
        for (Eigen::Index i = 0; i < number_of_indices; ++i) {
            weights(i) = p_weights(indices.derived()[i]);
        }
        return weights;
    }

    inline auto aa() const {
        return 55;
    }
    inline auto weights() const {
        return p_weights;
    }
    inline auto all_weights() const {
        return p_weights;
    }
    // ===================== Weights ===============================

    // ===================== knots  ===============================
    inline auto size_knot_1() const {
        return p_knot_1.size();
    }

    inline auto size_knot_2() const {
        return p_knot_2.size();
    }

    inline auto knot_1() const {
        return p_knot_1;
    }

    inline auto knot_2() const {
        return p_knot_2;
    }
    // ===================== knots ===============================

    /*! Swap the data of two meshes */
    friend void swap(SplinePatch & first, SplinePatch& second) noexcept
    {
        // enable ADL
        using std::swap;

        swap(first.p_nodes, second.p_nodes);
        swap(first.p_weights, second.p_weights);
        swap(first.p_indices, second.p_indices);
        swap(first.p_knot_1, second.p_knot_1);
        swap(first.p_knot_2, second.p_knot_2);
        swap(first.p_knotranges, second.p_knotranges);

    }

private:

    /**
     * Container of the mesh node positions
     */
    NodeContainer_t p_nodes;

    /**
     * Container of the mesh node weights
     */
    DynVector p_weights;

    ElementsIndices p_indices;

    // Store U direction knot vector
    DynVector p_knot_1;

    // Store U direction knot vector
    DynVector p_knot_2;

    /// knot_ranges containing the each elements knot range in both parametric directions
    ElementsKnotrange p_knotranges;


};
} // namespace caribou::topology


namespace caribou::topology {
    // -------------------------------------------------------------------------------------
    //                                    IMPLEMENTATION
    // -------------------------------------------------------------------------------------
    // Implementation of methods that can be specialized later for an explicit container type
    // -------------------------------------------------------------------------------------


    // -------------------------------------------------------------------------------------
    //                               TEMPLATE DEDUCTION GUIDES
    // -------------------------------------------------------------------------------------
    // Template deduction guides that can help the compiler to automatically determine the
    // template parameters of the Mesh from the constructor used.
    // -------------------------------------------------------------------------------------
    template <int WorldDimension, typename NodeIndex = UNSIGNED_INTEGER_TYPE, typename Real>
    SplinePatch(const std::vector<Eigen::Matrix<Real, 1, WorldDimension>> &) ->
    SplinePatch <
        WorldDimension, NodeIndex,
        EigenSplineNodesHolder<Eigen::Matrix<Real, Eigen::Dynamic, WorldDimension, (WorldDimension>1?Eigen::RowMajor:Eigen::ColMajor)>>
    >;

    template <typename Derived, typename NodeIndex = UNSIGNED_INTEGER_TYPE>
    SplinePatch(const Eigen::MatrixBase<Derived> &) ->
    SplinePatch <
        Eigen::MatrixBase<Derived>::ColsAtCompileTime, NodeIndex,
        EigenSplineNodesHolder<
            Eigen::Matrix<
                typename Eigen::MatrixBase<Derived>::Scalar,
                Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                (Eigen::MatrixBase<Derived>::ColsAtCompileTime>1?Eigen::RowMajor:Eigen::ColMajor)
            >
        >
    >;
} // namespace caribou::topology
