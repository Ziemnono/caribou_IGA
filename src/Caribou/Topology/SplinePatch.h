#pragma once

#include <Caribou/config.h>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Topology/BaseSplinePatch.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Geometry/BezierSurf.h>

#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <array>
#include <algorithm>

namespace caribou::topology {

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
        unsigned int WorldDimension,
        typename NodeContainerType = EigenSplineNodesHolder<Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, WorldDimension, (WorldDimension>1?Eigen::RowMajor:Eigen::ColMajor)>>
    >
    class SplinePatch : public BaseSplinePatch {
    public:

        using NodeIndex = UNSIGNED_INTEGER_TYPE;

        using Element = geometry::BezierSurf<WorldDimension>;
        //
        using ElementIndices = Eigen::Matrix<NodeIndex, geometry::traits<Element>::NumberOfNodesAtCompileTime, 1>;
        using ElementsIndices = Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::CanonicalDimension, Eigen::RowMajor>; // 4 -> [u1, v1, u2, v2]
        using ElementsKnotrange = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using ElementsExtraction = std::vector<ElementsKnotrange>;


        /**
         * Holder type that contains the nodes of the SplinePatch
         */
        using NodeContainer_t = std::decay_t<NodeContainerType>;
        /**
         * Eigen vector to store the node weights of the SplinePatch
         */
        using WeightsContainer = Eigen::Matrix<FLOATING_POINT_TYPE, 1, Eigen::Dynamic, Eigen::RowMajor>;

        static_assert(WorldDimension == 1 or WorldDimension == 2 or WorldDimension == 3, "The world dimension must be 1, 2 or 3.");

        using Self = SplinePatch<WorldDimension, NodeContainer_t>;
        static constexpr INTEGER_TYPE Dimension = WorldDimension;
        using Real = typename NodeContainer_t::Scalar;
        using WorldCoordinates = Eigen::Matrix<Real, 1, Dimension>;

        /*!
         * Default constructor for an empty mesh.
         */
        SplinePatch() : p_nodes {}, p_weights {}, p_buffer{}, p_knotranges{}, p_extraction{} {}

        /*!
         * Virtual default destructor
         */
        virtual ~SplinePatch() = default;

        /**
         * Construct the unstructured mesh with a set of point positions from a position vector.
         * @param positions A reference to a std::vector containing the position of the nodes
         *
         * @note A copy is made of all nodes position vectors is made.
         */
        explicit SplinePatch(const std::vector<WorldCoordinates> & positions, const std::vector<FLOATING_POINT_TYPE> & weights)
        {
            const auto n = positions.size();
            const auto n_weights = weights.size();
            if (n != n_weights){
                throw("Number of nodes should match with number of weights");
            }
            p_nodes.resize(n);
            p_weights.resize(n);
            for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
                auto node = this->p_nodes.node(i);
                p_weights(i) = weights[i];
                for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                    node[j] = positions[i][j];
                }
            }
        }

        explicit SplinePatch(const std::vector<WorldCoordinates> & positions, const WeightsContainer & weights)
        {
            const auto n = positions.size();
            const auto n_weights = weights.size();
            if (n != n_weights){
                throw("Number of nodes should match with number of weights");
            }
            p_nodes.resize(n);
            p_weights.resize(n);
            for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
                auto node = this->p_nodes.node(i);
                p_weights(i) = weights(i);
                for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                    node[j] = positions[i][j];
                }
            }
        }

        /**
         * Construct the unstructured mesh with an Eigen matrix containing the position vector
         * (NxD with N nodes of D world dimension).
         *
         * @param positions A reference to a NxD matrix containing the position vector
         */
        template <typename Derived>
        explicit SplinePatch(const Eigen::MatrixBase<Derived> & positions, const WeightsContainer & weights)
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
        }

        // ================== Our specialization starts =========================

        explicit SplinePatch(const NodeContainer_t & positions, const WeightsContainer & weights,
                             const ElementsIndices & indices, const ElementsKnotrange & knotrange, const ElementsExtraction & extraction) :
            p_nodes (positions), p_weights(weights), p_buffer(indices), p_knotranges(knotrange), p_extraction(extraction)
        { }


        explicit SplinePatch(const std::vector<WorldCoordinates> & positions, const WeightsContainer & weights,
                             const ElementsIndices & indices, const ElementsKnotrange & knotrange, const ElementsExtraction & extraction)
        {
            const auto n = positions.size();
            const auto n_weights = weights.size();
            if (n != n_weights){
                throw("Number of nodes should match with number of weights");
            }
            p_nodes.resize(n);
            p_weights.resize(n);
            for (std::size_t i = 0; i < static_cast<std::size_t> (n); ++i) {
                auto node = this->p_nodes.node(i);
                p_weights(i) = weights(i);
                for (std::size_t j = 0; j < static_cast<std::size_t>(Dimension); ++j) {
                    node[j] = positions[i][j];
                }
            }
            p_buffer = indices;
            p_knotranges = knotrange;
            p_extraction = extraction;
        }

        template <typename Derived>
        explicit SplinePatch(const Eigen::MatrixBase<Derived> & positions, const WeightsContainer & weights,
                             const ElementsIndices & indices, const ElementsKnotrange & knotrange, const ElementsExtraction & extraction)
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
            p_buffer = indices;
            p_knotranges = knotrange;
            p_extraction = extraction;
        }

        // ================== Our specialization ends ===========================


        /*! Copy constructor */
        SplinePatch(const SplinePatch & other)
        : p_nodes (other.p_nodes), p_weights (other.p_weights), p_buffer(other.p_buffer), p_knotranges(other.p_knotranges), p_extraction(other.p_extraction){ }

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
        inline auto number_of_nodes_per_elements() const -> UNSIGNED_INTEGER_TYPE {
                return p_buffer.cols();
        }

        [[nodiscard]]
        inline auto number_of_elements() const -> UNSIGNED_INTEGER_TYPE {
            return p_buffer.rows();
        }


        // =============================== Domain helpers Start ===========================

        inline auto element_indices(const UNSIGNED_INTEGER_TYPE & index) const {
            caribou_assert(index < number_of_elements() and "Trying to get an element that does not exists.");

            return Eigen::Map<const Eigen::Matrix<NodeIndex, geometry::traits<Element>::NumberOfNodesAtCompileTime, 1>, Eigen::Unaligned, Eigen::Stride<1, Eigen::Dynamic>>(
                p_buffer.row(index).data(), {1, p_buffer.innerStride()}
            );
        }


        inline auto element(const UNSIGNED_INTEGER_TYPE & element_id) const -> Element {
            caribou_assert(element_id < number_of_elements());

            using NodeMatrix = typename geometry::Element<Element>::template Matrix<geometry::traits<Element>::NumberOfNodesAtCompileTime, Dimension>;
            NodeMatrix node_positions;
            if constexpr (geometry::traits<Element>::NumberOfNodesAtCompileTime == caribou::Dynamic) {
                node_positions.resize(number_of_nodes_per_elements(), Dimension);
            }

            const auto node_indices = element_indices(element_id);
            for (std::size_t i = 0; i < static_cast<std::size_t>(node_indices.size()); ++i) {
                auto p1 = node_positions.row(i);
                auto p2 = this->mesh()->node(node_indices[i]);
                p1[0] = p2[0];
                if constexpr (Dimension > 1) {
                    p1[1] = p2[1];
                }
                if constexpr(Dimension > 2) {
                    p1[2] = p2[2];
                }
            }

            return Element(node_positions);
        }

        inline auto element_knotranges(const UNSIGNED_INTEGER_TYPE & index) const {
            caribou_assert(index < number_of_elements() and "Trying to get an element that does not exists.");

            return Eigen::Map<const Eigen::Matrix<NodeIndex, 2*geometry::traits<Element>::CanonicalDimension, 1>, Eigen::Unaligned, Eigen::Stride<1, Eigen::Dynamic>>(
                p_knotranges.row(index).data(), {1, p_buffer.innerStride()}
            );
        }


        inline auto element_extraction(const UNSIGNED_INTEGER_TYPE & index) const {
            caribou_assert(index < number_of_elements() and "Trying to get an element that does not exists.");
            return p_extraction[index];
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
         * Example:
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
            Eigen::Matrix<FLOATING_POINT_TYPE, 1, N, Eigen::RowMajor> weights;
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
            Eigen::Matrix<FLOATING_POINT_TYPE, 1, Eigen::Dynamic, Eigen::RowMajor> weights;
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
            Eigen::Matrix<FLOATING_POINT_TYPE, 1, Eigen::Dynamic, Eigen::RowMajor> weights;
            weights.resize(number_of_indices);
            for (Eigen::Index i = 0; i < number_of_indices; ++i) {
                weights(i) = p_weights(indices.derived()[i]);
            }
            return weights;
        }

        // ===================== Weights ===============================

        /*! Swap the data of two meshes */
        friend void swap(SplinePatch & first, SplinePatch& second) noexcept
        {
            // enable ADL
            using std::swap;

            swap(first.p_nodes, second.p_nodes);
            swap(first.p_weights, second.p_weights);
            swap(first.p_buffer, second.p_buffer);
            swap(first.p_knotranges, second.p_knotranges);
            swap(first.p_extraction, second.p_extraction);

        }

    private:

        /**
         * Container of the mesh node positions
         */
        NodeContainer_t p_nodes;

        /**
         * Container of the mesh node weights
         */
        WeightsContainer p_weights;

        /// stored outside of the domain instance), then this buffer is empty.
        ElementsIndices p_buffer;

        /// knot_ranges containing the each elements knot range in both parametric directions
        ElementsKnotrange p_knotranges;

        /// extraction containing the each elements extraction matrix
        ElementsExtraction p_extraction;

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
    template <int WorldDimension, typename Real>
    SplinePatch(const std::vector<Eigen::Matrix<Real, 1, WorldDimension>> &) ->
    SplinePatch <
        WorldDimension,
        EigenSplineNodesHolder<Eigen::Matrix<Real, Eigen::Dynamic, WorldDimension, (WorldDimension>1?Eigen::RowMajor:Eigen::ColMajor)>>
    >;

    template <typename Derived>
    SplinePatch(const Eigen::MatrixBase<Derived> &) ->
    SplinePatch <
        Eigen::MatrixBase<Derived>::ColsAtCompileTime,
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
