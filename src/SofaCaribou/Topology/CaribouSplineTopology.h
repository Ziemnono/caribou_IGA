#pragma once

#include <SofaCaribou/config.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Topology/SplinePatch.h>
#include <unsupported/Eigen/CXX11/Tensor>
DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/State.h>
#include <sofa/core/topology/Topology.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/BaseMatrix.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
#include <sofa/helper/fixed_array.h>
#include <sofa/helper/vector.h>
namespace sofa::type {
template <typename T> using vector = ::sofa::helper::vector<T>;
template <typename T, std::size_t N> using fixed_array = ::sofa::helper::fixed_array<T, N>;
}
#endif

namespace SofaCaribou::topology {
// Traits to get the Sofa vector type from the dimension
template <std::size_t Dim> struct SofaVecType {};
template <> struct SofaVecType<1> { using Type = sofa::defaulttype::Vec1Types; };
template <> struct SofaVecType<2> { using Type = sofa::defaulttype::Vec2Types; };
template <> struct SofaVecType<3> { using Type = sofa::defaulttype::Vec3Types; };

/**
 * The CaribouSplineTopology is a wrapper object over an instance of caribou::topology::Domain.
 *
 * It is templated over a caribou::geometry::Element parameter, which means that the topology
 * covered by this component is independent of the element type. There are two ways to construct
 * its topology:
 *
 * 1. By giving it a set of element nodes indices using the 'indices' data parameter. Here,
 *    internal caribou::topology::Mesh and caribou::topology::Domain instances will
 *    be automatically created over the mechanical state linked through the 'state' data
 *    parameter and the given indices, respectively.
 * 2. Using the CaribouSplineTopology::setDomain method. Here, the topology will be created from
 *    an already existing Domain. The 'indices' data parameter will be automatically filled-in
 *    from this Domain. Note that this Domain instance must exists throughout the lifespan of
 *    the topology, i.e. the pointer to the Domain must remain valid until this component is
 *    destroyed.
 *
 * @tparam Element The element type of this topology. Will be used to determine the structure
 *                 type that will hold the indices of the nodes. It must inherits from
 *                 caribou::geometry::Element.
 */
template <typename Element>
class CaribouSplineTopology : public sofa::core::objectmodel::BaseObject {
public:
    SOFA_CLASS(SOFA_TEMPLATE(CaribouSplineTopology, Element), sofa::core::objectmodel::BaseObject);

    // Type definitions
    using DataTypes = typename SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename Coord::value_type;
    using PointID  = sofa::core::topology::Topology::PointID;


    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<CaribouSplineTopology<Element>, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    template <typename T>
    using Data = sofa::core::objectmodel::Data<T>;

    using mat =sofa::defaulttype::BaseMatrix;

    template<int nRows, int nColumns>
    using Matrix =Eigen::Matrix<Real, nRows, nColumns, Eigen::RowMajor>;


    template<int nRows, int nColumns>
    using Tensor =Eigen::Matrix<Real, nRows, nColumns, Eigen::RowMajor>;

//    using Domain = typename caribou::topology::SplinePatch<Dimesion, PointID>;
    using LocalCoordinates = typename caribou::geometry::Element<Element>::LocalCoordinates;
    using WorldCoordinates = typename caribou::geometry::Element<Element>::WorldCoordinates;
    static constexpr INTEGER_TYPE Dimension = caribou::geometry::traits<Element>::Dimension;
    static constexpr INTEGER_TYPE NumberOfNodes = caribou::geometry::traits<Element>::NumberOfNodesAtCompileTime;
    static constexpr INTEGER_TYPE KnotDimension = 4;
    // SplinePatch typedef
    using SplinePatch = typename caribou::topology::SplinePatch<Dimension, PointID>;

    // Public methods
    CaribouSplineTopology();
    void init() override;

    [[nodiscard]] auto
    getTemplateName() const -> std::string override {
        return templateName(this);
    }
    static auto templateName(const CaribouSplineTopology<Element>* = nullptr) -> std::string {
        return "Unknown";
    }

    /**
     * Construct an element of the domain using the positions vector of the associated mechanical state at rest (rest_position).
     * @param element_id The id of the element in this topology.
     * @return A new Element instance from the domain.
     */
    [[nodiscard]] auto
    element(const UNSIGNED_INTEGER_TYPE & element_id) const noexcept -> Element {
        return this->p_patch->element(element_id);
    }

    /**
     * Get the pointer to the internal domain instance.
     */
    inline auto splinepatch() const noexcept -> const SplinePatch * {return p_patch;}

    /**
     * Get the number of elements contained inside this topology.
     */
    [[nodiscard]] inline auto number_of_elements() const noexcept {
        return p_patch->number_of_elements();
    }

    /**
     * Attach the underlying Domain to the given one and fill in the data attribute 'indices'.
     *
     * When a set of element indices is already there, it is overridden by the indices
     * of the given domain.
     *
     * \note No validation is performed to make sure the element indices are valid
     *       with respect to the position vector of a mechanical state.
     * \note The given domain instance must stay valid throughout the lifespan of this
     *       component.
     */
    void attachSplinePatch(const SplinePatch * p_patch);
    void attachSplinePatch(const SplinePatch & p_patch);

//    /**
//     * Create the underlying Domain by coping the indices from the data attribute 'indices'.
//     *
//     * When a set of element indices is already there, it is overridden by the indices
//     * of the given domain. If a Domain is attached to this topology (see CaribouSplineTopology::attachDomain()), its
//     * reference will be removed.
//     *
//     * Internal Mesh and Domain instances will be created and filled with the position vector
//     * of the mechanical state using the data attribute 'state', and the indices vector using
//     * the data attribute 'indices', respectively.
//     *
//     * \note No validation is performed to make sure the element indices are valid
//     *       with respect to the position vector of the mechanical state.
//     */
//    void initializeFromIndices();

    /**
     * Return true if the SOFA mesh topology is compatible with the type Element of this component.
     *
     * This method can be used when we wish to make sure a SOFA topology (such as a topology container or mesh topology)
     * is compatible with this topology, i.e. it can be used to populate its topology indices vector.
     *
     * This method should be template specialized for every Element type.
     */
    inline
    static auto mesh_is_compatible(const sofa::core::topology::BaseMeshTopology *) -> bool {
        return false;
    }

    /**
     * Return a pointer to the Data member containing the indices vector of the given SOFA topology. For example,
     * passing an HexahedronSetTopologyContainer as the arguement will return its "hexahedra" data.
     *
     * This method should be template specialized for every Element type.
     */
    inline
    static auto get_indices_data_from(const sofa::core::topology::BaseMeshTopology *) -> sofa::core::objectmodel::BaseData * {
        return nullptr;
    }

private:
    // Data members

    /// Position vector of the domain's nodes.
    Data<VecCoord> d_position;

    /// Node indices (w.r.t the position vector) of each elements.
    Data<sofa::type::vector<sofa::type::fixed_array<PointID, NumberOfNodes>>> d_indices;

//    /// Knot ranges
    Data<sofa::type::vector<sofa::type::fixed_array<Real, KnotDimension>>> d_knot_spans;

//    /// Weights
    Data<sofa::type::vector<Real>> d_weights;
    // U direction knot vector
    Data<sofa::type::vector<Real>> d_knot_1;
    // V direction knot vector
    Data<sofa::type::vector<Real>> d_knot_2;

    /// Extraction
//    Data<sofa::type::vector< sofa::type::vector< sofa::type::vector<Real> > > > d_extractions;
//    Data<Matrix<NumberOfNodes, NumberOfNodes> > d_extractions;

//    Data<sofa::type::vector<sofa::type::fixed_array<sofa::type::fixed_array<Real, NumberOfNodes>, NumberOfNodes>>> d_extractions;
//    Data<sofa::type::vector<Eigen::Matrix<Real, NumberOfNodes, NumberOfNodes, Eigen::RowMajor>>> d_extractions;
    /// Pointer to the Domain representing this topology of elements.
//    const Domain * p_domain {nullptr};
    // Pointer to spinepatch
    const SplinePatch * p_patch ;
    const SplinePatch p_patch_n ;
//    std::unique_ptr< caribou::topology::SplinePatch<Dimension, PointID> > p_patch;
};

} // namespace SofaCaribou::topology
