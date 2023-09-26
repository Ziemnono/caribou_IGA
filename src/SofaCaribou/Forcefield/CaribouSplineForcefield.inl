#pragma once

#include <SofaCaribou/Forcefield/CaribouSplineForcefield.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/visual/VisualParams.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
namespace sofa::type {
using RGBAColor = ::sofa::helper::types::RGBAColor;
template <typename Real>
using TBoundingBox = ::sofa::defaulttype::TBoundingBox<Real> ;
}
#endif

namespace SofaCaribou::forcefield {

template<typename Element>
CaribouSplineForcefield<Element>::CaribouSplineForcefield()
        : d_topology_container(initLink(
        "topology",
        "Topology container containing the elements on which this forcefield will be applied.")),
          d_drawScale(initData(&d_drawScale,
                               0.85,
                               "draw_scale",
                               "Scaling factor for the drawing of elements (between 0 and 1). The factor allows to shrink "
                               "the element relative to its center point when drawing it.")) {}

template<typename Element>
void CaribouSplineForcefield<Element>::init() {
    using sofa::core::topology::BaseMeshTopology;
    using sofa::core::objectmodel::BaseContext;
    using CaribouSplineTopology = SofaCaribou::topology::CaribouSplineTopology<Element>;

    Inherit1::init();

    auto *context = this->getContext();

    if (!this->mstate) {
        msg_warning() << "No mechanical object found in the current context node. The data parameter "
                      << "'" << this->mstate.getName() << "' can be use to set the path to a mechanical "
                      << "object having a template of '" << DataTypes::Name() << "'";
        return;
    }

    // If not topology is specified, try to find one automatically in the current context
    if (not d_topology_container.get()) {
        // No topology specified. Try to find one suitable.
        auto caribou_containers = context->template getObjects<CaribouSplineTopology>(
                BaseContext::Local);
        auto sofa_containers = context->template getObjects<BaseMeshTopology>(BaseContext::Local);
        std::vector<BaseMeshTopology *> sofa_compatible_containers;
        for (auto container : sofa_containers) {
            if (CaribouSplineTopology::mesh_is_compatible(container)) {
                sofa_compatible_containers.push_back(container);
            }
        }
        if (caribou_containers.empty() and sofa_compatible_containers.empty()) {
            msg_warning() << "Could not find a topology container in the current context. "
                          << "Please add a compatible one in the current context or set the "
                          << "container's path using the '" << d_topology_container.getName()
                          << "' data parameter.";
        } else {
            if (caribou_containers.size() + sofa_compatible_containers.size() > 1) {
                msg_warning() << "Multiple topologies were found in the context node. "
                              << "Please specify which one contains the elements on "
                              << "which this force field will be applied "
                              << "by explicitly setting the container's path in the  '"
                              << d_topology_container.getName() << "' data parameter.";
            } else {
                // Prefer caribou's containers first
                if (not caribou_containers.empty()) {
                    d_topology_container.set(caribou_containers[0]);
                } else {
                    d_topology_container.set(sofa_compatible_containers[0]);
                }
            }

            msg_info() << "Automatically found the topology '" << d_topology_container.get()->getPathName()
                       << "'.";
        }
    }

    // Create a caribou internal Domain over the topology
    if (d_topology_container.get()) {
        auto sofa_topology = dynamic_cast<BaseMeshTopology *>(d_topology_container.get());
        auto caribou_topology = dynamic_cast<CaribouSplineTopology *>(d_topology_container.get());
        if (sofa_topology) {
            // Initialize a new caribou topology from the SOFA topology
//            p_topology = sofa::core::objectmodel::New<CaribouSplineTopology>();
//            p_topology->findData("indices")->setParent(CaribouSplineTopology::get_indices_data_from(sofa_topology));
//            p_topology->findData("position")->setParent(this->getMState()->findData("position"));
//            p_topology->init();
            msg_warning() << "No caribou topology is exist. \n";
        } else {
            // A Caribou topology already exists in the scene
            p_topology = caribou_topology;
        }

        if (number_of_elements() == 0) {
            msg_warning() << "No element found in the topology '" << d_topology_container.get()->getPathName() << "'";
        }
    }
}

template<typename Element>
void CaribouSplineForcefield<Element>::computeBBox(const sofa::core::ExecParams *, bool onlyVisible) {
    using namespace sofa::core::objectmodel;

    if (!onlyVisible) return;
    if (!this->mstate) return;

    sofa::helper::ReadAccessor<Data < VecCoord>>
    x = this->mstate->read(sofa::core::VecCoordId::position());

    static const Real max_real = std::numeric_limits<Real>::max();
    static const Real min_real = std::numeric_limits<Real>::lowest();
    Real maxBBox[3] = {min_real, min_real, min_real};
    Real minBBox[3] = {max_real, max_real, max_real};
    for (size_t i = 0; i < x.size(); i++) {
        for (int c = 0; c < 3; c++) {
            if (x[i][c] > maxBBox[c]) maxBBox[c] = static_cast<Real>(x[i][c]);
            else if (x[i][c] < minBBox[c]) minBBox[c] = static_cast<Real>(x[i][c]);
        }
    }

    this->f_bbox.setValue(sofa::type::TBoundingBox<Real>(minBBox, maxBBox));
}

template<typename Element>
template<typename Derived>
auto CaribouSplineForcefield<Element>::canCreate(Derived *o,
                                           sofa::core::objectmodel::BaseContext *context,
                                           sofa::core::objectmodel::BaseObjectDescription *arg) -> bool {
    using namespace sofa::core::objectmodel;
    using CaribouSplineTopology = SofaCaribou::topology::CaribouSplineTopology<Element>;

    std::string requested_element_type = arg->getAttribute( "template", "");
    std::string this_element_type = Derived::templateName(o);

    // to lower
    std::string requested_element_type_lower = requested_element_type, this_element_type_lower = this_element_type;
    std::transform(requested_element_type.begin(), requested_element_type.end(), requested_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });
    std::transform(this_element_type.begin(), this_element_type.end(), this_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });

    if (requested_element_type_lower == this_element_type_lower) {
        return Inherit1::canCreate(o, context, arg);
    }

    if (not requested_element_type.empty()) {
        arg->logError("Requested element type is not '"+this_element_type+"'.");
        return false;
    }

    std::string topology_path = arg->getAttribute("topology", "");
    if (not topology_path.empty()) {
        topology_path = topology_path.substr(1); // removes the "@"
        // Make sure the specified topology has elements of type Element
        auto topology = context->get<BaseObject>(topology_path);
        auto caribou_topology = dynamic_cast<CaribouSplineTopology *>(topology);
        auto sofa_topology = dynamic_cast<sofa::core::topology::BaseMeshTopology *>(topology);
        if (not caribou_topology and (not sofa_topology or not CaribouSplineTopology::mesh_is_compatible(sofa_topology))) {
            arg->logError("Cannot deduce the element type from the specified mesh topology '" + topology_path + "'.");
            return false;
        }
    } else {
        // Try to find a compatible topology in the current context
        BaseObject * topology = nullptr;
        auto objects = context->getObjects<BaseObject>(BaseContext::SearchDirection::Local);
        for (auto * object : objects) {
            auto caribou_topology = dynamic_cast<const CaribouSplineTopology *>(object);
            auto sofa_topology = dynamic_cast<const sofa::core::topology::BaseMeshTopology *>(object);
            if (caribou_topology  or (sofa_topology and CaribouSplineTopology::mesh_is_compatible(sofa_topology))) {
                topology = object;
                break;
            }
        }
        if (not topology) {
            arg->logError("Cannot find a topology in the current context from which the template '"+this_element_type+"' can be deduced.");
            return false;
        }

        if (Inherit1::canCreate(o, context, arg)) {
            arg->setAttribute("topology", "@" + topology->getPathName());
            return true;
        } else {
            return false;
        }
    }

    return Inherit1::canCreate(o, context, arg);
}

template<typename Element>
void CaribouSplineForcefield<Element>::draw(const sofa::core::visual::VisualParams *vparams) {
    using Color = sofa::type::RGBAColor;
    static constexpr auto CanonicalDimension = caribou::geometry::traits<Element>::CanonicalDimension;

    using LocalCoordinate = typename Element::LocalCoordinates;
    using Scalar = typename Element::Scalar;

    [[maybe_unused]]
    static const unsigned long long kelly_colors_hex[] =
    {
        0xFFFFB300, // Vivid Yellow
        0xFF803E75, // Strong Purple
        0xFFFF6800, // Vivid Orange
        0xFFA6BDD7, // Very Light Blue
        0xFFC10020, // Vivid Red
        0xFFCEA262, // Grayish Yellow
        0xFF817066, // Medium Gray

        // The following don't work well for people with defective color vision
        0xFF007D34, // Vivid Green
        0xFFF6768E, // Strong Purplish Pink
        0xFF00538A, // Strong Blue
        0xFFFF7A5C, // Strong Yellowish Pink
        0xFF53377A, // Strong Violet
        0xFFFF8E00, // Vivid Orange Yellow
        0xFFB32851, // Strong Purplish Red
        0xFFF4C800, // Vivid Greenish Yellow
        0xFF7F180D, // Strong Reddish Brown
        0xFF93AA00, // Vivid Yellowish Green
        0xFF593315, // Deep Yellowish Brown
        0xFFF13A13, // Vivid Reddish Orange
        0xFF232C16, // Dark Olive Green
    };

    if (!vparams->displayFlags().getShowForceFields())
        return;

    const auto nb_elements = number_of_elements();

    if (nb_elements == 0)
        return;

    vparams->drawTool()->saveLastState();
    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);
    vparams->drawTool()->disableLighting();

//    const double & scale = d_drawScale.getValue();
    const VecCoord& sofa_x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.data()->data(),  sofa_x.size(), Dimension);

    // Subdivision of an element

    const int ndivs = 2; // divisions in a row.
    const int ncells = ndivs * ndivs; // Total number of cells.
    const int tpnts = (ndivs+1) * (ndivs+1);

    int tcells = ncells * nb_elements; // Total number of cells in whole geometry.

    // For 2D and 3D, we draw triangles.
    // We create groups of faces to draw. Each group will be assigned with a color.
    // For 2D, there is only one group (one "visible" face per element...).
    // For 3D, each face of an element is assigned to a group, hence will have a different color.
    std::vector<std::vector<sofa::type::Vector3>> triangle_groups;
    if constexpr(CanonicalDimension == 2) {
        triangle_groups.resize(1);
        auto nnodes = 10; // 4 Nodes per every quad elements.
        triangle_groups[0].reserve(tcells*nnodes); // Numb
    } else {
        std::cout << "ONLY SURFACES CAN BE VISUALISED\n";
    }

    // ------------------ NURBS HELPERS START --------------------

    const Scalar limit = -0.999999;
    const Scalar increment = -2*limit/ndivs;
    std::vector<LocalCoordinate> grid_points;
    grid_points.resize(tpnts);
    int count = 0;
    for( int vi = 0; vi <= ndivs; vi++){
        for( int ui = 0; ui <= ndivs; ui++){
            grid_points[count] = LocalCoordinate(limit + increment*ui, limit + increment*vi);
            count++;
        }
    }
    // Storing NURBS quad node connection data. This will be used in triangularizing the grid.
    Eigen::Matrix<UNSIGNED_INTEGER_TYPE, -1, 4> tn_inds;
    tn_inds.resize(ncells, 4);
    int val;
    count = 0;
    for (int cell = 0; cell < ncells; ++cell) {
        if ((cell != 0) && ((cell)%ndivs == 0)){ count++; }
        val = cell+count;
        tn_inds.row(cell) << val, val+1, val+ndivs+1, val+ndivs+2;
    }

    // ------------------ NURBS HELPERS END ----------------------

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        // Fetch the node indices of the element
        auto node_indices = this->topology()->splinepatch()->element_indices(element_id);
        Matrix<NumberOfNodesPerElement, Dimension> element_nodes_position;

        if constexpr(NumberOfNodesPerElement == caribou::Dynamic) {
            element_nodes_position.resize(node_indices.rows(), Dimension);
        }

        // Fetch the initial positions of the element's nodes
        for (Eigen::Index node_id = 0; node_id < element_nodes_position.rows(); ++node_id) {
            element_nodes_position.row(node_id) = X.row(node_indices[node_id]);
        }

        // Element NURBS data
        auto degrees = this->topology()->splinepatch()->get_degrees();               // degrees
        auto knot1 = this->topology()->splinepatch()->knot_1();                      // knot1
        auto knot2 = this->topology()->splinepatch()->knot_2();                      // knot2
        auto weights = this->topology()->splinepatch()->weights(node_indices);       // weights
        auto knot_span = this->topology()->splinepatch()->knot_range(element_id);    // knot_span

        // Create an Element instance from the node positions
        const Element scaled_element = Element(degrees, element_nodes_position, knot1, knot2, weights, knot_span); // changed e to scaled_element.

        // 2D and 3D elements: we draw triangles
        if constexpr(CanonicalDimension == 2) {
            // 2D elements: we triangulate only one face

            for (int cell = 0; cell < ncells; ++cell) {

                auto triangle_1 = std::vector {scaled_element.world_coordinates(grid_points[tn_inds(cell, 0)]),
                                               scaled_element.world_coordinates(grid_points[tn_inds(cell, 1)]),
                                               scaled_element.world_coordinates(grid_points[tn_inds(cell, 2)])};

                auto triangle_2 = std::vector {scaled_element.world_coordinates(grid_points[tn_inds(cell, 1)]),
                                               scaled_element.world_coordinates(grid_points[tn_inds(cell, 2)]),
                                               scaled_element.world_coordinates(grid_points[tn_inds(cell, 3)])};
                for (const auto &n : triangle_1) {
                    triangle_groups[0].emplace_back(n[0], n[1], 0);
                }
                for (const auto &n : triangle_2) {
                    triangle_groups[0].emplace_back(n[0], n[1], 0);
                }
            }
        }
    }

    for (std::size_t group_id = 0; group_id < triangle_groups.size(); ++group_id) {
        const auto & triangles = triangle_groups[group_id];
        const auto &hex_color = kelly_colors_hex[group_id % 20];
        const Color group_color(
            static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(16))) / 255.),
            static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(8))) / 255.),
            static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(0))) / 255.),
            static_cast<float> (1)
        );
        vparams->drawTool()->drawTriangles(triangles, group_color);
    }

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

} // namespace SofaCaribou::forcefield
