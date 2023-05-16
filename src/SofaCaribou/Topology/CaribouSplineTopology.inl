#pragma once

#include <SofaCaribou/Topology/CaribouSplineTopology.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201299)
namespace sofa { using Index = unsigned int; }
#endif

namespace SofaCaribou::topology {

template <typename Element>
CaribouSplineTopology<Element>::CaribouSplineTopology ()
: d_position(initData(&d_position,
    "position",
    "Position vector of the patch nodes."))
, d_indices(initData(&d_indices,
    "indices",
    "Node indices (w.r.t the position vector) of each elements."))
, d_knot_spans(initData(&d_knot_spans,
    "knots",
    "Knot span of each elements"))
, d_weights(initData(&d_weights,
    "weights",
    "weight vector of the patch nodes."))
, d_knot_1(initData(&d_knot_1,
    "knot_1",
    "weight vector of the patch nodes."))
, d_knot_2(initData(&d_knot_2,
    "knot_2",
    "weight vector of the patch nodes."))
//, d_extractions(initData(&d_extractions,
//    "extractions",
//    "vector of extraction_matrix of each element"))
{}

template <typename Element>
void CaribouSplineTopology<Element>::attachSplinePatch(const caribou::topology::SplinePatch<Dimension, PointID> * patch) {
    this->p_patch = patch;
    // creating d_indices pointer to write into it.
    using namespace sofa::helper;
    auto positions = WriteOnlyAccessor<Data<VecCoord>>(d_position);
    auto indices = WriteOnlyAccessor<Data<sofa::type::vector<sofa::type::fixed_array<PointID, NumberOfNodes>>>>(d_indices);
    auto knot_spans = WriteOnlyAccessor<Data<sofa::type::vector<sofa::type::fixed_array<Real, KnotDimension>>>>(d_knot_spans);
    auto weights = WriteOnlyAccessor<Data<sofa::type::vector<Real>>>(d_weights);
    auto knot_1 = WriteOnlyAccessor<Data<sofa::type::vector<Real>>>(d_knot_1);
    auto knot_2 = WriteOnlyAccessor<Data<sofa::type::vector<Real>>>(d_knot_2);

    const auto number_of_elements = patch->number_of_elements();
    indices.resize(number_of_elements);
    knot_spans.resize(number_of_elements);
   // extractions.resize(number_of_elements);
    for (sofa::Index element_id = 0; element_id < number_of_elements; ++element_id) {
        const auto element_indices = patch->element_indices(element_id);
        const auto element_knots = patch->element_knotranges(element_id);

        for (sofa::Index node_id = 0; node_id < NumberOfNodes; ++node_id) {
            indices[element_id][node_id] = static_cast<PointID>(element_indices[node_id]);
        }

        for (sofa::Index node_id = 0; node_id < KnotDimension; ++node_id) {
            knot_spans[element_id][node_id] = static_cast<Real>(element_knots[node_id]);
        }
    }

    const auto patch_knot_1 = patch->knot_1();
    const auto patch_knot_2 = patch->knot_2();

    knot_1.resize(patch_knot_1.size());
    knot_2.resize(patch_knot_2.size());

    for (sofa::Index i = 0; i < patch_knot_1.size(); ++i) {
        knot_1[i] = static_cast<Real>(patch_knot_1(i));
    }

    for (sofa::Index i = 0; i < patch_knot_2.size(); ++i) {
        knot_2[i] = static_cast<Real>(patch_knot_2(i));
    }


    const auto total_nodes = patch->number_of_nodes();
    weights.resize(total_nodes);
    positions.resize(total_nodes);
    for (sofa::Index node_id = 0; node_id < total_nodes; ++node_id) {
        weights[node_id] = static_cast<Real>(patch->weight(node_id)); // Weights
        // Assigning node data.
        const auto node = patch->position(node_id);
        for (sofa::Index i = 0; i < 2; i++) {
            positions[node_id][i] = static_cast<Real>(node[i]);
        }
    }
}


template<typename Element>
void CaribouSplineTopology<Element>::init() {
    using namespace sofa::core::objectmodel;

//    if (d_indices.isSet()) {
//        if (p_domain != nullptr) {
//            msg_warning()
//                    << "Initializing the topology from a set of indices, while a domain has already been attached. The "
//                       "latter will be overridden.";
//        }

//        // If some indices are set, but no mechanical state is linked, let's try to find one in the current context node
//        if (not d_indices.getValue().empty() and not d_position.isSet()) {
//            auto * context = this->getContext();
//            auto state = context->template get<sofa::core::State<DataTypes>>(BaseContext::Local);
//            if (state and state->findData("rest_position")) {
//                d_position.setParent(state->findData("rest_position"));
//                msg_info() << "Automatically found the nodes positions from'" << state->findData("rest_position")->getLinkPath() << "'.";
//            }
//        }

//        this->initializeFromIndices();
//    }


}

} // namespace SofaCaribou::topology
