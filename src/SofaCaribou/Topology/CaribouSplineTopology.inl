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
: d_degrees(initData(&d_degrees,
    "degrees",
    "Degrees of each parametric direction."))
, d_position(initData(&d_position,
    "position",
    "Position vector of the patch nodes."))
, d_indices(initData(&d_indices,
    "indices",
    "Node indices (w.r.t the position vector) of each elements."))
, d_knot_spans(initData(&d_knot_spans,
    "knot_spans",
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
{}

template <typename Element>
void CaribouSplineTopology<Element>::attachSplinePatch(const caribou::topology::SplinePatch<Dimension, PointID> * patch) {
    this->p_patch = patch;
    // creating d_indices pointer to write into it.
    using namespace sofa::helper;
    auto degrees = WriteOnlyAccessor<Data<sofa::type::vector<UNSIGNED_INTEGER_TYPE>>>(d_degrees);
    auto positions = WriteOnlyAccessor<Data<VecCoord>>(d_position);
//    auto indices = WriteOnlyAccessor<Data<sofa::type::vector<sofa::type::vector<PointID>>>>(d_indices);
    auto indices = WriteOnlyAccessor<Data<sofa::type::vector<sofa::type::fixed_array<PointID, 9>>>>(d_indices);
    auto knot_spans = WriteOnlyAccessor<Data<sofa::type::vector<sofa::type::fixed_array<Real, KnotDimension>>>>(d_knot_spans);
    auto weights = WriteOnlyAccessor<Data<sofa::type::vector<Real>>>(d_weights);
    auto knot_1 = WriteOnlyAccessor<Data<sofa::type::vector<Real>>>(d_knot_1);
    auto knot_2 = WriteOnlyAccessor<Data<sofa::type::vector<Real>>>(d_knot_2);

    degrees.resize(2);
    degrees[0] = static_cast<UNSIGNED_INTEGER_TYPE>(patch->degree_1());
    degrees[1] = static_cast<UNSIGNED_INTEGER_TYPE>(patch->degree_2());
    const auto nodes_per_element = patch->number_of_nodes_per_elements();
    const auto number_of_elements = patch->number_of_elements();
    indices.resize(number_of_elements);
    knot_spans.resize(number_of_elements);
   // extractions.resize(number_of_elements);
    for (sofa::Index element_id = 0; element_id < number_of_elements; ++element_id) {
        const auto element_indices = patch->element_indices(element_id);
        const auto element_knots = patch->element_knotranges(element_id);
//        indices[element_id].resize(nodes_per_element);
        for (sofa::Index node_id = 0; node_id < nodes_per_element; ++node_id) {
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
void CaribouSplineTopology<Element>::intialise_from_scene(void){
    using namespace sofa::helper;

    auto degrees = ReadAccessor<Data<sofa::type::vector<UNSIGNED_INTEGER_TYPE>>>(d_degrees);

    if (degrees.empty()) {
        msg_warning() << "Degrees are empty " << d_degrees.getName() << "";
//        return;
    }

    // Sanity checks
    auto indices = ReadAccessor<Data<sofa::type::vector<sofa::type::fixed_array<PointID, 9>>>>(d_indices);
//    auto indices = ReadAccessor<Data<sofa::type::vector<sofa::type::vector<PointID>>>>(d_indices);
//    auto indices = ReadAccessor<Data<sofa::type::vector<sofa::type::vector<PointID>>>>(d_indices);

    if (indices.empty()) {
        msg_warning() << "Indices are empty " << d_indices.getName() << "";
//        return;
    }
    else {
        std::cout << "Indices are : " << indices[0] << "\n";
        std::cout << "Indices are : " << indices[1] << "\n";
    }

    auto positions = ReadAccessor<Data<VecCoord>>(d_position);
    if (positions.empty()) {
        msg_warning() << "Positions are not assigned - "<< d_position.getName() << " ";
//        return;
    }

    auto weights = ReadAccessor<Data<sofa::type::vector<Real>>>(d_weights);
    auto knot_spans = ReadAccessor<Data<sofa::type::vector<sofa::type::fixed_array<Real, KnotDimension>>>>(d_knot_spans);
    auto knot_1 = ReadAccessor<Data<sofa::type::vector<Real>>>(d_knot_1);
    auto knot_2 = ReadAccessor<Data<sofa::type::vector<Real>>>(d_knot_2);

    if (weights.empty()) {
        msg_warning() << "Weights are not assigned";
    }
    if (knot_spans.empty()) {
        msg_warning() << "Knot span is not assigned";
    }
    if (knot_1.empty()) {
        msg_warning() << "Knot 1 is not assigned";
    }
    if (knot_2.empty()) {
        msg_warning() << "Knot 2 is not assigned";
    }

    if (positions.empty() or indices.empty() or weights.empty() or
        knot_spans.empty() or knot_1.empty() or knot_2.empty()) {
        msg_warning() << "In big if condition";

    }
    else{
        msg_warning() << "In big Else";
        const auto number_of_nodes = positions.size();
        const auto number_of_weights = weights.size();

        const auto number_of_elements = indices.size();
//        const auto nodes_per_element = indices[0].size();

        const auto knot_1_size = knot_1.size();
        const auto knot_2_size = knot_2.size();


        if (number_of_nodes != number_of_weights){
            msg_warning() << "Nodes and weights are not matching";
        }

        msg_warning() << "Positions size  " << positions.size() << " ";
        msg_warning() << "Indices size    " << indices.size() << " ";
        msg_warning() << "Weights size    " << weights.size() << " ";
        msg_warning() << "Knot spans size " << knot_spans.size() << " ";
        msg_warning() << "Knot 1 size     " << knot_1.size() << " ";
        msg_warning() << "Knot 2 size     " << knot_2.size() << " ";
        msg_warning() << "Nodes per elem  " << 9 << " ";
        Double_Matrix s_initial_positions(positions.size(), Dimension);
        Double_Vector s_weights(weights.size());
        USInt_Matrix s_indices(indices.size(),9);
        Double_Matrix s_knot_ranges(knot_spans.size(),KnotDimension);
        Double_Vector s_knot1(knot_1_size);
        Double_Vector s_knot2(knot_2_size);
        USInt_Vector s_degrees(2);

        // Degrees
        s_degrees << degrees[0], degrees[1];
        // Position and Weights
        for (std::size_t i = 0; i < number_of_nodes; ++i) {
            const auto & n = positions[i];
            // 2 dimensional
            s_initial_positions(i, 0) = n[0];
            s_initial_positions(i, 1) = n[1];
            s_weights(i) = weights[i];
        }

        // Indices and knot spans
        for (std::size_t i = 0; i < number_of_elements; i++){
            const auto & ind = indices[i];
            for (int j = 0; j < static_cast<int>(9); ++j) {
                s_indices(i, j) = ind[j];
            }
            const auto & k_span = knot_spans[i];
            for (int j = 0; j < KnotDimension; ++j) {
                s_knot_ranges(i, j) = k_span[j];
            }
        }

        // Knot 1
        for (std::size_t i = 0; i < knot_1_size; i++){
            s_knot1(i) = knot_1[i];
        }

        // knot 2
        for (std::size_t i = 0; i < knot_2_size; i++){
            s_knot2(i) = knot_2[i];
        }

        std::cout << "================= END IN CST ====================\n";
        std::cout << "Position \n" << s_initial_positions << "\n";
        std::cout << "indices \n" << s_indices << "\n";
        std::cout << "weights \n" << s_weights << "\n";
        std::cout << "spans \n" << s_knot_ranges << "\n";
        std::cout << "knot 1 \n" << s_knot1 << "\n";
        std::cout << "knot 2 \n" << s_knot2 << "\n";
        std::cout << "================= START IN CST ====================\n";
        this->p_patch_n = SplinePatch(s_degrees, s_initial_positions, s_weights, s_indices, s_knot1, s_knot2, s_knot_ranges);
        this->p_patch = new SplinePatch(s_degrees, s_initial_positions, s_weights, s_indices, s_knot1, s_knot2, s_knot_ranges);
        p_patch->print_spline_patch();
        return;

    }




}

template<typename Element>
void CaribouSplineTopology<Element>::init() {
    using namespace sofa::core::objectmodel;

    std::cout << "Hello I am in Caribout Spline Topology\n";
    if (d_indices.isSet()) {
        msg_warning()<< "In indicies if";
        if (p_patch != nullptr) {
            msg_warning()
                    << "Initializing the topology from a set of indices, while a domain has already been attached. The "
                       "latter will be overridden.";
        }

        // If some indices are set, but no mechanical state is linked, let's try to find one in the current context node
        if (not d_indices.getValue().empty() and not d_position.isSet()) {
            auto * context = this->getContext();
            auto state = context->template get<sofa::core::State<DataTypes>>(BaseContext::Local);
            if (state and state->findData("rest_position")) {
                d_position.setParent(state->findData("rest_position"));
                msg_info() << "Automatically found the nodes positions from'" << state->findData("rest_position")->getLinkPath() << "'.";
            }
        }

        this->intialise_from_scene();
    }


}

} // namespace SofaCaribou::topology
