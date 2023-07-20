#include <SofaCaribou/config.h>
#include <SofaCaribou/Forcefield/TractionSplineForcefield.h>
#include <SofaCaribou/Forcefield/CaribouSplineForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/helper/AdvancedTimer.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 210600)
#include <sofa/helper/ScopedAdvancedTimer.h>
#endif
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION <= 201299)
namespace sofa::helper::visual { using DrawTool = sofa::core::visual::DrawTool; }
#endif

namespace SofaCaribou::forcefield {

template <typename Element>
TractionSplineForcefield<Element>::TractionSplineForcefield()
    : d_traction(initData(&d_traction,
               "traction",
               "Loading on a boundary with directions."))
    , d_boundary(initData(&d_boundary,
               "boundary",
               "Where traction is being applied."))
{
}

template <typename Element>
void TractionSplineForcefield<Element>::init()
{
    Inherit::init();
    std::cout << "\nHello I am in Traction Spline\n";
    // Initializing the elements;
    initialize_elements();
    apply_load();
}



template <typename Element>
void TractionSplineForcefield<Element>::apply_load()
{
    std::cout << "\nHello I am in Traction Spline : Load\n";
    sofa::helper::WriteAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;

    const auto rest_positions = this->getMState()->readRestPositions();
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0 (rest_positions.ref().data()->data(),  rest_positions.size(), Dimension);

    const auto nb_nodes = X0.rows();

    std::cout << "Number of nodes : " << nb_nodes << "\n";

    nodal_forces.resize(nb_nodes);

//    for (int i=0; i < nb_nodes; i++){
//        std::cout << "Printing nodal_forces : " << " i : " << nodal_forces[i]<< "\n";
//    }


    const int & boundary = d_boundary.getValue();

    size_t nb_elements;
    if (boundary == 1 || boundary == 3){
        nb_elements = this->topology()->splinepatch()->number_of_elems_on_boundary(1);
    }
    else if (boundary == 1 || boundary == 3){
        nb_elements = this->topology()->splinepatch()->number_of_elems_on_boundary(2);
    }
    else{
        throw std::invalid_argument("Please provide valid argument : 1-Down, 2-Left, 3-Up, 4-Right");
    }

    Deriv traction = d_traction.getValue();

//    std::cout << " +++++++++++++ Traction Force +++++++++++++ ";

    std::cout << "Boundary is : " << boundary << "\n";
    std::cout << "traction is : " << traction << "\n";

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Fetch the node indices of the element
        auto node_indices = this->topology()->splinepatch()->element_boundary_nodes(boundary, element_id);

        // Fetch the initial positions of the element's nodes
        Matrix<NumberOfNodesPerElement, Dimension> initial_nodes_position;

        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            initial_nodes_position.row(i).noalias() = X0.row(node_indices[i]);
        }
        Real J2 = p_jacobian_pp[element_id];
//        int count = 0;
        // Integration of the traction increment over the element.
        for (GaussNode &gauss_node : p_elements_quadrature_nodes[element_id]) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Shape values at each nodes evaluated at the gauss point position
            const auto N = gauss_node.N;

            // Traction evaluated at the gauss point position
            const auto F = traction * w * J2 * detJ;

//            std::cout << "------------ Gauss Point : " << count << " -- ";
//            std::cout << "Weight : " << w << "\n";
//            std::cout << "J1 : " << detJ << "\n";
//            std::cout << "J2 : " << J2 << "\n";
//            std::cout << "Basis : \n" << N << "\n";

            // Tractive forces w.r.t the gauss node applied on each nodes
            for (size_t i = 0; i < NumberOfNodesPerElement; ++i) {
                nodal_forces[node_indices[i]] += F*N[i];
            }

        }
    }

//    for (int i=0; i < nb_nodes; i++){
//        std::cout << "Printing nodal_forces : " << i << " : " << nodal_forces[i]<< "\n";
//    }
}

template <typename Element>
void TractionSplineForcefield<Element>::addForce(const sofa::core::MechanicalParams* mparams, Data<VecDeriv>& d_f, const Data<VecCoord>& d_x, const Data<VecDeriv>& /*d_v*/)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_x);
    sofa::helper::AdvancedTimer::stepBegin("TractionSplineForce::addForce");
    sofa::helper::ReadAccessor<Data<VecDeriv>> nodal_forces = d_nodal_forces;
    sofa::helper::WriteAccessor<Data<VecDeriv>> f = d_f;

//    std::cout << "\n============================ 1 ADD Forece ============================\n";
//    for (size_t i=0; i < f.size(); i++){
//        std::cout << "Printing nodal_forces : " << i << " : " << f[i] << "\n";
//    }
//    std::cout << "\n============================ 1 ADD Forece ============================\n";

    for (size_t i = 0; i < f.size(); ++i)
        f[i] += nodal_forces[i];

//    std::cout << "\n============================ 2 ADD Forece ============================\n";
//    for (size_t i=0; i < f.size(); i++){
//        std::cout << "Printing nodal_forces : " << i << " : " << f[i] << "\n";
//    }
//    std::cout << "\n============================ 2 ADD Forece ============================\n";
    sofa::helper::AdvancedTimer::stepEnd("TractionSplineForce::addForce");
}

template<typename Element>
void TractionSplineForcefield<Element>::initialize_elements() {
    using namespace sofa::core::objectmodel;
    std::cout << "\nHello I am in Traction Spline : Init Elements\n";
    sofa::helper::ScopedAdvancedTimer _t_ ("TractionSplineForcefield::initialize_elements");

    if (!this->mstate)
        return;

    const int & boundary = d_boundary.getValue();

    size_t nb_elements;
    if (boundary == 1 || boundary == 3){
        nb_elements = this->topology()->splinepatch()->number_of_elems_on_boundary(1);
    }
    else if (boundary == 1 || boundary == 3){
        nb_elements = this->topology()->splinepatch()->number_of_elems_on_boundary(2);
    }
    else{
        throw std::invalid_argument("Please provide valid argument : 1-Down, 2-Left, 3-Up, 4-Right");
    }

    // Resize the container of elements'quadrature nodes
    if (p_elements_quadrature_nodes.size() != nb_elements) {
        p_elements_quadrature_nodes.resize(nb_elements);
        p_jacobian_pp.resize(nb_elements);
    }

    // Translate the Sofa's mechanical state vector to Eigen vector type
    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), sofa_x0.size(), Dimension);

    // Loop on each element and compute the shape functions and their derivatives for every of their integration points
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Get an Element instance from the Domain
        const auto initial_element = this->topology()->splinepatch()->boundary_element(boundary, element_id);

        // Fill in the Gauss integration nodes for this element
        p_elements_quadrature_nodes[element_id] = get_gauss_nodes(element_id, initial_element);
        p_jacobian_pp[element_id] = initial_element.jacobian_papa();
    }

}

template<typename Element>
auto TractionSplineForcefield<Element>::get_gauss_nodes(const size_t & /*element_id*/, const boundry_Element &element) const -> GaussContainer {
    GaussContainer gauss_nodes {};
    if constexpr (NumberOfGaussNodesPerElement == caribou::Dynamic) {
        gauss_nodes.resize(element.number_of_gauss_nodes());
    }

    const auto nb_of_gauss_nodes = gauss_nodes.size();
    for (std::size_t gauss_node_id = 0; gauss_node_id < nb_of_gauss_nodes; ++gauss_node_id) {
        const auto & g = element.gauss_node(gauss_node_id);

        const auto J = element.jacobian(g.position);
        Real detJ;
        if constexpr (boundry_Element::Dimension == 2 || boundry_Element::Dimension == 3) {
            detJ = std::abs(J.norm());
        } else {
            detJ = std::abs(J.determinant());
        }

        // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
        const Vector<NumberOfNodesPerElement> N = element.L(g.position);


        GaussNode & gauss_node = gauss_nodes[gauss_node_id];
        gauss_node.weight               = g.weight;
        gauss_node.jacobian_determinant = detJ;
        gauss_node.N                    = N;
    }

    return gauss_nodes;
}

} // namespace SofaCaribou::forcefield
