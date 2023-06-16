#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Forcefield/ElasticSplineForcefield.h>
#include <SofaCaribou/Forcefield/CaribouSplineForcefield.inl>
#include <SofaCaribou/Topology/CaribouSplineTopology.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/MechanicalParams.h>
DISABLE_ALL_WARNINGS_END

#include <Caribou/Mechanics/Elasticity/Strain.h>
#ifdef CARIBOU_WITH_OPENMP
#include <omp.h>
#endif

namespace SofaCaribou::forcefield {

template <typename Element>
ElasticSplineForcefield<Element>::ElasticSplineForcefield()
: d_material(initLink(
    "material",
    "Material used to compute the hyperelastic force field."))
, d_enable_multithreading(initData(&d_enable_multithreading,
    false,
    "enable_multithreading",
    "Enable the multithreading computation of the stiffness matrix. Only use this if you have a "
    "very large number of elements, otherwise performance might be worse than single threading."
    "When enabled, use the environment variable OMP_NUM_THREADS=N to use N threads."))
{
}

template <typename Element>
void ElasticSplineForcefield<Element>::init()
{
    using sofa::core::topology::BaseMeshTopology;
    using sofa::core::objectmodel::BaseContext;
    Inherit::init();

    // No material set, try to find one in the current context
    if (not d_material.get()) {
        auto materials = this->getContext()->template getObjects<material::ElasticMaterial<DataTypes>>(BaseContext::Local);
        if (materials.empty()) {
            msg_warning() << "Could not find an hyperelastic material in the current context.";
        } else if (materials.size() > 1) {
            msg_warning() << "Multiple materials were found in the context node. "   <<
                             "Please specify which one should be use by explicitly " <<
                             "setting the material's path in the '" << d_material.getName() << "' parameter.";
        } else {
            d_material.set(materials[0]);
            msg_info() << "Automatically found the material '" << d_material.get()->getPathName() << "'.";
        }
    }

    Real E0  = 10000;
    Real nu0 = 0.3;
    p_C <<  1,    nu0,    0,
          nu0,      1,    0,
            0,      0, (1-nu0)/2;

    p_C = (E0/(1-nu0*nu0)) * p_C;

    // Compute and store the shape functions and their derivatives for every integration points
    initialize_elements();

    // Assemble the initial stiffness matrix
    assemble_stiffness();
}

template<typename Element>
void ElasticSplineForcefield<Element>::addForce(const sofa::core::MechanicalParams *mparams, sofa::core::MultiVecDerivId fId) {
    if (mparams) {
        // Stores the identifier of the x position vector for later use in the stiffness matrix assembly.
        p_X_id = mparams->x();
        std::cout << "p_X_id\n";
        std::cout << p_X_id << "\n";
    }
    std::cout << "\n1st Force function " << fId << " \n";
//    std::cout << fId;
    Inherit::addForce(mparams, fId);
}

template <typename Element>
void ElasticSplineForcefield<Element>::addForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::core::objectmodel::Data<VecDeriv>& d_f,
    const sofa::core::objectmodel::Data<VecCoord>& d_x,
    const sofa::core::objectmodel::Data<VecDeriv>& d_v)
{
    using namespace sofa::core::objectmodel;
    using namespace sofa::helper;
    std::cout << "\n2nd Force function\n";

    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_v);

    if (!this->mstate)
        return;

    const auto material = d_material.get();
    if (!material) {
        return;
    }

    // Update material parameters in case the user changed it
    material->before_update();

    ReadAccessor<Data<VecCoord>> sofa_x = d_x;
    WriteAccessor<Data<VecDeriv>> sofa_f = d_f;

    if (sofa_x.size() != sofa_f.size())
        return;
    const auto nb_nodes = sofa_x.size();
    const auto nb_elements = this->number_of_elements();

    if (nb_nodes == 0 || nb_elements == 0)
        return;

    if (p_elements_quadrature_nodes.size() != nb_elements)
        return;

    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.ref().data()->data(),  nb_nodes, Dimension);
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>> forces  (&(sofa_f[0][0]),  nb_nodes, Dimension);

//    std::cout << "X position\n" << X << "\n";
    std::cout << "Forces\n" << forces << "\n";
    sofa::helper::AdvancedTimer::stepBegin("ElasticSplineForcefield::addForce");

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Fetch the node indices of the element
        auto node_indices = this->topology()->splinepatch()->element_indices(element_id);

        // Fetch the initial and current positions of the element's nodes
        Matrix<NumberOfNodesPerElement, Dimension> current_nodes_position;

        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            current_nodes_position.row(i).noalias() = X.row(node_indices[i]);
        }
        // Compute the nodal forces
        Matrix<NumberOfNodesPerElement, Dimension> nodal_forces;
        nodal_forces.fill(0);
        // Jacobian para to parent
        const auto J2 = p_jacobian_pp[element_id];
        for (GaussNode &gauss_node : p_elements_quadrature_nodes[element_id]) {

            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto & dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Deformation tensor at gauss node
            const Mat22 F = current_nodes_position.transpose()*dN_dx;
            const auto J = F.determinant();
            // Right Cauchy-Green strain tensor at gauss node
            const Mat22 C = F.transpose() * F;

            // Second Piola-Kirchhoff stress tensor at gauss node
            const Mat22 S = material->PK2_stress(J, C);
//            std::cout << "\n-----------------START-------------------\n";
//            std::cout << "dN_dx \n" << dN_dx << "\n";
//            std::cout << "Deformation tensor F \n" << F << "\n";
//            std::cout << "Deformation jacobian J \n" << J << "\n";
//            std::cout << "C-G strain tensor \n" << C << "\n";
//            std::cout << "PK2 stress tensor  \n" << S << "\n";
//            std::cout << "\n------------------END-------------------\n";
            // Elastic forces w.r.t the gauss node applied on each nodes
            for (size_t i = 0; i < NumberOfNodesPerElement; ++i) {
                const auto dx = dN_dx.row(i).transpose();
                const Vector<Dimension> f_ = (detJ * w * J2) * F*S*dx;
                for (size_t j = 0; j < Dimension; ++j) {
                    nodal_forces(i, j) += f_[j];
                }
            }
        }
        std::cout << "\nNodal Force \n" << nodal_forces << "\n";
        for (size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            for (size_t j = 0; j < Dimension; ++j) {
                sofa_f[node_indices[i]][j] -= nodal_forces(i,j);
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("ElasticSplineForcefield::addForce");

    // This is the only I found to detect when a stiffness matrix reassembly is needed for calls to addDForce
    K_is_up_to_date = false;
    eigenvalues_are_up_to_date = false;
}

template <typename Element>
void ElasticSplineForcefield<Element>::addDForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::core::objectmodel::Data<VecDeriv>& d_df,
    const sofa::core::objectmodel::Data<VecDeriv>& d_dx)
{
    using namespace sofa::core::objectmodel;
    std::cout << "\n3rd Force function\n";

    if (not K_is_up_to_date) {
        assemble_stiffness();
    }

    auto kFactor = static_cast<Real> (mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue()));
    sofa::helper::ReadAccessor<Data<VecDeriv>> sofa_dx = d_dx;
    sofa::helper::WriteAccessor<Data<VecDeriv>> sofa_df = d_df;

    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> DX   (&(sofa_dx[0][0]), sofa_dx.size()*3);
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>>       DF   (&(sofa_df[0][0]), sofa_df.size()*3);

    sofa::helper::AdvancedTimer::stepBegin("ElasticSplineForcefield::addDForce");

    for (int k = 0; k < p_K.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<Real>::InnerIterator it(p_K, k); it; ++it) {
            const auto i = it.row();
            const auto j = it.col();
            const auto v = -1 * it.value() * kFactor;
            if (i != j) {
                DF[i] += v*DX[j];
                DF[j] += v*DX[i];
            } else {
                DF[i] += v*DX[i];
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("ElasticSplineForcefield::addDForce");
}

template <typename Element>
void ElasticSplineForcefield<Element>::addKToMatrix(
    sofa::defaulttype::BaseMatrix * matrix,
    SReal kFact, unsigned int & offset)
{
    if (not K_is_up_to_date) {
        assemble_stiffness();
    }
    std::cout << "\n4th Force function\n";

    sofa::helper::AdvancedTimer::stepBegin("ElasticSplineForcefield::addKToMatrix");

    // K is symmetric, so we only stored "one side" of the matrix.
    // But to accelerate the computation, coefficients were not
    // stored only in the upper or lower triangular part, but instead
    // in whatever triangular part (upper or lower) the first node
    // index of the element was. This means that a coefficient (i,j)
    // might be in the lower triangular part, while (k,l) is in the
    // upper triangular part. But no coefficient will be both in the
    // lower AND the upper part.

    for (int k = 0; k < p_K.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<Real>::InnerIterator it(p_K, k); it; ++it) {
            const auto i = it.row();
            const auto j = it.col();
            const auto v = -1 * it.value() * kFact;
            if (i != j) {
                matrix->add(offset+i, offset+j, v);
                matrix->add(offset+j, offset+i, v);
            } else {
                matrix->add(offset+i, offset+i, v);
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("ElasticSplineForcefield::addKToMatrix");
}

template <typename Element>
SReal ElasticSplineForcefield<Element>::getPotentialEnergy (
    const sofa::core::MechanicalParams* mparams,
    const sofa::core::objectmodel::Data<VecCoord>& d_x) const {
    using namespace sofa::core::objectmodel;

    SOFA_UNUSED(mparams);
    std::cout << "\n Potential energy\n";

    if (!this->mstate)
        return 0.;

    const auto material = d_material.get();
    if (!material) {
        return 0;
    }

    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x = d_x;
    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();

    if (sofa_x.size() != sofa_x0.size() )
        return 0.;

    const auto nb_nodes = sofa_x.size();
    const auto nb_elements = this->number_of_elements();

    if (nb_nodes == 0 || nb_elements == 0)
        return 0;

    if (p_elements_quadrature_nodes.size() != nb_elements)
        return 0;

    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.ref().data()->data(),  nb_nodes, Dimension);
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), nb_nodes, Dimension);

    SReal Psi = 0.;

    sofa::helper::AdvancedTimer::stepBegin("ElasticSplineForcefield::getPotentialEnergy");

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        // Fetch the node indices of the element
        auto node_indices = this->topology()->splinepatch()->element_indices(element_id);

        // Fetch the initial and current positions of the element's nodes
        Matrix<NumberOfNodesPerElement, Dimension> initial_nodes_position;
        Matrix<NumberOfNodesPerElement, Dimension> current_nodes_position;

        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            initial_nodes_position.row(i).noalias() = X0.row(node_indices[i]);
            current_nodes_position.row(i).noalias() = X.row(node_indices[i]);
        }

        // Compute the nodal displacement
        Matrix<NumberOfNodesPerElement, Dimension> U {};
        for (size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            const auto u = sofa_x[node_indices[i]] - sofa_x0[node_indices[i]];
            for (size_t j = 0; j < Dimension; ++j) {
                U(i, j) = u[j];
            }
        }

        // Compute the nodal forces

        for (const GaussNode & gauss_node : p_elements_quadrature_nodes[element_id]) {

            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto & detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto & dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto & w = gauss_node.weight;

            // Deformation tensor at gauss node
            const auto & F = caribou::mechanics::elasticity::strain::F(dN_dx, U);
            const auto J = F.determinant();

            // Strain tensor at gauss node
            const Mat22 C = F.transpose() * F;

            // Add the potential energy at gauss node
            Psi += (detJ * w) *  material->strain_energy_density(J, C);
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("ElasticSplineForcefield::getPotentialEnergy");

    return Psi;
}

template <typename Element>
void ElasticSplineForcefield<Element>::initialize_elements()
{
    using namespace sofa::core::objectmodel;

    sofa::helper::AdvancedTimer::stepBegin("ElasticSplineForcefield::initialize_elements");

    if (!this->mstate)
        return;
    // Resize the container of elements'quadrature nodes
    const auto nb_elements = this->number_of_elements();
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
        const auto initial_element = this->topology()->element(element_id);
        // Fill in the Gauss integration nodes for this element
        p_elements_quadrature_nodes[element_id] = get_gauss_nodes(element_id, initial_element);
        p_jacobian_pp[element_id] = initial_element.jacobian_papa();
    }
    // Compute the volume
    Real v = 0.;
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        for (const auto & gauss_node : gauss_nodes_of(element_id)) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            v += detJ*w;
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("ElasticSplineForcefield::initialize_elements");
}

template <typename Element>
void ElasticSplineForcefield<Element>::assemble_stiffness()
{
    assemble_stiffness(*this->mstate->read (p_X_id.getId(this->mstate)));
}

template<typename Element>
void ElasticSplineForcefield<Element>::assemble_stiffness(const sofa::core::objectmodel::Data<VecCoord> & x) {
    using namespace sofa::core::objectmodel;

    const sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x= x;
    const auto nb_nodes = sofa_x.size();
    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.ref().data()->data(),  nb_nodes, Dimension);

    assemble_stiffness(X);
}

template<typename Element>
template<typename Derived>
void ElasticSplineForcefield<Element>::assemble_stiffness(const Eigen::MatrixBase<Derived> & x) {
    const auto material = d_material.get();

    [[maybe_unused]]
    const auto enable_multithreading = d_enable_multithreading.getValue();
    if (!material) {
        return;
    }

    // Update material parameters in case the user changed it
    material->before_update();

//    static const auto Id = Mat22::Identity();
    const auto nb_elements = this->number_of_elements();
    const auto nb_nodes = x.rows();
    const auto nDofs = nb_nodes*Dimension;
    p_K.resize(nDofs, nDofs);

    ///< Triplets are used to store matrix entries before the call to 'compress'.
    /// Duplicates entries are summed up.
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(nDofs*24*2);

    sofa::helper::AdvancedTimer::stepBegin("ElasticSplineForcefield::update_stiffness");
#pragma omp parallel for if (enable_multithreading)
    for (int element_id = 0; element_id < static_cast<int>(nb_elements); ++element_id) {
        // Fetch the node indices of the element
        auto node_indices = this->topology()->splinepatch()->element_indices(element_id);

        // Fetch the current positions of the element's nodes
        Matrix<NumberOfNodesPerElement, Dimension> current_nodes_position;

        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            current_nodes_position.row(i).noalias() = x.row(node_indices[i]).template cast<Real>();
        }

        using Stiffness = Eigen::Matrix<FLOATING_POINT_TYPE, NumberOfNodesPerElement*Dimension, NumberOfNodesPerElement*Dimension, Eigen::RowMajor>;
        Stiffness Ke = Stiffness::Zero();

        // Jacobian para to parent
        const auto J2 = p_jacobian_pp[element_id];
        int count = 1;
        for (const auto & gauss_node : gauss_nodes_of(element_id)) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto detJ = gauss_node.jacobian_determinant;

            // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
            const auto dN_dx = gauss_node.dN_dx;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            Matrix<3, NumberOfNodesPerElement*Dimension> B;
            B.setZero();
            for (int i = 0; i<NumberOfNodesPerElement; i++){
                B(0, i*Dimension) = dN_dx(i, 0);
                B(1, i*Dimension + 1) = dN_dx(i, 1);
                B(2, i*Dimension) = dN_dx(i, 1);
                B(2, i*Dimension + 1) = dN_dx(i, 0);
            }

//            std::cout << "\n --------------------Start-------------------------------\n";
//            std::cout << "Gauss point : " << count << "\n";
//            std::cout << "B Matrix \n" << B << "\n";
//            std::cout << "B Matrix size\n" << dN_dx << "\n";
//            std::cout << "Elasticity Matrix  : \n" << p_C << "\n";
//            std::cout << "Normal jacob : " << detJ << "\n";
//            std::cout << "IGA jacob : " << J2 << "\n";
//            std::cout << "Weight : " << w << "\n";
            Ke = Ke + B.transpose() * p_C * B * detJ * J2 * w;
            count = count + 1;
        }
//        std::cout << "\n Element stiffness \n" << Ke << "\n";


#pragma omp critical
        for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
            // Node index of the ith node in the global stiffness matrix
            const auto x = static_cast<int>(node_indices[i]*Dimension);
            for (int m = 0; m < Dimension; ++m) {
                for (int n = m; n < Dimension; ++n) {
                    triplets.emplace_back(x+m, x+n, Ke(i*Dimension+m,i*Dimension+n));
                }
            }

            for (std::size_t j = i+1; j < NumberOfNodesPerElement; ++j) {
                // Node index of the jth node in the global stiffness matrix
                const auto y = static_cast<int>(node_indices[j]*Dimension);
                for (int m = 0; m < Dimension; ++m) {
                    for (int n = 0; n < Dimension; ++n) {
                        triplets.emplace_back(x+m, y+n, Ke(i*Dimension+m,j*Dimension+n));
                    }
                }
            }
        }
    }
    p_K.setFromTriplets(triplets.begin(), triplets.end());
    sofa::helper::AdvancedTimer::stepEnd("ElasticSplineForcefield::update_stiffness");

    K_is_up_to_date = true;
    eigenvalues_are_up_to_date = false;
}

template <typename Element>
auto ElasticSplineForcefield<Element>::get_gauss_nodes(const std::size_t & /*element_id*/, const Element & element) const -> GaussContainer {
    GaussContainer gauss_nodes {};
    if constexpr (NumberOfGaussNodesPerElement == caribou::Dynamic) {
        gauss_nodes.resize(element.number_of_gauss_nodes());
    }
    const auto nb_of_gauss_nodes = gauss_nodes.size();
    for (std::size_t gauss_node_id = 0; gauss_node_id < nb_of_gauss_nodes; ++gauss_node_id) {
        const auto & g = element.gauss_node(gauss_node_id);

        const auto J = element.jacobian(g.position);
        const Mat22 Jinv = J.inverse();
        const auto detJ = std::abs(J.determinant());

        // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
        const Matrix<NumberOfNodesPerElement, Dimension> dN_dx =
            (Jinv.transpose() * element.dL(g.position).transpose()).transpose();

        GaussNode & gauss_node = gauss_nodes[gauss_node_id];
        gauss_node.weight               = g.weight;
        gauss_node.jacobian_determinant = detJ;
        gauss_node.dN_dx                = dN_dx;
    }
    return gauss_nodes;
}

template <typename Element>
auto ElasticSplineForcefield<Element>::eigenvalues() -> const Vector<Eigen::Dynamic> & {
    if (not eigenvalues_are_up_to_date) {
#ifdef EIGEN_USE_LAPACKE
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> k (K());
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> eigensolver(k, Eigen::EigenvaluesOnly);
#else
        Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<Real>> eigensolver(K(), Eigen::EigenvaluesOnly);
#endif
        if (eigensolver.info() != Eigen::Success) {
            msg_error() << "Unable to find the eigen values of K.";
        }

        p_eigenvalues = eigensolver.eigenvalues();
        eigenvalues_are_up_to_date = true;
    }

    return p_eigenvalues;
}

template <typename Element>
auto ElasticSplineForcefield<Element>::cond() -> Real {
    const auto & values = eigenvalues();
    const auto min = values.minCoeff();
    const auto max = values.maxCoeff();

    return min/max;
}

} // namespace SofaCaribou::forcefield
