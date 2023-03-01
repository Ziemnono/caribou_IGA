#pragma once

#include <Caribou/config.h>

#include <string>
#include <memory>
#include <Eigen/Eigen>

namespace caribou::topology {

    /**
     * A mesh is a discrete view of a space by the mean of nodes. It can contain one or many domains, which are
     * responsible for maintaining the topological relation between the nodes of a subspace of the mesh.
     */
    struct BaseSplinePatch {
        /** Destructor */
        virtual ~BaseSplinePatch() = default;

        /*!
         * Get the dimension of the mesh, ie, the number of coordinates of a point.
         */
        [[nodiscard]] virtual auto dimension() const -> UNSIGNED_INTEGER_TYPE = 0;

        /*!
         * Get the number of nodes of the mesh.
         */
        [[nodiscard]] virtual auto number_of_nodes() const -> UNSIGNED_INTEGER_TYPE = 0;

        [[nodiscard]] virtual auto canonical_dimension() const -> UNSIGNED_INTEGER_TYPE= 0;

        [[nodiscard]] virtual auto number_of_nodes_per_elements() const -> UNSIGNED_INTEGER_TYPE= 0;

        [[nodiscard]] virtual auto number_of_elements() const -> UNSIGNED_INTEGER_TYPE= 0;


        /*!
         * Fetch the coordinates (x, y and z) of the given node.
         *
         * \note This method will always return a 3D vector, even if the Mesh instance is in 1D or 2D. In these
         *       cases, the z (and y in 1D) components of the position vector will be 0.
         * \note This method is very inefficient. For time sensitive applications, consider using directly the
         *       Mesh::node() method of the derived Mesh class.
         */
        [[nodiscard]] virtual auto node(const UNSIGNED_INTEGER_TYPE & node_id) const -> Eigen::Vector3d  = 0;

};

}
