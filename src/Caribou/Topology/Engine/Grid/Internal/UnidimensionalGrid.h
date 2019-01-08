#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_UNIDIMENSIONALGRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_UNIDIMENSIONALGRID_H

#include <Caribou/Topology/Engine/Grid/Internal/Grid.h>

namespace caribou {
namespace topology {
namespace engine {
namespace internal {

/**
 * Simple representation of an unidimensional (1D) Grid in space.
 *
 * ** Do not use this class directly. Use instead caribou::topology::engine::Grid. **
 *
 * The functions declared in this class can be used with any type of grids (static grid, container grid, etc.).
 *
 * Do to so, it uses the Curiously recurring template pattern (CRTP) :
 *    https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 *
 * A Grid is a set of multiple cell entities (same-length lines in 1D, rectangles
 * in 2D and rectangular hexahedrons in 3D) aligned in the x, y and z axis.
 *
 * @tparam GridType_ Type of the derived grid class that will implement the final functions.
 */
template <class GridType_>
struct BaseUnidimensionalGrid : public BaseGrid<1, GridType_>
{
    static constexpr size_t Dimension = 1;

    using GridType = GridType_;

    using Int = size_t;
    using Float = FLOATING_POINT_TYPE;

    using Index = Int;
    using NodeIndex = Int;
    using CellIndex = Int;
    using RelativePosition = Float;
    using Size = Float;

    using BaseGrid<Dimension, GridType_>:: BaseGrid;


private:
    inline constexpr
    const GridType &
    Self() const
    {
        return static_cast<const GridType &> (*this);
    }

};

} // namespace internal
} // namespace engine
} // namespace topology
} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_INTERNAL_UNIDIMENSIONALGRID_H
