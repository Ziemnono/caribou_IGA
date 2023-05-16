#pragma once

#ifdef LEGACY_CXX
#include <experimental/filesystem>
namespace fs = ::std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = ::std::filesystem;
#endif

#include <bitset>
#include<vector>
#include<string>

#include <Caribou/Topology/IO/NURBSReader.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/BezierCrv.h>
#include <Caribou/Geometry/BezierSurf.h>


namespace caribou::topology::io {

template<UNSIGNED_INTEGER_TYPE Dimension>
auto nurbs_extract_axes_from_3D_vectors(const Double_Matrix & input_points, const int & number_of_points) -> std::array<UNSIGNED_INTEGER_TYPE, Dimension>;

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
NURBSReader<Dimension, NodeIndex>::NURBSReader(std::string filepath, coreNurbs<NodeIndex> * reader, PatchType * patch)
: p_filepath(std::move(filepath)), p_reader(std::move(reader)), p_patch(std::move(patch)){ }

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
auto NURBSReader<Dimension, NodeIndex>::Read(const std::string &filepath) -> NURBSReader<Dimension, NodeIndex> {
    // Make sure filepath is a valid path
    if (not fs::exists(filepath)) {
        throw std::invalid_argument("File '" + filepath + "' does not exists or cannot be read.");
    }
    // Get all data from the file
    coreNurbs<NodeIndex> * reader = new coreNurbs<NodeIndex>();
//    std::unique_ptr<coreNurbs> reader = std::make_unique<coreNurbs>();
    reader->SetFileName(filepath);
    reader->Update();

    // ------------------------  SplinePatch pointer -------------------------
    Double_Matrix nodes = reader->GetPoints();
    // Import Weights
    Double_Vector weights = reader->GetWeights();
    const auto number_of_elements = reader->GetNumberOfElements(); // Total Number of elements
    const auto number_of_nodes_per_element = reader->GetNumberOfElementPoints(); // Number of nodes per element
    // Indices of the elements
    Matrix<NodeIndex> indices;
    indices.resize(number_of_elements, number_of_nodes_per_element);
    indices = reader->GetIndices();

    // Knot range of elements
    Double_Matrix knot_ranges;
    knot_ranges.resize(number_of_elements, 4);
    knot_ranges = reader->GetKnotRanges();
    Double_Matrix a;
    DynVector knot_1(reader->get_knot_u_size());
    DynVector knot_2(reader->get_knot_v_size());
    knot_1 = reader->get_knot_u();
    knot_2 = reader->get_knot_v();

    PatchType * m =  new PatchType(nodes, weights, indices, knot_1, knot_2, knot_ranges);
//    const auto axes = nurbs_extract_axes_from_3D_vectors<Dimension>(reader->GetPoints(), reader->GetNumberOfPoints());
    NURBSReader<Dimension, NodeIndex> nr(filepath, reader, m);
//    free reader;
    return nr;
}

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
auto NURBSReader<Dimension, NodeIndex>::patch () const -> PatchType{
//    using WorldCoordinates = typename SplinePatch<Dimension>::WorldCoordinates;
//    using WeightContainer = typename SplinePatch<Dimension>::WeightsContainer;

    PatchType m;
    const auto number_of_nodes = p_reader->GetNumberOfPoints();
    if (number_of_nodes == 0) {
        return m;
    }

//    // Import nodes
//    std::vector<WorldCoordinates> nodes;
//    nodes.resize(p_reader->GetNumberOfPoints());

//    for (int i = 0; i < static_cast<int>(number_of_nodes); ++i) {
//        auto v = p_reader->GetPoint(i);
//        for (int axis = 0; axis < static_cast<int>(Dimension); ++axis) {
//            nodes[i][axis] = v(p_axes[axis]);
//        }
//    }

//    for (int i = 0; i < static_cast<int>(number_of_nodes); ++i) {
//        auto v = p_reader->GetPoint(i);
//        for (int axis = 0; axis < static_cast<int>(Dimension); ++axis) {
//            nodes[i][axis] = v(axis);
//        }
//    }

    Double_Matrix nodes = p_reader->GetPoints();
    // Import Weights
    Double_Matrix weights = p_reader->GetWeights();
    // Import elements
    const auto number_of_elements = p_reader->GetNumberOfElements(); // Total Number of elements

    const auto number_of_nodes_per_element = p_reader->GetNumberOfElementPoints(); // Number of nodes per element
    // Indices of the elements
    Matrix<NodeIndex> indices;
    indices.resize(number_of_elements, number_of_nodes_per_element);
    indices = p_reader->GetIndices();

    // Knot range of elements
    Double_Matrix knot_ranges;
    knot_ranges.resize(number_of_elements, 4);
    knot_ranges = p_reader->GetKnotRanges();
    Double_Matrix a;

    DynVector knot_1(p_reader->get_knot_u_size());
    DynVector knot_2(p_reader->get_knot_v_size());
    knot_1 = p_reader->get_knot_u();
    knot_2 = p_reader->get_knot_v();

    m = PatchType(nodes, weights, indices, knot_1, knot_2, knot_ranges);
    return m;
}



template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
auto NURBSReader<Dimension, NodeIndex>::patch_ptr () const -> const PatchType * {
    return p_patch.get();
}

/**
 * Extract the axes from an array of 3D coordinates.
 *
 * This method will first check which axis (x, y or z) has all of its coordinates equals.
 *
 * For example, the input
 *    [[24, 0, 12],
 *     [64, 0, 22],
 *     [51, 0, 9]]
 * will produce the following axis for a 2D mesh
 *    [0, 2] // x is the 0-axis, y is the 2-axis
 * and will raise an error for a 1D mesh
 *
 * The input
 *    [[5, 0, 12],
 *     [5, 0, 22],
 *     [5, 0, 9]]
 * will produce the following axis for a 1D mesh
 *    [2] // x is the 2-axis
 * and will raise an error for a 2D mesh
 *
 * @tparam Dimension The dimension of the output coordinates]
 * @param input_points Input array of 3D vectors
 * @param number_of_points Number of points in the input array
 */
template<UNSIGNED_INTEGER_TYPE Dimension>
auto nurbs_extract_axes_from_3D_vectors(const Double_Matrix & input_points, const int & number_of_points) -> std::array<UNSIGNED_INTEGER_TYPE, Dimension>
{
    // Find the good axes (only for dimensions 1 and 2)
    std::array<UNSIGNED_INTEGER_TYPE, Dimension> axes {};
    if (Dimension == 3) {
      axes[0] = 0;
      axes[1] = 1;
      axes[2] = 2;
    } else {
        auto w = input_points.row(0);
        std::bitset<3> is_all_the_same(0b111); // Set all to '1' at first
        std::array<double, 3> last_value {w(0), w(1), w(2)};
        for (int i = 1; i < number_of_points; ++i) {
            auto v = input_points.row(i);
            for (int axis = 0; axis < 3; ++axis) {
                if (last_value[axis] != v(axis))
                    is_all_the_same[axis] = false;
            }
        }

        if (is_all_the_same.flip().count() != Dimension) {
            throw std::runtime_error("Unable to convert a 3D field to a " + std::to_string(Dimension) +
                                     "D field. Their is " + std::to_string(is_all_the_same.count()) +
                                     " axes from the input mesh that have the same value.");
        }

        // Flip back
        is_all_the_same.flip();

        for (int axis = 0, c = 0; axis < 3; ++axis) {
            if (not is_all_the_same[axis]) {
                axes[c++] = axis;
            }
        }
    }

    return axes;
}

//template class NURBSReader<1, unsigned long long int>;
//template class NURBSReader<1, unsigned long int>;
//template class NURBSReader<1, unsigned int>;

template class NURBSReader<2, unsigned long long int>;
template class NURBSReader<2, unsigned long int>;
template class NURBSReader<2, unsigned int>;

template class NURBSReader<3, unsigned long long int>;
template class NURBSReader<3, unsigned long int>;
template class NURBSReader<3, unsigned int>;
}
