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


#include <Caribou/Topology/IO/CoreNURBS.h>

#include <Caribou/Topology/IO/IGACellType.h>

namespace caribou::topology::io {

template<UNSIGNED_INTEGER_TYPE Dimension>
auto nurbs_extract_axes_from_3D_vectors(const Double_Matrix & input_points, const int & number_of_points) -> std::array<UNSIGNED_INTEGER_TYPE, Dimension>;

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
NURBSReader<Dimension, NodeIndex>::NURBSReader(std::string filepath, coreNurbs * reader, std::array<UNSIGNED_INTEGER_TYPE, Dimension> axes)
: p_filepath(std::move(filepath)), p_reader(std::move(reader)), p_axes(axes)
{
    // Segments
    register_element_type<geometry::BezierCrv<Dimension>>(BezierCrv);

    if constexpr (Dimension > 1) {
        // Quads
        register_element_type<geometry::BezierSurf<Dimension>>(BezierSurf);
    }
}


template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
auto NURBSReader<Dimension, NodeIndex>::Read(const std::string &filepath) -> NURBSReader<Dimension, NodeIndex> {
    // Make sure filepath is a valid path
    if (not fs::exists(filepath)) {
        throw std::invalid_argument("File '" + filepath + "' does not exists or cannot be read.");
    }

    // Get all data from the file
    coreNurbs * reader;
    reader->SetFileName(filepath);
    reader->Update();

    const auto axes = nurbs_extract_axes_from_3D_vectors<Dimension>(reader->GetPoints(), reader->GetNumberOfPoints());


    return NURBSReader<Dimension, NodeIndex>(filepath, reader, axes);
}

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
auto NURBSReader<Dimension, NodeIndex>::mesh () const -> SplineMesh<Dimension> {
    using WorldCoordinates = typename SplineMesh<Dimension>::WorldCoordinates;
    using WeightContainer = typename SplineMesh<Dimension>::WeightsContainer;

    SplineMesh<Dimension> m;
    const auto number_of_nodes = p_reader->GetNumberOfPoints();
    if (number_of_nodes == 0) {
        return m;
    }

    // Import nodes
    std::vector<WorldCoordinates> nodes;
    nodes.resize(p_reader->GetNumberOfPoints());


    for (size_t i = 0; i < static_cast<size_t>(number_of_nodes); ++i) {
        auto v = p_reader->GetPoint(i);
        for (std::size_t axis = 0; axis < Dimension; ++axis) {
            nodes[i][axis] = v(p_axes[axis]);
        }
    }

    // Import Weights
    WeightContainer weights = p_reader->GetWeights();

    m = SplineMesh<Dimension> (nodes, weights);

    // Import elements
    const auto type = p_reader->GetCellType();
    const auto number_of_elements = p_reader->GetNumberOfElements(); // Total Number of elements

    const auto number_of_nodes_per_element = p_reader->GetNumberOfElementPoints(); // Number of nodes per element
    // Indices of the elements
    Eigen::Matrix<UNSIGNED_INTEGER_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> indices;
    indices.resize(number_of_elements, number_of_nodes_per_element);
    indices = p_reader->GetIndices();

    // Knot range of elements
    Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> knot_ranges;
    knot_ranges.resize(number_of_elements, 4);
    knot_ranges = p_reader->GetKnotRanges();

    // Extraction matrix
    using ElementExtraction = std::vector<Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;
    ElementExtraction element_extractions(number_of_elements);
//    const auto extraction_size = p_reader->GetExtractionSize();
    for (int i = 0; i < number_of_elements; ++i) {
//        element_extractions[i].resize(extraction_size, extraction_size);
        element_extractions[i] = p_reader->GetExtraction(i);
    }


//    const auto & domain_builder = p_domain_builders.at(static_cast<IGACellType>(type));
//    domain_builder(m, indices, knot_ranges, element_extractions);
    std::string name = "domain_bzr";
    m.template add_domain<geometry::BezierSurf<Dimension>, NodeIndex>(name, indices, knot_ranges, element_extractions);

    return m;
}


//template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
//void NURBSReader<Dimension, NodeIndex>::print (std::ostream & out) const {
//    vtkUnstructuredGrid * output = p_reader->GetOutput();

//    out << "input has " << output->GetNumberOfPoints() << " points.\n";
//    out << "input has " << output->GetNumberOfCells() << " cells.\n";

//    vtkSmartPointer <vtkCellTypes> types = vtkSmartPointer <vtkCellTypes>::New();
//    output->GetCellTypes(types);
//    out << types->GetNumberOfTypes() << " types\n";
//    for (unsigned int i = 0; i < types->GetNumberOfTypes(); ++i) {
//        const auto type = types->GetCellType(i);
//        auto cells = vtkSmartPointer <vtkIdTypeArray>::New();
//        output->GetIdsOfCellsOfType(type, cells);
//        out << vtkCellTypes::GetClassNameFromTypeId(type) << " : " << cells->GetDataSize() << "\n";
//    }

//    vtkIdType numberOfCellArrays = output->GetCellData()->GetNumberOfArrays();
//    if (numberOfCellArrays > 0) {
//        out << "Cell data arrays:\n";
//        for (vtkIdType i = 0; i < numberOfCellArrays; i++) {
//            const auto dataTypeID = output->GetCellData()->GetArray(i)->GetDataType();
//            const auto dataTypeIDStr = output->GetCellData()->GetArray(i)->GetDataTypeAsString();
//            out << "  Array " << i << ": "
//                << output->GetCellData()->GetArrayName(i)
//                << "  (type: " << dataTypeIDStr << " - " << dataTypeID << ")\n";
//        }
//    }

//    vtkIdType numberOfFieldArrays = output->GetFieldData()->GetNumberOfArrays();
//    if (numberOfFieldArrays) {
//        out << "Field data arrays:\n";
//        for (vtkIdType i = 0; i < numberOfFieldArrays; i++) {
//            const auto dataTypeID = output->GetFieldData()->GetArray(i)->GetDataType();
//            const auto dataTypeIDStr = output->GetFieldData()->GetArray(i)->GetDataTypeAsString();
//            out << "  Array " << i << ": "
//                << output->GetFieldData()->GetArrayName(i)
//                << "  (type: " << dataTypeIDStr << " - " << dataTypeID << ")\n";
//        }
//    }
//}

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
        auto v = input_points.row(0);
        std::bitset<3> is_all_the_same(0b111); // Set all to '1' at first
        std::array<double, 3> last_value {v[0], v[1], v[2]};
        for (size_t i = 1; i < number_of_points; ++i) {
            v = input_points.row(i);
            for (std::size_t axis = 0; axis < 3; ++axis) {
                if (last_value[axis] != v[axis])
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

        for (std::size_t axis = 0, c = 0; axis < 3; ++axis) {
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

//template class NURBSReader<2, unsigned long long int>;
//template class NURBSReader<2, unsigned long int>;
template class NURBSReader<2, unsigned int>;

//template class NURBSReader<3, unsigned long long int>;
//template class NURBSReader<3, unsigned long int>;
template class NURBSReader<3, unsigned int>;
}
