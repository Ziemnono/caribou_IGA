#pragma once

#include <string>
#include <iostream>
#include <utility>
#include <functional>
#include <unordered_map>

#include <Caribou/config.h>
#include <Caribou/Topology/SplinePatch.h>

#include <Caribou/Topology/IO/CoreNURBS.h>
//#include <Caribou/Topology/IO/IGACellType.h>

namespace caribou::topology::io {

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex = UNSIGNED_INTEGER_TYPE>
class NURBSReader {
    static_assert(Dimension == 1 or Dimension == 2 or Dimension == 3, "The NURBSReader can only read 1D, 2D or 3D fields.");
public:
    using PatchType = SplinePatch<Dimension>;
    // Indices storage format.
    using ElementsIndices = Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    // Knot ranges storage format.
    using Double_Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    // Extraction Matrices storage format.
    using ElementsExtraction = std::vector<Double_Matrix>; //

    /** Build a new NURBSReader instance by reading a txt file. */
    static auto Read(const std::string & filepath) -> NURBSReader;

    /** Print information about the current txt file. */
    void print (std::ostream &out) const;

    /** Build the actual unstructured mesh from the txt file. */
    [[nodiscard]]
    auto patch() const -> PatchType;

    ~NURBSReader(){
        free(p_reader);
    }

private:
    NURBSReader(std::string filepath, coreNurbs * reader);

    const std::string p_filepath;
    coreNurbs * p_reader;
    std::array<UNSIGNED_INTEGER_TYPE, Dimension> p_axes;
};

//extern template class NURBSReader<1, unsigned long long int>;
//extern template class NURBSReader<1, unsigned long int>;
//extern template class NURBSReader<1, unsigned int>;

extern template class NURBSReader<2, unsigned long long int>;
extern template class NURBSReader<2, unsigned long int>;
extern template class NURBSReader<2, unsigned int>;

extern template class NURBSReader<3, unsigned long long int>;
extern template class NURBSReader<3, unsigned long int>;
extern template class NURBSReader<3, unsigned int>;

template<UNSIGNED_INTEGER_TYPE Dimension>
auto operator<<(std::ostream& os, const caribou::topology::io::NURBSReader<Dimension> & t) -> std::ostream&
{
    t.print(os);
    return os;
}

} /// namespace caribou::topology::io
