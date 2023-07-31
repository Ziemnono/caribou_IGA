#pragma once

#include <Caribou/config.h>

#include<fstream>
#include<vector>
#include <string>
#include <iostream>
#include <utility>
#include <functional>
#include <unordered_map>

#include <Eigen/Core>

//#include <Caribou/Topology/IO/IGACellType.h>

namespace caribou::topology::io {

template<typename dtype>
using Vector = Eigen::Matrix<dtype, Eigen::Dynamic, 1>;

template<typename dtype>
using Matrix = Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic>;

using Double_Vector = Vector<FLOATING_POINT_TYPE>;
using Int_Vector = Vector<int>;
using Double_Matrix = Matrix<FLOATING_POINT_TYPE>;
using Int_Matrix = Matrix<int>;
using USInt_Matrix = Matrix<UNSIGNED_INTEGER_TYPE>;
using USIVector = Eigen::Matrix<UNSIGNED_INTEGER_TYPE, Eigen::Dynamic, 1>;

namespace utils{

// NURBS topology storage
template <typename NodeIndex = UNSIGNED_INTEGER_TYPE>
class para_topo{
    public:
//        using RangeMatrix = Matrix<double>;
//        using ConnectMatrix = Matrix<int>;
        para_topo(const Double_Matrix & elrange, const Matrix<NodeIndex> & connect){
            p_elrange.resize(elrange.rows(), elrange.cols());
            p_elrange = elrange;
            p_connect.resize(connect.rows(), connect.cols());
            p_connect = connect; // Element connectivity
        }

        Double_Matrix get_elrange() const{ return p_elrange;}

        Matrix<NodeIndex> get_elconn() const { return p_connect;}

        void display(void)
        {
            std::cout << "Elrange matrix\n" << p_elrange;
            std::cout << "Connectivity matrix\n" << p_connect;
        }
    private:
        Double_Matrix p_elrange; // Element wise knot range. Each row having element [U1, V1, U2, V2].
        Matrix<NodeIndex> p_connect; // Element wise connectivity.
};

/*
INPUT  :: Accepts a vector
OUTPUT :: Number of non-zero values in the input vector.
*/
auto no_nonzeros(const Double_Vector & vect) -> int {
    int count = 0;

    for (size_t i = 0; i< static_cast<size_t>(vect.size()); i++)
    {
        if (vect[i] == 0.0){continue;}
        count = count + 1;
    }
    return count;
}

/*
INPUT  :: Accepts a vector (knot vector).
OUTPUT :: Number of possible elements from the vector.
*/
auto num_elements(const Double_Vector & knot_vector) -> int{
    int n_elems = 0; // Number of elements
    for (size_t i = 1; i < static_cast<size_t>(knot_vector.size()); i++)
    {
        if (knot_vector[i-1] != knot_vector[i]){
            n_elems = n_elems+1;
        }
    }
    return n_elems;
}
template <typename NodeIndex = UNSIGNED_INTEGER_TYPE>
para_topo<NodeIndex> gen_topo(const Double_Vector & knotVec, const int & no_elems, const int & util_p){
    Double_Matrix util_elRange(no_elems, 2);
    Int_Matrix util_elKnotIndices(no_elems, 2);
    Matrix<NodeIndex> util_elConn(no_elems, util_p+1);

    int elem;
    elem = 0;
    float currentKnotVal, previousKnotVal;
    previousKnotVal = 0.0;

    for (int i = 0; i < knotVec.size(); i++)
    {
        currentKnotVal = knotVec(i);
        if (knotVec(i) != previousKnotVal){
            util_elRange.row(elem) << previousKnotVal, currentKnotVal;
            util_elKnotIndices.row(elem) << i-1, i;
            elem = elem + 1;
        }
        previousKnotVal = currentKnotVal;
    }

    int numRepeatedKnots  = 0;

    Double_Vector previousKnotVals(util_p), currentKnotVals(util_p);
    Int_Vector indices(util_p);
    Double_Vector ones = Double_Vector::Constant(util_p, 1.0);

    for (int e = 0; e < no_elems; e++)
    {
        for (size_t j = 0; j < static_cast<size_t>(util_p); j++)
        {
            indices(j) = util_elKnotIndices(e,0)-util_p+1+j;
        }

        for (size_t j = 0; j < static_cast<size_t>(util_p); j++)
        {
            previousKnotVals(j) = knotVec(indices(j));
            currentKnotVals(j) = knotVec(util_elKnotIndices(e,0));
        }

        if ((previousKnotVals == currentKnotVals) && (1 < no_nonzeros(previousKnotVals))){
            numRepeatedKnots = numRepeatedKnots + 1;
        }
        for (int j = 0; j <= util_p; j++)
        {
            util_elConn(e,j) = util_elKnotIndices(e,0)-util_p+j;
        }

    }

    para_topo<NodeIndex> topo(util_elRange, util_elConn);
    return topo;
}

}  // Utils namespace

template <typename NodeIndex = UNSIGNED_INTEGER_TYPE>
struct coreNurbs{

public:
    void SetFileName(const std::string & filename){
        p_filename = filename;
    };

    void SetFileName(const char * filename){
        p_filename = filename;
    };

    void Update(void){
        std::ifstream file;
        file.open(p_filename, std::ios::in);
        std::string s;

        for (int i = 0; i < 4; ++i) {
            getline(file, s);
        }
        file >> pdim;
        file >> rdim;

        getline(file, s);
        getline(file, s);

        // NURBS degrees p and q.
        file >> p;
        file >> q;
        // No of control points in each parametric directions.
        file >> cp_u;
        file >> cp_v;
        // Skipping a line.
        getline(file, s);
        // Knot vecto length of each parametric dir.
        len_knot_u = cp_u+p+1;
        len_knot_v = cp_v+q+1;
        // Knot vectors.
        knot_u.resize(len_knot_u);
        knot_v.resize(len_knot_v);
        // Writing knot_U data.
        for (int i = 0; i < len_knot_u; i++){ file >> knot_u(i); }
        // Writing knot_v data.
        for (int i = 0; i < len_knot_v; i++){ file >> knot_v(i); }

        // Control point information
        t_cp = cp_u*cp_v; // Total control points.
        pnts.resize(t_cp, rdim); // Control points array.
        wgts.resize(t_cp); // Weights array.

        // Points
        for (int j = 0; j < rdim; j++)
        {
            for (int i = 0; i < t_cp; i++)
            {
                file >> pnts(i,j);
            }
        }
        // Weights
        for (int i = 0; i < t_cp; i++){ file >> wgts(i); }

        file.close();

        // elements
        nelems_u = utils::num_elements(knot_u); // u elements
        nelems_v = utils::num_elements(knot_v); // v elements
        t_nelems = nelems_u * nelems_v; // Total number of elements.

        // Set Global Indices
        global_indices.resize(cp_v, cp_u);
        for (int j = 0; j < cp_v; j++){
            for (int i = 0; i < cp_u; i++){
                global_indices(j,i) = i + j*cp_u;
            }
        }


        // generate Topology
        init_topo();
    }

    void init_topo(void){
        // Connectivity Generation
        utils::para_topo<NodeIndex> topo_u = utils::gen_topo<NodeIndex>(knot_u, nelems_u, p);
        utils::para_topo<NodeIndex> topo_v = utils::gen_topo<NodeIndex>(knot_v, nelems_v, q);

        Double_Matrix elrange_u = topo_u.get_elrange();
        Matrix<NodeIndex> elconn_u = topo_u.get_elconn();
        Double_Matrix elrange_v = topo_v.get_elrange();
        Matrix<NodeIndex> elconn_v = topo_v.get_elconn();

        elRange.resize(t_nelems, 4);
        elConn.resize(t_nelems, (p+1)*(q+1));
        int e = 0, c;
        for (int v = 0; v < nelems_v; v++){
            for (int u = 0; u < nelems_u; u++){
                elRange.row(e) << elrange_u(u,0), elrange_v(v,0), elrange_u(u,1), elrange_v(v,1);
                c = 0;
                for (int i = 0; i < q+1; i++){
                    for (int j = 0; j < p+1; j++){
                        elConn(e,c) = global_indices(elconn_v(v,i), elconn_u(u,j));
                        c = c + 1;
                    }
                }
                e = e+1;
            }
        }
    }


//    IGACellType GetCellType(void) const {return BezierSurf;};
    int GetNodesPerElement(void) const {return (p+1)*(q+1);};
    UNSIGNED_INTEGER_TYPE GetP(void) const {return p;};                                  // Degree p
    UNSIGNED_INTEGER_TYPE GetQ(void) const {return q;};                                  // Degree q
    int GetReadDim(void) {return rdim;};                               // Dimensions of control points
    int GetNumberOfElementPoints(void) const {return (p+1)*(q+1);};    // Number of nodes per cell
    int get_no_pnts_u(void) const {return cp_u;};                      // Control points u
    int get_no_pnts_v(void) const {return cp_v;};                      // Control points v
    int GetNumberOfPoints(void) const {return t_cp;};
    int get_no_elems_u(void) const {return nelems_u;};                 // Number of elements u
    int get_no_elems_v(void) const {return nelems_v;};                 // Number of elements v
    int GetNumberOfElements(void) const {return t_nelems;};            // Total no of elements
    Double_Vector get_knot_u(void) const {return knot_u;};             // Knot vector u
    Double_Vector get_knot_v(void) const {return knot_v;};             // Knot vector v
    int get_knot_u_size(void) const {return knot_u.size();};             // Knot vector u
    int get_knot_v_size(void) const {return knot_v.size();};             // Knot vector v
    Double_Matrix GetPoints() const {return pnts;};                    // Control points
    Double_Vector GetPoint(const UNSIGNED_INTEGER_TYPE & i) const {return pnts.row(i);};     // Control points
    Double_Vector GetWeights(void) const {return wgts;};               // Weights
    FLOATING_POINT_TYPE GetWeight(const UNSIGNED_INTEGER_TYPE & i) const {return wgts(i);};
    Double_Matrix GetKnotRanges(void) const {return elRange;};          // Element parametric range
    Matrix<NodeIndex> GetIndices(void) const {return elConn;};                 // Element connectivity
    void print_core_nurbs(void){
        std::cout << "P and q : " << p << " " << q << "\n";
        std::cout << "Knot U  : " << knot_u.transpose() << "\n";
        std::cout << "Knot V  : " << knot_v.transpose() << "\n";
        std::cout << "Points in U & V  : " << cp_u << " " << cp_v << "\n";
        std::cout << "Elems in U & V  : " << nelems_u << " " << nelems_v << "\n";
        std::cout << "Total Elems : " << t_nelems << "\n";
        std::cout << "Points \n" << pnts.transpose() << "\n";
        std::cout << "Weights \n" << wgts.transpose() << "\n";
    }
private:
    std::string p_filename;
    int pdim, rdim;
    int p, q, cp_u, cp_v, len_knot_u, len_knot_v, t_cp, nelems_u, nelems_v, t_nelems;
    Double_Vector knot_u;
    Double_Vector knot_v;
    Double_Matrix pnts;
    Double_Vector wgts;
    Int_Matrix global_indices;
    Double_Matrix elRange;
    Matrix<NodeIndex> elConn;
};


} /// namespace caribou::topology::io::nurbs
