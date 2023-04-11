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
using Vector = Eigen::Matrix<dtype, 1, Eigen::Dynamic, Eigen::RowMajor>;

template<typename dtype>
using Matrix = Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Double_Vector = Vector<double>;
using Int_Vector = Vector<int>;
using Double_Matrix = Matrix<double>;
using Int_Matrix = Matrix<int>;
using USInt_Matrix = Matrix<UNSIGNED_INTEGER_TYPE>;


// NURBS topology storage
class para_topo{
    public:
//        using RangeMatrix = Matrix<double>;
//        using ConnectMatrix = Matrix<int>;
        para_topo(const Double_Matrix & elrange, const USInt_Matrix & connect){
            p_elrange.resize(elrange.rows(), elrange.cols());
            p_elrange = elrange;
            p_connect.resize(connect.rows(), connect.cols());
            p_connect = connect; // Element connectivity
        }

        Double_Matrix get_elrange() const{ return p_elrange;}

        USInt_Matrix get_elconn() const { return p_connect;}

        void display(void)
        {
            std::cout << "Elrange matrix\n" << p_elrange;
            std::cout << "Connectivity matrix\n" << p_connect;
        }
    private:
        Double_Matrix p_elrange; // Element wise knot range. Each row having element [U1, V1, U2, V2].
        USInt_Matrix p_connect; // Element wise connectivity.
};


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
        nelems_u = num_elements(knot_u); // u elements
        nelems_v = num_elements(knot_v); // v elements
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
        init_extraction();
    }

    void init_topo(void){
        // Connectivity Generation
        para_topo topo_u = gen_topo(knot_u, nelems_u, p);
        para_topo topo_v = gen_topo(knot_v, nelems_v, q);

        Double_Matrix elrange_u = topo_u.get_elrange();
        USInt_Matrix elconn_u = topo_u.get_elconn();
        Double_Matrix elrange_v = topo_v.get_elrange();
        USInt_Matrix elconn_v = topo_v.get_elconn();

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

    void init_extraction(void){
        C1 = extraction(p, nelems_u, knot_u);
        C2 = extraction(q, nelems_v, knot_v);
    }

    Double_Matrix GetExtraction(const int& elem_index){
        if (elem_index > t_nelems)
            throw std::invalid_argument("Element index shouldn't exceed total number of elements");

        int index_u = elem_index%(p+1);
        int index_v = elem_index/(p+1);
        return kron(C2[index_v], C1[index_u]);
    }

//    IGACellType GetCellType(void) const {return BezierSurf;};
    int GetExtractionSize(void) const {return (p+1)*(q+1);};
    int GetP(void) const {return p;};                                  // Degree p
    int GetQ(void) const {return q;};                                  // Degree q
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
    Double_Matrix GetPoints() const {return pnts;};                    // Control points
    Double_Vector GetPoint(const UNSIGNED_INTEGER_TYPE & i) const {return pnts.row(i);};     // Control points
    Double_Vector GetWeights(void) const {return wgts;};               // Weights
    FLOATING_POINT_TYPE GetWeight(const UNSIGNED_INTEGER_TYPE & i) const {return wgts(i);};
    Double_Matrix GetKnotRanges(void) const {return elRange;};          // Element parametric range
    USInt_Matrix GetIndices(void) const {return elConn;};                 // Element connectivity

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
    USInt_Matrix elConn;
    std::vector<Double_Matrix> C1, C2;
    auto no_nonzeros(const Double_Vector & vect) -> int;
    auto num_elements(const Double_Vector & knot_vector) -> int;
    para_topo gen_topo(const Double_Vector & knotVec, const int & no_elems, const int & p);
    std::vector<Double_Matrix> extraction(const int& utils_p, const int& num_elems, const Double_Vector & knot);
    Double_Matrix kron(const Double_Matrix & A, const Double_Matrix & B);

};

/*
INPUT  :: Accepts a vector
OUTPUT :: Number of non-zero values in the input vector.
*/
auto coreNurbs::no_nonzeros(const Double_Vector & vect) -> int {
    int count = 0;

    for (size_t i = 0; i< static_cast<size_t>(vect.cols()); i++)
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
auto coreNurbs::num_elements(const Double_Vector & knot_vector) -> int{
    int n_elems = 0; // Number of elements
    for (size_t i = 1; i < static_cast<size_t>(knot_vector.cols()); i++)
    {
        if (knot_vector[i-1] != knot_vector[i]){
            n_elems = n_elems+1;
        }
    }
    return n_elems;
}

para_topo coreNurbs::gen_topo(const Double_Vector & knotVec, const int & no_elems, const int & util_p){
    Double_Matrix util_elRange(no_elems, 2);
    Int_Matrix util_elKnotIndices(no_elems, 2);
    USInt_Matrix util_elConn(no_elems, p+1);

    int elem;
    elem = 0;
    float currentKnotVal, previousKnotVal;
    previousKnotVal = 0.0;

    for (int i = 0; i < knotVec.cols(); i++)
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

    para_topo topo(util_elRange, util_elConn);
    return topo;
}


// Extraction operation
std::vector<Double_Matrix> coreNurbs::extraction(const int& utils_p, const int& num_elems, const Double_Vector & knot){
    int knot_len = knot.cols();
    float nume, alpha;
    int i, j, r, s, k, save, mult;
    int m = knot_len-utils_p-1;
    // int num_elems = std::set<Scalar>(knot.array().data(), knot.array().data()+knot.array().size()).size() - 1;
    int a = utils_p+1;
    int b = a+1;
    int nb = 0;
    std::vector<Double_Matrix> C(num_elems);
    Double_Matrix ones = Double_Matrix::Identity(utils_p+1, utils_p+1);
    C[0] = ones;
    Double_Vector alphas(utils_p+1);
    while (b <= m) {
        C[nb+1] = ones;
        i = b;
        while ((b <= m) && (knot[b] == knot[b-1])) { b = b+1; }
        mult = b - i + 1;
        if (mult < utils_p){
            nume = knot[b-1] - knot[a-1];
            for (j = utils_p; j > mult; j--) {
                alphas[j-mult] = nume/(knot[a+j-1]-knot[a-1]);
            }
            r = utils_p-mult;
            for (j = 1; j < r+1; ++j) {
                save = r-j+1;
                s = mult + j;
                for (k = utils_p+1; k > s; k--) {
                    alpha = alphas[k-s];
                    C[nb].col(k-1) = alpha * C[nb].col(k-1) + (1 - alpha) * C[nb].col(k-2);
                }
                if (b <= m) {
                    // C[nb+1](Eigen::seq(save-1, save+j-1), save-1) = C[nb](Eigen::seq(p-j, p), p);
                    // instead of using "Eiigen::seq" function below loop is used.
                    for (int itr = 0; itr <= j; itr++){
                        C[nb+1](save-1+itr, save-1) = C[nb](utils_p-j+itr, utils_p);
                    }
                }
            }
            nb = nb+1;
            if (b<= m){
                a = b;
                b = b+1;
            }
        }
        else if (mult == utils_p){
            if (b <= m){
                nb = nb+1;
                a = b;
                b = b + 1;
            }
        }
    }
    return C;
}

// Matrix Tensor Product
Double_Matrix coreNurbs::kron(const Double_Matrix & A, const Double_Matrix & B) {
    /*
    a = []; b = []

    a11 * b   a12 * b . . .

    a21 * b   a22 * b . . .
    .            .    .
    .            .      .
    .            .        .

    */
    const int Ar = A.rows();
    const int Ac = A.cols();
    const int Br = B.rows();
    const int Bc = B.cols();

    Double_Matrix M(Ar*Br, Ar*Br);

    for (int k = 0; k < Ar; k++){
        for (int l = 0; l < Ac; l++){
            for (int i = 0; i < Br; i++){
                for (int j = 0; j < Bc; j++){
                    // std::cout << i + k*Br << "-" << j + l * Bc << "\n";
                    M(i + k*Br, j + l * Bc) = B(i,j) * A(k,l);
    }}}}
    return M;
}



} /// namespace caribou::topology::io::nurbs
