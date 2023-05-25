#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseNurbsSurf.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct NurbsSurf;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<NurbsSurf <_Dimension>> {
     // static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; // NurbsSurf degree
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 2; // Two parametric direction.
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 3*3; // NurbsSurf control points (p+1)*(q+1)
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 3*3; // Gauss points per element
 };

 /**
  * Bezier Surface two canonicalDimensions
  *
  * \verbatim
  *        v
  *        ^
  *        |
  *  6-----7-----8
  *  |     |     |
  *  |     |     |
  *  3     4---- 5 --> u
  *  |           |
  *  |           |
  *  0-----1-----2
  *  \endverbatim
  *
  * @tparam _Dimension The world coordinates dimension
  */

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct NurbsSurf: public BaseNurbsSurf<NurbsSurf <_Dimension>> {
    // Types

    using Scalar = FLOATING_POINT_TYPE;

    using Base = BaseNurbsSurf<NurbsSurf <_Dimension>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    using Dyn_Matrix = typename Base::Dyn_Matrix;
    using Dyn_Vector = typename Base::Dyn_Vector;

    static constexpr UNSIGNED_INTEGER_TYPE Degree_1 = 2; //NurbsSurf
    static constexpr UNSIGNED_INTEGER_TYPE Degree_2 = 2; //NurbsSurf
    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime; //$$$NurbsSurf
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime; //$$$NurbsSurf

    static constexpr FLOATING_POINT_TYPE canonical_nodes [NumberOfNodesAtCompileTime][CanonicalDimension] = {
    //    u,  v
        {-1, -1}, // Node 0
        {+0, -1}, // Node 1
        {+1, -1}, // Node 2
        {-1, +0}, // Node 3
        {+0, +0}, // Node 4
        {+1, +0}, // Node 5
        {-1, +1}, // Node 6
        {+0, +1}, // Node 7
        {+1, +1}, // Node 8

    };

    // Constructors
    using Base::Base;
    NurbsSurf() : Base() {
        if constexpr (Dimension == 2) {
            // WorldCoordinates(x, y, w)
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1]); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1]); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1]); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1]); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(canonical_nodes[4][0], canonical_nodes[4][1]); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(canonical_nodes[5][0], canonical_nodes[5][1]); // Node 50.5 *
            this->p_nodes.row(6) = WorldCoordinates(canonical_nodes[6][0], canonical_nodes[6][1]); // Node 6
            this->p_nodes.row(7) = WorldCoordinates(canonical_nodes[7][0], canonical_nodes[7][1]); // Node 7
            this->p_nodes.row(8) = WorldCoordinates(canonical_nodes[8][0], canonical_nodes[8][1]); // Node 8
        }
        else {
            // WorldCoordinates(x, y, z, w)
            this->p_nodes.row(0) = WorldCoordinates(canonical_nodes[0][0], canonical_nodes[0][1], 0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates(canonical_nodes[1][0], canonical_nodes[1][1], 0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates(canonical_nodes[2][0], canonical_nodes[2][1], 0); // Node 2
            this->p_nodes.row(3) = WorldCoordinates(canonical_nodes[3][0], canonical_nodes[3][1], 0); // Node 3
            this->p_nodes.row(4) = WorldCoordinates(canonical_nodes[4][0], canonical_nodes[4][1], 0); // Node 4
            this->p_nodes.row(5) = WorldCoordinates(canonical_nodes[5][0], canonical_nodes[5][1], 0); // Node 5
            this->p_nodes.row(6) = WorldCoordinates(canonical_nodes[6][0], canonical_nodes[6][1], 0); // Node 6
            this->p_nodes.row(7) = WorldCoordinates(canonical_nodes[7][0], canonical_nodes[7][1], 0); // Node 7
            this->p_nodes.row(8) = WorldCoordinates(canonical_nodes[8][0], canonical_nodes[8][1], 0); // Node 8
        };
    }

    // Function to calculate knot span of a parametric point
    int findspan(const float& u, const int& p, const Dyn_Vector& U) const {
        // n -> Number of control points.
        // u    -> Parametric points.
        // U -> Knot vector.

        int n = U.size() - p - 1; // No of control points
        if (u == U[n+1])
        {
            return n;
        }
        int low, high, mid;
        low = p;
        high = n+1;
        mid = floor((low + high)/2);
        while (u < U[mid] || u >= U[mid+1]) {
            if (u < U[mid]) {
                high = mid;
            }
            else {
                low = mid;
            }
            mid = floor((low+high)/2);
        }

        return mid;
    }

    // Function to calculate element-wise basis functions at given parametric point.
    Dyn_Vector bspbasisfun(const Scalar& u,const int& p,const Dyn_Vector& U) const {
        Dyn_Vector N;   N.resize(p+1);  N.setZero();
        if (u == 1.0){
            N[p] = 1.0;
            return N;
        }

        int span = findspan(u, p, U);

        Dyn_Vector left;    left.resize(p+1);   left.setZero();
        Dyn_Vector right;   right.resize(p+1);  right.setZero();

        N[0] = 1.0;
        float saved, temp;
        for (int j = 1; j <= p; j++){
            left[j] = u - U[span+1-j];
            right[j] = U[span+j] - u;
            saved = 0.0;
            for (int r = 0; r<j; r++){
                temp = N[r]/(right[r+1]+left[j-r]);
                N[r] = saved + right[r+1]*temp;
                saved = left[j-r]*temp;
            }
            N[j] = saved;
        }
        return N;
    }



    // Function to calculate element-wise shape function derivative.
    Dyn_Vector basis_deri(const Scalar& u, const int & p, Dyn_Vector U) const {
        int Numbasis = p+1;
        // Copying the knotvector to U for easy usage.
//            int npts = U.size() - p - 1;
        static constexpr int n = 1; // Derivative order.
        int i = findspan(u, p, U);
        Dyn_Matrix ndu; ndu.resize(Numbasis, Numbasis); ndu.setZero();
        Dyn_Matrix ders; ders.resize(n+1, Numbasis); ders.setZero();
        Dyn_Matrix a;   a.resize(2, Numbasis); a.setZero();

        ndu(0,0) = 1.0;

        Dyn_Vector left; left.resize(Numbasis);
        Dyn_Vector right; right.resize(Numbasis);
        Scalar saved, temp;

        int j;
        for (j = 1; j <= p; j++)
        {
            left[j] = u - U[i+1-j];
            right[j] = U[i+j] - u;
            saved = 0.0;
            for (int r = 0; r < j; r++)
            {
                ndu(j,r) = right[r+1] + left[j-r];
                temp = ndu(r,j-1)/ndu(j,r);
                ndu(r,j) = saved + right[r+1]*temp;
                saved = left[j-r] * temp;
            }
            ndu(j,j) = saved;
        }

        for (j = 0; j <= p; j++)
        {
            ders(0,j) = ndu(j,p);
        }

        int s1, s2, rk, pk, j1, j2;
        Scalar d;
        for (int r = 0; r <= p; r++) {
            s1 = 0; s2 = 1;
            a(0,0) = 1.0;
            for (int k = 1; k <=n ; k++) {
                d = 0.0;
                rk = r-k;
                pk = p - k;
                if (r >= k) {
                    a(s2, 0) = a(s1,0)/ndu(pk+1, rk);
                    d = a(s2, 0)*ndu(rk, pk);
                }
                if (rk >= -1) {
                    j1 = 1;
                }
                else {
                    j1 = -rk;
                }

                if ((r-1) <= pk) {
                    j2 = k-1;
                }
                else {
                    j2 = p-r;
                }

                for (j = j1; j <= j2; j++) {
                    a(s2, j) = (a(s1, j)-a(s1, j-1))/ndu(pk+1, rk+j);
                    d = d + a(s2, j) * ndu(rk+j, pk);
                }

                if (r <= pk) {
                    a(s2, k) = -1*a(s1, k-1)/ndu(pk+1, r);
                    d = d + a(s2, k) * ndu(r, pk);
                }
                ders(k, r) = d;
                j = s1;
                s1 = s2;
                s2 = j;
            }
        }
        int r = p;
        for (int k = 1; k <= n; k++) {
            for (j = 0; j <= p; j++) {
                ders(k, j) = ders(k, j) * r;
            }
            r = r * (p-k);
        }
        return ders.row(1);
    }

    void print(void){
        std::cout << "Control Points : \n" << this->p_nodes.transpose() << "\n";
        std::cout << "Weight Points  : \n" << this->p_weights.transpose() << "\n";
        std::cout << "Knot Vector 1  : \n" << this->p_knot_1.transpose() << "\n";
        std::cout << "Knot Vector 2  : \n" << this->p_knot_2.transpose() << "\n";
        std::cout << "Knot span      : \n" << this->p_knot_span.transpose() << "\n";

    }

    inline auto jacobian_papa() -> Scalar {
        auto Jxi = 0.5 * (this->p_knot_span[2] - this->p_knot_span[0]);
        auto Jeta = 0.5 * (this->p_knot_span[3] - this->p_knot_span[1]);
        return Jxi * Jeta;
    }


private:
    // Implementations
    friend struct Element<NurbsSurf <_Dimension>>;
    friend struct BaseNurbsSurf<NurbsSurf <_Dimension>>;


    auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        // U direction basis
        Scalar u = 0.5 * (xi[0]*(this->p_knot_span[2] - this->p_knot_span[0]) + (this->p_knot_span[2] + this->p_knot_span[0]));
        Scalar v = 0.5 * (xi[1]*(this->p_knot_span[3] - this->p_knot_span[1]) + (this->p_knot_span[3] + this->p_knot_span[1]));

        auto knot1 = this->p_knot_1;
        auto basis_u = bspbasisfun(u, Degree_1, knot1);

        auto knot2 = this->p_knot_2;
        auto basis_v = bspbasisfun(v, Degree_2, knot2);

        Dyn_Matrix m = basis_v * basis_u.transpose();
        Eigen::Map<Vector<NumberOfNodesAtCompileTime>> basis(m.data(), m.size());
        basis = basis.array() * this->p_weights.array();
        basis = basis * (1/basis.sum());
        return basis;

    };

    inline auto get_dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {
        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> dL;
        Scalar u = 0.5 * (xi[0]*(this->p_knot_span[2] - this->p_knot_span[0]) + (this->p_knot_span[2] + this->p_knot_span[0]));
        Scalar v = 0.5 * (xi[1]*(this->p_knot_span[3] - this->p_knot_span[1]) + (this->p_knot_span[3] + this->p_knot_span[1]));

//        std::cout << "\nu - " <<  u << "\n";
//        std::cout << "v - " <<  v << "\n";
        auto knot1 = this->p_knot_1;
        auto knot2 = this->p_knot_2;

        Dyn_Vector Nbasis_u(Degree_1+1);
        Dyn_Vector Nbasis_v(Degree_2+1);

        Dyn_Vector Nderis_u(Degree_1+1);
        Dyn_Vector Nderis_v(Degree_2+1);

        Dyn_Matrix Nbasis(2,3); // Parametric direction wise basis functions   U & V
        Dyn_Matrix Nderis(2,3); // Parametric direction wise basis derivatives U & V
        // U direction basis & derivative
        Nbasis_u = bspbasisfun(u, Degree_1, knot1);
        Nderis_u = basis_deri(u, Degree_1, knot1);

        // V direction basis & derivative
        Nbasis_v = bspbasisfun(v, Degree_2, knot2);
        Nderis_v = basis_deri(v, Degree_2, knot2);

        Dyn_Matrix global_storage(3,3);
        global_storage =  Nbasis_v * Nbasis_u.transpose();
//        std::cout << "Golbal basis Matrix \n" <<  global_storage << "\n";

        Eigen::Map<Dyn_Vector> basis(global_storage.data(), global_storage.size());
//        std::cout << "Golbal basis \n" <<  basis.transpose() << "\n";

        Dyn_Matrix deris(9,2);
        // Bv_dBu
        Dyn_Matrix deri_m1(3,3);
        deri_m1 = Nbasis_v * Nderis_u.transpose();
        deris.col(0) = Eigen::Map<Dyn_Vector>(deri_m1.data(), deri_m1.size());

        // dBv_Bu
        Dyn_Matrix deri_m2(3,3);
        deri_m2 = Nderis_v * Nbasis_u.transpose();
        deris.col(1) = Eigen::Map<Dyn_Vector>(deri_m2.data(), deri_m2.size());

//        std::cout << "First deris \n" <<  deris << "\n";

        // Weights * Basis functions component-wise multiplication.
        Dyn_Vector B1 = this->p_weights.array() * basis.array();
        double w1 = B1.sum();

        // Calculating dB_du and dB_dv
        Dyn_Vector dB_du_B2 = this->p_weights.array() * deris.col(0).array();
        double dB_du_w2 = dB_du_B2.sum();
        deris.col(0) = dB_du_B2/w1 - B1 * dB_du_w2/(w1*w1);

        // Calculating  dB_dv
        Dyn_Vector dB_dv_B2 = this->p_weights.array() * deris.col(1).array();
        double dB_dv_w2 = dB_dv_B2.sum();
        deris.col(1) = dB_dv_B2/w1 - B1 * dB_dv_w2/(w1*w1);

        return deris;
    };



    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0), -1*sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 0
            GaussNode {LocalCoordinates(0,                -1*sqrt(3.0/5.0)), (8.0/9.0) * (5.0/9.0)}, // Node 1
            GaussNode {LocalCoordinates(sqrt(3.0/5.0),    -1*sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 2

            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0),                0), (5.0/9.0) * (8.0/9.0)}, // Node 3
            GaussNode {LocalCoordinates(0,                               0), (8.0/9.0) * (8.0/9.0)}, // Node 4
            GaussNode {LocalCoordinates(sqrt(3.0/5.0),                   0), (5.0/9.0) * (8.0/9.0)}, // Node 5

            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0),    sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 6
            GaussNode {LocalCoordinates(0,                   sqrt(3.0/5.0)), (8.0/9.0) * (5.0/9.0)}, // Node 7
            GaussNode {LocalCoordinates(sqrt(3.0/5.0),       sqrt(3.0/5.0)), (5.0/9.0) * (5.0/9.0)}, // Node 8

        };
        return gauss_nodes;
    }
};

}
