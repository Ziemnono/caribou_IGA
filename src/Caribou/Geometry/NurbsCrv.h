#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseNurbsCrv.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct NurbsCrv;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<NurbsCrv <_Dimension>> {
     // static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; // NurbsCrv degree
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 1; // 1 parametric direction.
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 3; // NurbsCrv control points (p+1)*(q+1)
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 3; // Gauss points per element
 };

 /**
  * Bezier Surface two canonicalDimensions
  *
  * \verbatim
  *       --> u
  *  0-----1-----2
  *  \endverbatim
  *
  * @tparam _Dimension The world coordinates dimension
  */

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct NurbsCrv: public BaseNurbsCrv<NurbsCrv <_Dimension>> {
    // Types

    using Scalar = FLOATING_POINT_TYPE;

    using Base = BaseNurbsCrv<NurbsCrv<_Dimension>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    using Dyn_Matrix = typename Base::Dyn_Matrix;
    using Dyn_Vector = typename Base::Dyn_Vector;

    static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; //NurbsCrv
    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime; //$$$NurbsCrv
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime; //$$$NurbsCrv

    // Constructors
    using Base::Base;
    NurbsCrv() : Base() {
        if constexpr (Dimension == 1) {
            this->p_nodes[0] = -1;
            this->p_nodes[1] = +0;
            this->p_nodes[2] = +1;
        } else if constexpr (Dimension == 2) {
            // WorldCoordinates(x, y, w)
            this->p_nodes.row(0) = WorldCoordinates(-1, 0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates( 0, 0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates( 1, 0); // Node 2

        } else {
            // WorldCoordinates(x, y, z, w)
            this->p_nodes.row(0) = WorldCoordinates(-1, 0, 0); // Node 0
            this->p_nodes.row(1) = WorldCoordinates( 0, 0, 0); // Node 1
            this->p_nodes.row(2) = WorldCoordinates( 1, 0, 0); // Node 2

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
        std::cout << "Knot Vector  : \n" << this->p_knots.transpose() << "\n";
        std::cout << "Knot span      : \n" << this->p_knot_span.transpose() << "\n";

    }

private:
    // Implementations
    friend struct Element<NurbsCrv <_Dimension>>;
    friend struct BaseNurbsCrv<NurbsCrv <_Dimension>>;


    auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {
        // U direction basis
        Scalar u = 0.5 * (xi[0]*(this->p_knot_span[1] - this->p_knot_span[0]) + (this->p_knot_span[1] + this->p_knot_span[0]));

        auto knots = this->p_knots;
        Vector<NumberOfNodesAtCompileTime> basis = bspbasisfun(u, Degree, knots);

        basis = basis.array() * this->p_weights.array();
        basis = basis * (1/basis.sum());
        return basis;

    };

    inline auto get_dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {
        Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> dL;
        Scalar u = 0.5 * (xi[0] * (this->p_knot_span[1] - this->p_knot_span[0]) + (this->p_knot_span[1] + this->p_knot_span[0]));

//        std::cout << "\nu - " <<  u << "\n";
//        std::cout << "v - " <<  v << "\n";
        auto knots = this->p_knots;

        Dyn_Vector N(Degree+1);
        Dyn_Vector Nderis(Degree+1);

        N = bspbasisfun(u, Degree, knots);
        Nderis = basis_deri(u, Degree, knots);

        Scalar w = (N.array() * this->p_weights.array()).sum();
        Scalar dwdxi = (Nderis.array() * this->p_weights.array()).sum();

        Dyn_Vector fac = this->p_weights.array()/(w*w);
        N = N.array() * fac.array() * w;
        Nderis = (Nderis*w - N*dwdxi).array() * fac.array();
//        for(int i = 0; i <= Degree; i++)
//        {


//            fac      = this->p_weights[i]/(w*w);
//            N[i]     = N[i]*fac*w;
//            Nderis[i] = (Nderis[i]*w - N[i]*dwdxi) * fac;
//        }
        return Nderis;
    };



    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0)),  (5.0/9.0)}, // Node 0

            GaussNode {LocalCoordinates(0), (8.0/9.0)}, // Node 1

            GaussNode {LocalCoordinates(sqrt(3.0/5.0)), (5.0/9.0)}, // Node 2

        };
        return gauss_nodes;
    }
};

}
