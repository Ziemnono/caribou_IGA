#pragma once

#include <Caribou/config.h>
#include <Caribou/Geometry/Base/BaseSpline.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Spline;

template<UNSIGNED_INTEGER_TYPE _Dimension>
 struct traits<Spline <_Dimension>> {
     // static constexpr UNSIGNED_INTEGER_TYPE Degree = 2; // Spline degree
     static constexpr UNSIGNED_INTEGER_TYPE CanonicalDimension = 1; // Single parametric direction.
     static constexpr UNSIGNED_INTEGER_TYPE Dimension = _Dimension;
     static constexpr INTEGER_TYPE NumberOfNodesAtCompileTime = 3; // Spline control points
     static constexpr INTEGER_TYPE NumberOfGaussNodesAtCompileTime = 3; // Gauss points per element
 };

/**
 * Single element Spline
 *
 * \verbatim
 * P1 : 0-----+-----1 --> u
 * \endverbatim
 *
 * @tparam _Dimension The world coordinates dimension
 */
template<UNSIGNED_INTEGER_TYPE _Dimension>
struct Spline: public BaseSpline<Spline <_Dimension>> {
    // Types

    using Scalar = FLOATING_POINT_TYPE;

    using Base = BaseSpline<Spline <_Dimension>>;
    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime; //$$$Spline
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime; //$$$Spline

    static constexpr auto Degree = 2; //Spline
//    static constexpr UNSIGNED_INTEGER_TYPE NumberOfCtrlPnts = Degree+1; // No of control points.
    static constexpr UNSIGNED_INTEGER_TYPE NumBasis = Degree+1; // Nmber of basis per element.
    const  Vector<NumberOfNodesAtCompileTime+Degree+1> KnotVector = open_knot_vect();

    static constexpr Scalar XiPara1 = 0.0; //KnotVector[2]
    static constexpr Scalar XiPara2 = 1.0; //KnotVector[3]


    // Constructors
    using Base::Base;
    Spline() : Base() {
        if constexpr (Dimension == 1) {
            this->p_nodes[0] = -1; //
            this->p_nodes[1] = 0;  //
            this->p_nodes[2] = +1; //
        } else if constexpr (Dimension == 2) {
            this->p_nodes.row(0) = WorldCoordinates(-1, 0);
            this->p_nodes.row(1) = WorldCoordinates(0, 0); //
            this->p_nodes.row(2) = WorldCoordinates(+1, 0); //
        } else {
            this->p_nodes.row(0) = WorldCoordinates(-1, 0, 0);
            this->p_nodes.row(1) = WorldCoordinates(0, 0, 0); //
            this->p_nodes.row(2) = WorldCoordinates(+1, 0, 0); //
        }
    };


private:
    // Implementations
    friend struct Element<Spline <_Dimension>>;
    friend struct BaseSpline<Spline <_Dimension>>;

    inline auto findspan(const Scalar& u) const -> UNSIGNED_INTEGER_TYPE{
        const auto & n = NumberOfNodesAtCompileTime;
        const Vector<NumberOfNodesAtCompileTime+Degree+1> & U = KnotVector;

        if (u == U[n+2])
        {
            return n;
        }
        UNSIGNED_INTEGER_TYPE low, high, mid;
        low = Degree;
        high = n+1;
        mid = (low + high)/2;
        while (u < U[mid] || u >= U[mid+1])
        {
            if (u < U[mid]){high = mid;}
            else {low = mid;}
            mid = (low+high)/2;
        }

        return mid;
    }

    inline auto open_knot_vect() const -> Vector<NumberOfNodesAtCompileTime+Degree+1> {
        Vector<NumberOfNodesAtCompileTime+Degree+1> U;
        for (int i=0; i <= NumberOfNodesAtCompileTime+Degree; i++)
        {
            if (i<=Degree){U[i] = 0;}
            else {U[i] = 1;}
        };
        return U;
    }

    inline auto get_L(const LocalCoordinates & xi) const -> Vector<NumberOfNodesAtCompileTime> {

        const Scalar u = 0.5 * (xi[0]*(XiPara2 - XiPara1) + (XiPara2 - XiPara1));
        const UNSIGNED_INTEGER_TYPE i = findspan(u);
        const Vector<NumberOfNodesAtCompileTime+Degree+1> & U = KnotVector; // U is alias of KnotVector

        Vector<NumBasis> N;
        Vector<NumBasis> left;
        Vector<NumBasis> right;
        N[0] = 1;
        Scalar saved, temp;
        for (int j = 1; j <= Degree; j++){
            left[j] = u - U[i+1-j];
            right[j] = U[i+j] - u;
            saved = 0.0;
            for (int r = 0; r<j; r++){
                temp = N[r]/(right[r+1]+left[j-r]);
                N[r] = saved + right[r+1]*temp;
                saved = left[j-r]*temp;
            }
            N[j] = saved;
        }
        return N;

//        return {
//            static_cast<FLOATING_POINT_TYPE>(1/2. * (1 - u)),
//            static_cast<FLOATING_POINT_TYPE>(1/2. * (1 + u))
//        };
    };

    inline auto get_dL(const LocalCoordinates & xi) const -> Matrix<NumberOfNodesAtCompileTime, CanonicalDimension> {

        const auto  & u = 0.5 * ((XiPara2 - XiPara1) * xi[0] + (XiPara2 - XiPara1));

        // Copying the knotvector to U for easy usage.
        const Vector<NumberOfNodesAtCompileTime+Degree+1> & U = KnotVector;
        static constexpr int n = 1; // Derivative order.
        const UNSIGNED_INTEGER_TYPE i = findspan(u);
        Matrix<NumBasis, NumBasis> ndu;
        Matrix<n+1, NumBasis> ders;
        Matrix<2, NumBasis> a;
        ndu(0,0) = 1.0;

        Vector<NumBasis> left;
        Vector<NumBasis> right;
        Scalar saved, temp;
        int j;
        for (j = 1; j <= Degree; j++)
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

        for (j = 0; j <= Degree; j++)
        {
            ders(0,j) = ndu(j,Degree);
        }

        int s1, s2, rk, pk, j1, j2;
        Scalar d;
        for (int r = 0; r <= Degree; r++)
        {
            s1 = 0; s2 = 1;
            a(0,0) = 1.0;
            for (int k = 1; k <=n ; k++)
            {
                d = 0.0;
                rk = r-k;
                pk = Degree - k;
                if (r >= k)
                {
                    a(s2, 0) = a(s1,0)/ndu(pk+1, rk);
                    d = a(s2, 0)*ndu(rk, pk);
                }
                if (rk >= -1)
                {
                    j1 = 1;
                }
                else
                {
                    j1 = -rk;
                }

                if ((r-1) <= pk)
                {
                    j2 = k-1;
                }
                else
                {
                    j2 = Degree-r;
                }

                for (j = j1; j <= j2; j++)
                {
                    a(s2, j) = (a(s1, j)-a(s1, j-1))/ndu(pk+1, rk+j);
                    d = d + a(s2, j) * ndu(rk+j, pk);
                }

                if (r <= pk)
                {
                    a(s2, k) = -1*a(s1, k-1)/ndu(pk+1, r);
                    d = d + a(s2, k) * ndu(r, pk);
                }
                ders(k, r) = d;
                j = s1;
                s1 = s2;
                s2 = j;
            }
        }
        int r = Degree;
        for (int k = 1; k <= n; k++)
        {
            for (j = 0; j <= Degree; j++)
            {
                ders(k, j) = ders(k, j) * r;
            }
            r = r * (Degree-k);
        }
        auto para_natural_jac = 0.5 * (XiPara2 - XiPara1);
        auto dl = ders.row(1)*para_natural_jac;
//        return ders.row(1)*para_natural_jac;
        return {
            static_cast<FLOATING_POINT_TYPE>(dl[0]),
            static_cast<FLOATING_POINT_TYPE>(dl[1]),
            static_cast<FLOATING_POINT_TYPE>(dl[2])
        };
    };



    inline auto get_gauss_nodes() const -> const auto & {
        static std::vector<GaussNode> gauss_nodes {
            GaussNode {LocalCoordinates(-1*sqrt(3.0/5.0)), 5.0/9.0}, // Node 0
            GaussNode {LocalCoordinates(0)             , 8.0/9.0}, // Node 1
            GaussNode {LocalCoordinates(sqrt(3.0/5.0)) , 5.0/9.0}, // Node 2
        };
        return gauss_nodes;
    }
};

}
