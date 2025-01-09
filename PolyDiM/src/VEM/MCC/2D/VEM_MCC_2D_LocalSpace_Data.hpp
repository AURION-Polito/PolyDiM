#ifndef __VEM_MCC_2D_LocalSpace_Data_HPP
#define __VEM_MCC_2D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "QuadratureData.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
struct VEM_MCC_2D_Polygon_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;

    Eigen::MatrixXd Vertices;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::Matrix3d> TriangulationVertices;
    Eigen::VectorXd EdgesLength;
    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
    Eigen::MatrixXd EdgesNormal;
};

struct VEM_MCC_2D_Velocity_LocalSpace_Data final
{
    /// Order
    unsigned int Order;

    /// Geometric dimension
    unsigned int Dimension;

    /// Number of basis functions corresponding to degrees of freedom internal to edges (2D only).
    unsigned int NumBoundaryBasisFunctions;

    /// Number of basis functions corresponding to degrees of freedom internal to element (nabla).
    unsigned int NumNablaInternalBasisFunctions;

    /// Number of basis functions corresponding to degrees of freedom internal to element (bigoplus).
    unsigned int NumBigOPlusInternalBasisFunctions;
    unsigned int NumInternalBasisFunctions; ///< Number of basis functions corresponding to internal moments.
    unsigned int NumBasisFunctions;         ///< Number of basis functions.

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadrature;

    /// Vandermonde matrix of the polynomial basis at internal quadrature points.
    Eigen::MatrixXd VanderInternal;
    Eigen::MatrixXd VanderInternalKp1;

    std::vector<Eigen::MatrixXd> FacesVanderInternal;

    /// \brief Coefficients of basis functions on the reference edge.
    Eigen::VectorXd EdgeBasisCoefficients;

    Eigen::MatrixXd Pi0k; ///< Matrix representing the \f$\Pi^0_{\mathrm{order}}\f$ operator.

    Eigen::MatrixXd Wmatrix;

    double Diameter;
    double Measure;
    Eigen::Vector3d Centroid;

    /// \brief Mass matrix of the polynomial basis.
    /// \details \f$ H_{ij} = \int_E m_i m_j \f$, where \f$ E \f$ is the input polygon.
    Eigen::MatrixXd Hmatrix;               ///< mass matrix of order order
    Eigen::LLT<Eigen::MatrixXd> H_km1_LLT; ///< LLT factorization of the mass matrix of order-1 monomials.

    /// \brief Vandermonde matrix of the polynomial basis at boundary quadrature points.
    Eigen::MatrixXd VanderBoundary;
    Eigen::MatrixXd VanderBoundaryKp1;

    unsigned int Nk;
    unsigned int NkNabla;

    Eigen::MatrixXd TkNabla;
    Eigen::MatrixXd TkBigOPlus;
    Eigen::MatrixXd GkVanderInternal;
    Eigen::MatrixXd GkVanderBoundaryTimesNormal;

    Eigen::MatrixXd QmatrixKp1;
    Eigen::MatrixXd QmatrixInvKp1;
    Eigen::MatrixXd QmatrixTkNablaInv;
    Eigen::MatrixXd QmatrixTkNabla;

    Eigen::MatrixXd TkNablaDof;

    Eigen::MatrixXd Bmatrix;
    Eigen::MatrixXd Dmatrix;
    Eigen::MatrixXd Gmatrix;

    Eigen::MatrixXd Vmatrix;

    std::vector<Eigen::MatrixXd> VanderBasisFunctionValuesOnFace;
};

struct VEM_MCC_2D_Pressure_LocalSpace_Data final
{
    unsigned int Order;     ///< Order
    unsigned int Dimension; ///< Geometric dimension

    unsigned int Nk;
    unsigned int Nkm1;

    unsigned int NumBasisFunctions; ///< Number of basis functions.

    double Diameter;
    Eigen::Vector3d Centroid;

    Gedim::Quadrature::QuadratureData InternalQuadrature; ///< Internal quadrature points and weights

    Eigen::MatrixXd VanderInternal; /// Vandermonde matrix of the polynomial basis at internal quadrature points.

    Eigen::MatrixXd Hmatrix; ///< mass matrix of order order
};

} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
