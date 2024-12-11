#ifndef __VEM_MCC_2D_Velocity_LocalSpace_Data_HPP
#define __VEM_MCC_2D_Velocity_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"
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
    const Eigen::MatrixXd& Vertices;
    const Eigen::Vector3d& Centroid;
    const double& Measure;
    const double& Diameter;
    const std::vector<Eigen::Matrix3d>& TriangulationVertices;
    const Eigen::VectorXd& EdgesLength;
    const std::vector<bool>& EdgesDirection;
    const Eigen::MatrixXd& EdgesTangent;
    const Eigen::MatrixXd& EdgesNormal;
};

struct VEM_MCC_3D_Polyhedron_Geometry final
{
    const Gedim::GeometryUtilities& GeometryUtility;

    const Eigen::MatrixXd& Vertices;
    const Eigen::Vector3d& Centroid;
    const double& Measure;
    const double& Diameter;
    const std::vector<Eigen::MatrixXd>& TetrahedronVertices;

    const std::vector<Eigen::Matrix3d>& FacesRotationMatrix;
    const std::vector<Eigen::Vector3d>& FacesTranslation;
    const std::vector<Eigen::Vector3d>& FacesNormal;
    const std::vector<bool>& FacesNormalDirection;
    const std::vector<bool>& FacesGlobalNormalDirection;
    const std::vector<double>& FacesMeasure;
    const std::vector<Eigen::Vector3d>& FacesCentroid2D;
    const std::vector<double>& FacesDiameter;
    const std::vector<std::vector<Eigen::Matrix3d>>& FacesTriangulationVertices2D;
};

struct VEM_MCC_Velocity_LocalSpace_Data final
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
    unsigned int NumInternalBasisFunctions;  ///< Number of basis functions corresponding to internal moments.
    unsigned int NumBasisFunctions; ///< Number of basis functions.

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

    /// \brief Mass matrix of the polynomial basis.
    /// \details \f$ H_{ij} = \int_E m_i m_j \f$, where \f$ E \f$ is the input polygon.
    Eigen::MatrixXd Hmatrix; ///< mass matrix of order order
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
    Eigen::MatrixXd StabMatrix;

    Eigen::MatrixXd Vmatrix;

    std::vector<Eigen::MatrixXd> VanderBasisFunctionValuesOnFace;
};
}
}
}

#endif
