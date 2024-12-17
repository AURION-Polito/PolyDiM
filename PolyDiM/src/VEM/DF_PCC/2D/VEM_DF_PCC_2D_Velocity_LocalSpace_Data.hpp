#ifndef __VEM_DF_PCC_2D_Velocity_LocalSpace_Data_HPP
#define __VEM_DF_PCC_2D_Velocity_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim {
namespace VEM {
namespace DF_PCC {

/// \brief Structure containing the geometric properties of the element
struct VEM_DF_PCC_2D_Polygon_Geometry final
{
    const Eigen::MatrixXd &Vertices;
    const Eigen::Vector3d &Centroid;
    const double &Measure;
    const double &Diameter;
    const std::vector<Eigen::Matrix3d> &TriangulationVertices;
    const Eigen::VectorXd &EdgesLength;
    const std::vector<bool> &EdgesDirection;
    const Eigen::MatrixXd &EdgesTangent;
    const Eigen::MatrixXd &EdgesNormal;
};

/// \brief Structure containing the local matrices and the main variables to compute the vritual element discrete matrices
struct VEM_DF_PCC_2D_Velocity_LocalSpace_Data final
{

    unsigned int Order; ///< Order
    unsigned int Dimension; ///< Geometric dimension

    unsigned int NKp1;
    unsigned int Nk;
    unsigned int Nkm1;
    unsigned int Nkm2;
    unsigned int Nkm3;
    unsigned int Nkm4;

    unsigned int NumVertexBasisFunctions;
    unsigned int NumEdgeBasisFunctions;

    unsigned int NumBoundaryBasisFunctions;

    unsigned int NumDivergenceInternalBasisFunctions;
    unsigned int NumBigOPlusInternalBasisFunctions;

    unsigned int NumInternalBasisFunctions;

    unsigned int NumBasisFunctions; ///< Number of basis functions.

    Gedim::Quadrature::QuadratureData InternalQuadrature; ///< Internal quadrature points and weights
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadrature; ///< Boundary quadrature points and weights
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData EdgesDOFs; ///< Boundary quadrature points and weights


    Eigen::MatrixXd VanderInternal;  /// Vandermonde matrix of the polynomial basis at internal quadrature points.
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives; /// Vandermonde matrices of the derivatives of the polynomial basis at internal quadrature points.
    Eigen::MatrixXd VanderBoundary; ///< Vandermonde matrix of the polynomial basis at boundary quadrature points.
    Eigen::MatrixXd VanderBoundaryKp1;
    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives;
    Eigen::MatrixXd VanderKp1EdgeProjections;

    std::vector<std::vector<Eigen::MatrixXd>> VectorDecompositionMatrices;
    std::vector<Eigen::MatrixXd> VanderGBigOPlus;
    std::vector<Eigen::MatrixXd> VanderGBigOPluskm2;

    /// \brief Mass matrix of the polynomial basis.
    /// \details \f$ H_{ij} = \int_E m_i m_j \f$, where \f$ E \f$ is the input polygon.
    Eigen::MatrixXd Hmatrix; ///< mass matrix of order order
    Eigen::MatrixXd HmatrixKp1;

    Eigen::VectorXd ReferenceEdgeInternalPoints;
    Eigen::VectorXd ReferenceEdgeDofInternalPoints;
    Eigen::VectorXd EdgeBasisCoefficients;

    std::vector<Eigen::MatrixXd> PiNabla;
    std::vector<Eigen::MatrixXd> Pi0km2;
    std::vector<Eigen::MatrixXd> Pi0k;
    std::vector<Eigen::MatrixXd> Pi0km1Der;

    Eigen::MatrixXd StabMatrix; ///< Matrix used for stabilizing elliptic bilinear forms.

    Eigen::MatrixXd Wmatrix;
    Eigen::MatrixXd Vmatrix;
    std::vector<Eigen::MatrixXd> Bmatrix;
    std::vector<Eigen::MatrixXd> Dmatrix;
    Eigen::MatrixXd Gmatrix;
    std::vector<Eigen::MatrixXd> Cmatrix;
    std::vector<Eigen::MatrixXd> Cmatrixkm2;
    std::vector<Eigen::MatrixXd> Ematrix;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
