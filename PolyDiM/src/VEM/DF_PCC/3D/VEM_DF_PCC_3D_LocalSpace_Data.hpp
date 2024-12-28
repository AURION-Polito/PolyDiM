#ifndef __VEM_DF_PCC_3D_LocalSpace_Data_HPP
#define __VEM_DF_PCC_3D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "GeometryUtilities.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

/// \brief Structure containing the geometric properties of the element
struct VEM_DF_PCC_3D_Polyhedron_Geometry final
{
    const double Tolerance1D;
    const double Tolerance2D;
    const double Tolerance3D;

    const Eigen::MatrixXd &Vertices;
    const Eigen::MatrixXi &Edges;
    const std::vector<Eigen::MatrixXi> &Faces;
    const Eigen::Vector3d &Centroid;
    const double &Measure;
    const double &Diameter;
    const std::vector<Eigen::MatrixXd> &TetrahedronVertices;

    const std::vector<Eigen::Matrix3d> &FacesRotationMatrix;
    const std::vector<Eigen::Vector3d> &FacesTranslation;
    const std::vector<Eigen::Vector3d> &FacesNormal;
    const std::vector<bool> &FacesNormalDirection;
    const std::vector<bool> &FacesNormalGlobalDirection;
    const std::vector<std::array<Eigen::Vector3d, 2>> &FacesTangents;
    const std::vector<std::array<bool, 2>> &FacesTangentsGlobalDirection;

    const std::vector<bool> &EdgesDirection;
    const Eigen::MatrixXd &EdgesTangent;
};

/// \brief Structure containing the local matrices and the main variables to compute the vritual element discrete
/// matrices
struct VEM_DF_PCC_3D_Velocity_LocalSpace_Data final
{
    unsigned int Order;     ///< Order
    unsigned int Dimension; ///< Geometric dimension

    unsigned int NKp1;
    unsigned int Nk;
    unsigned int Nkm1;
    unsigned int Nkm2;
    unsigned int Nkm3;
    unsigned int Nkm4;

    std::vector<PCC::VEM_PCC_2D_LocalSpace_Data> facesLocalSpace;

    unsigned int NumVertexBasisFunctions;
    unsigned int NumEdgeBasisFunctions;
    unsigned int NumFaceBasisFunctions;
    unsigned int NumNormalBasisFunctions;
    unsigned int NumTangentsBasisFunctions;
    unsigned int NumBoundaryBasisFunctions;
    unsigned int NumDivergenceInternalBasisFunctions;
    unsigned int NumBigOPlusInternalBasisFunctions;
    unsigned int NumInternalBasisFunctions;
    unsigned int NumBasisFunctions; ///< Number of basis functions.

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC BoundaryQuadrature;
    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC BoundaryQuadratureKL;

    double Diameter;
    double Measure;
    Eigen::Vector3d Centroid;

    Eigen::MatrixXd VanderInternal; ///< Vandermonde matrix of the polynomial basis at internal quadrature points.
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives; ///< Vandermonde matrices of the derivatives of the
                                                            ///< polynomial basis at internal quadrature points.
    Eigen::MatrixXd VanderBoundary; ///< \brief Vandermonde matrix of the polynomial basis at boundary quadrature
                                    ///< points.
    Eigen::MatrixXd VanderBoundaryKL;

    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives;

    std::vector<std::vector<Eigen::MatrixXd>> VectorDecompositionMatrices;
    std::vector<Eigen::MatrixXd> VanderGBigOPlus;
    std::vector<Eigen::MatrixXd> VanderGBigOPluskm2;

    /// \brief Mass matrix of the polynomial basis.
    /// \details \f$ H_{ij} = \int_E m_i m_j \f$, where \f$ E \f$ is the input polygon.
    Eigen::MatrixXd Hmatrix; ///< mass matrix of order order
    Eigen::MatrixXd HmatrixKp1;

    std::vector<Eigen::MatrixXd> VanderFaceProjectionsKm1;
    Eigen::MatrixXd VanderEdgeDofs;
    Eigen::MatrixXd VanderFaceProjectionsKp1TimesNormal;
    std::vector<Eigen::MatrixXd> FaceScaledMomentsBasis;
    std::vector<Eigen::MatrixXd> ScaledHmatrixOnBoundary;

    std::vector<Eigen::MatrixXd> PiNabla;
    std::vector<Eigen::MatrixXd> Pi0km2;
    std::vector<Eigen::MatrixXd> Pi0k;
    std::vector<Eigen::MatrixXd> Pi0km1Der;

    Eigen::MatrixXd Wmatrix;
    Eigen::MatrixXd Vmatrix;
    std::vector<Eigen::MatrixXd> Bmatrix;
    std::vector<Eigen::MatrixXd> Dmatrix;
    Eigen::MatrixXd Gmatrix;
    std::vector<Eigen::MatrixXd> Cmatrix;
    std::vector<Eigen::MatrixXd> Cmatrixkm2;
    std::vector<Eigen::MatrixXd> Ematrix;
};

struct VEM_DF_PCC_3D_Pressure_LocalSpace_Data final
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

} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
