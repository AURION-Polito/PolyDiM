#ifndef __VEM_PCC_3D_LocalSpace_Data_HPP
#define __VEM_PCC_3D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

struct VEM_PCC_3D_Polyhedron_Geometry final
{
    const Gedim::GeometryUtilities &GeometryUtility;

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

    const std::vector<bool> &EdgesDirection;
    const Eigen::MatrixXd &EdgesTangent;
};

struct VEM_PCC_3D_LocalSpace_Data final
{
    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC BoundaryQuadrature;

    std::vector<VEM_PCC_2D_LocalSpace_Data> facesLocalSpace;

    /// Geometrical dimension
    unsigned int Dimension;
    /// Order of the space
    unsigned int Order;

    /// Number of basis functions corresponding to degrees of freedom on vertices.
    unsigned int NumVertexBasisFunctions;
    /// Number of basis functions corresponding to degrees of freedom internal to edges.
    unsigned int NumEdgeBasisFunctions;
    /// Number of basis functions corresponding to degrees of freedom internal to faces.
    unsigned int NumFaceBasisFunctions;
    /// Number of basis functions corresponding to internal moments.
    unsigned int NumInternalBasisFunctions;
    unsigned int NumBasisFunctions; ///< Number of basis functions.
    unsigned int NumEdgeDofs;       ///< Num of dofs on each edge
    unsigned int NumFaceDofs;       ///< Num of dofs on each face
    unsigned int NumBoundaryBasisFunctions;
    unsigned int NumProjectorBasisFunctions; ///< dimension of the polynomial basis used for projectors.

    unsigned int Nkm1; ///< Dimension of the polynomial space of degree order-1.
    unsigned int Nkm2; ///< Dimension of the polynomial space of degree order-2.

    /// Vandermonde matrix of the polynomial basis at internal quadrature points.
    Eigen::MatrixXd VanderInternal;
    /// Vandermonde matrices of the derivatives of the polynomial basis at internal quadrature
    /// points.
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives;

    /// \brief Coefficients of basis functions on the reference edge.
    Eigen::VectorXd EdgeBasisCoefficients;
    Eigen::MatrixXd PiNabla; ///< Matrix representing the \f$\Pi^\nabla_{\mathrm{order}}\f$ operator.
    Eigen::MatrixXd Pi0km1;  ///< Matrix representing the \f$\Pi^0_{\mathrm{order}-1}\f$ operator.
    Eigen::MatrixXd Pi0k;    ///< Matrix representing the \f$\Pi^0_{\mathrm{order}}\f$ operator.

    /// Vector of matrices representing the \f$\Pi^0_{\mathrm{order}-1}\f$ operator applied to
    /// derivatives of basis functions.
    std::vector<Eigen::MatrixXd> Pi0km1Der;
    Eigen::MatrixXd StabMatrix;     ///< Matrix used for stabilizing elliptic bilinear forms.
    Eigen::MatrixXd StabMatrixPi0k; ///< Matrix used for stabilizing elliptic bilinear forms.
    /// \brief Mass matrix of the polynomial basis.
    /// \details \f$ H_{ij} = \int_E m_i m_j \f$, where \f$ E \f$ is the input polygon.
    Eigen::MatrixXd Hmatrix;               ///< mass matrix of order order
    Eigen::LLT<Eigen::MatrixXd> H_km1_LLT; ///< LLT factorization of the mass matrix of order-1 monomials.
    Eigen::MatrixXd Cmatrix;               ///< C matrix
    Eigen::MatrixXd Bmatrix;               ///< B matrix
    Eigen::MatrixXd Gmatrix;               ///< G matrix
    Eigen::MatrixXd Dmatrix;               ///< G matrix
    std::vector<Eigen::MatrixXd> Ematrix;

    /// \brief Vandermonde matrix of the polynomial basis at boundary quadrature points.
    Eigen::MatrixXd VanderBoundary;

    Eigen::MatrixXd VanderEdgeDofs;
    Eigen::MatrixXd ScaledHmatrixOnBoundary;

    Eigen::MatrixXd VanderFaceProjections;
    std::vector<Eigen::MatrixXd> FaceScaledMomentsBasis;
    Eigen::MatrixXd PointEdgeDofsCoordinates;

    /// \brief Vandermonde matrices of the derivatives of the polynomial basis at boundary
    /// quadrature points.
    /// \sa \ref vanderBoundary.
    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives;

    Eigen::MatrixXd Qmatrix; // change of basis matrix: pV = mV*Q'
    Eigen::MatrixXd QmatrixInv;
    Eigen::MatrixXd Qmatrixkm1;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
