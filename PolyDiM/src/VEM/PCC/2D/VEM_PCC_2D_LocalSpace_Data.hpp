#ifndef __VEM_PCC_2D_LocalSpace_Data_HPP
#define __VEM_PCC_2D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

/// \brief Structure containing the geometric properties of the element
struct VEM_PCC_2D_Polygon_Geometry final
{
    const double Tolerance1D;
    const double Tolerance2D;

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

struct VEM_PCC_2D_Inertia_Data final
{
    Eigen::MatrixXd Vertices;                           ///< cell2D vertices coordinates
    Eigen::MatrixXd OrderedVertices;                    ///< cell2D vertices coordinates
    Eigen::Vector3d Centroid;                           ///< cell2D centroids
    double Measure;                                     ///< cell2D areas
    double Diameter;                                    ///< cell2D diameters
    std::vector<Eigen::Matrix3d> TriangulationVertices; ///< cell2D triangulations
    Eigen::VectorXd EdgesLength;                        ///< cell2D edge lengths
    std::vector<bool> EdgesDirection;                   ///< cell2D edge directions
    Eigen::MatrixXd EdgesTangent;                       ///< cell2D edge tangents
    Eigen::MatrixXd EdgesNormal;                        ///< cell2D edge normals

    Eigen::Matrix3d Fmatrix;
    Eigen::Matrix3d FmatrixInv;
    Eigen::Vector3d translation;
    double absDetFmatrix;
    double signDetQ;

    double constantStiff;
    double constantMass;
};

/// \brief Structure containing the local matrices and the main variables to compute the vritual element discrete
/// matrices
struct VEM_PCC_2D_LocalSpace_Data final
{
    unsigned int Dimension;               ///< Geometrical dimension
    unsigned int Order;                   ///< Order of the space
    unsigned int NumVertexBasisFunctions; ///< Number of basis functions corresponding to degrees of freedom on
                                          ///< vertices.
    unsigned int NumEdgeBasisFunctions;   ///< Number of basis functions corresponding to degrees of freedom internal to
                                          ///< edges.
    unsigned int NumBoundaryBasisFunctions;  ///< Sum of \ref
                                             ///< VEM::PCC::VEM_PCC_2D_LocalSpace_Data::NumVertexBasisFunctions and \ref
                                             ///< VEM::PCC::VEM_PCC_2D_LocalSpace_Data::NumEdgeBasisFunctions.
    unsigned int NumInternalBasisFunctions;  ///< Number of basis functions corresponding to internal moments.
    unsigned int NumBasisFunctions;          ///< Number of basis functions.
    unsigned int NumEdgeDofs;                ///< Num of dofs on each edge
    unsigned int NumProjectorBasisFunctions; ///< dimension of the polynomial basis used for projectors.
    unsigned int Nkm1;                       ///< Dimension of the polynomial space of degree order-1.
    unsigned int Nkm2;                       ///< Dimension of the polynomial space of degree order-2.
    unsigned int Nklm1;

    Gedim::Quadrature::QuadratureData InternalQuadrature;                   ///< Internal quadrature points and weights
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadrature; ///< Boundary quadrature points and weights

    Gedim::Quadrature::QuadratureData InternalQuadratureKL; ///< Internal quadrature points and weights for e2
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadratureKL; ///< Boundary quadrature points and
                                                                              ///< weights for e2

    double Diameter;
    double Measure;
    Eigen::Vector3d Centroid;

    Eigen::MatrixXd VanderInternal;   ///< Vandermonde matrix of the polynomial basis at internal quadrature points.
    Eigen::MatrixXd VanderInternalKL; ///< Vandermonde matrix of the polynomial basis at internal quadrature points.
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives; ///< Vandermonde matrices of the derivatives of the
                                                            ///< polynomial basis at internal quadrature points.
    Eigen::MatrixXd VanderBoundary; ///< Vandermonde matrix of the polynomial basis at boundary quadrature points.
    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives; ///< Vandermonde matrices of the derivatives of the
                                                            ///< polynomial basis at boundary quadrature points.

    Eigen::MatrixXd PiNabla; ///< Matrix representing the \f$\Pi^\nabla_{\mathrm{order}}\f$ operator.
    Eigen::MatrixXd Pi0km1;  ///< Matrix representing the \f$\Pi^0_{\mathrm{order}-1}\f$ operator.
    Eigen::MatrixXd Pi0k;    ///< Matrix representing the \f$\Pi^0_{\mathrm{order}}\f$ operator.
    Eigen::MatrixXd Pi0klm1;
    std::vector<Eigen::MatrixXd> Pi0km1Der; ///< Vector of matrices representing the \f$\Pi^0_{\mathrm{order}-1}\f$
                                            ///< operator applied to derivatives of basis functions.

    Eigen::MatrixXd Hmatrix; ///< Mass matrix of the polynomial basis: \f$ H_{ij} = \int_E m_i m_j \f$, where \f$ E \f$
                             ///< is the input polygon.
    Eigen::MatrixXd H_klm1_matrix;
    Eigen::LLT<Eigen::MatrixXd> H_km1_LLT;  ///< LLT factorization of the mass matrix of order-1 monomials.
    Eigen::LLT<Eigen::MatrixXd> H_klm1_LLT; ///< LLT factorization of the mass matrix of order+l-1 monomials.
    Eigen::MatrixXd Cmatrix;                ///< C matrix
    Eigen::MatrixXd Bmatrix;                ///< B matrix
    Eigen::MatrixXd Gmatrix;                ///< G matrix
    Eigen::MatrixXd Dmatrix;                ///< D matrix
    std::vector<Eigen::MatrixXd> Ematrix;   ///< E matrices

    Eigen::MatrixXd Qmatrix;    ///< change of basis matrix: pV = mV*Q'
    Eigen::MatrixXd QmatrixInv; ///< inverse of \a Qmatrix
    Eigen::MatrixXd Qmatrixkm1; ///< change of basis matrix of order (order-1)

    VEM_PCC_2D_Inertia_Data inertia_data;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
