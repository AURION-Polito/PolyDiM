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

struct VEM_PCC_2D_Polygon_Geometry final
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

struct VEM_PCC_2D_Inertia_Data final
{
    Eigen::MatrixXd Vertices;
    Eigen::MatrixXd OrderedVertices;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::Matrix3d> TriangulationVertices;
    Eigen::VectorXd EdgesLength;
    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
    Eigen::MatrixXd EdgesNormal;

    Eigen::Matrix3d Fmatrix;
    Eigen::Matrix3d FmatrixInv;
    Eigen::Vector3d translation;
    double absDetFmatrix;
    double signDetQ;

    double constantStiff;
    double constantMass;
};

struct VEM_PCC_2D_LocalSpace_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumVertexBasisFunctions;

    unsigned int NumEdgeBasisFunctions;

    unsigned int NumBoundaryBasisFunctions;

    unsigned int NumInternalBasisFunctions;
    unsigned int NumBasisFunctions;
    unsigned int NumEdgeDofs;
    unsigned int NumProjectorBasisFunctions;
    unsigned int Nkm1;
    unsigned int Nkm2;
    unsigned int Nklm1;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadrature;

    Gedim::Quadrature::QuadratureData InternalQuadratureKL;
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadratureKL;

    double Diameter;
    double Measure;
    Eigen::Vector3d Centroid;

    Eigen::MatrixXd VanderInternal;
    Eigen::MatrixXd VanderInternalKL;
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives;

    Eigen::MatrixXd VanderBoundary;
    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives;

    Eigen::MatrixXd PiNabla;
    Eigen::MatrixXd Pi0km1;
    Eigen::MatrixXd Pi0k;
    Eigen::MatrixXd Pi0klm1;
    std::vector<Eigen::MatrixXd> Pi0km1Der;

    Eigen::MatrixXd Hmatrix;

    Eigen::MatrixXd H_klm1_matrix;
    Eigen::LLT<Eigen::MatrixXd> H_km1_LLT;
    Eigen::LLT<Eigen::MatrixXd> H_klm1_LLT;
    Eigen::MatrixXd Cmatrix;
    Eigen::MatrixXd Bmatrix;
    Eigen::MatrixXd Gmatrix;
    Eigen::MatrixXd Dmatrix;
    std::vector<Eigen::MatrixXd> Ematrix;

    Eigen::MatrixXd Qmatrix;
    Eigen::MatrixXd QmatrixInv;
    Eigen::MatrixXd Qmatrixkm1;

    VEM_PCC_2D_Inertia_Data inertia_data;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
