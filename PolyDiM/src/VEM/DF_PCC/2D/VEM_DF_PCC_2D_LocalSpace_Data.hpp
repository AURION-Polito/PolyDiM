#ifndef __VEM_DF_PCC_2D_LocalSpace_Data_HPP
#define __VEM_DF_PCC_2D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

struct VEM_DF_PCC_2D_Polygon_Geometry final
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

struct VEM_DF_PCC_2D_Velocity_LocalSpace_Data final
{
    unsigned int Order;
    unsigned int Dimension;

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

    unsigned int NumBasisFunctions;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadrature;
    Quadrature::VEM_Quadrature_2D::Edges_QuadratureData EdgesDOFs;

    double Diameter;
    double Measure;
    Eigen::Vector3d Centroid;

    Eigen::MatrixXd VanderInternal;
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives;

    Eigen::MatrixXd VanderBoundary;
    Eigen::MatrixXd VanderBoundaryKp1;
    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives;
    Eigen::MatrixXd VanderKp1EdgeProjections;

    std::vector<std::vector<Eigen::MatrixXd>> VectorDecompositionMatrices;
    std::vector<Eigen::MatrixXd> VanderGBigOPlus;
    std::vector<Eigen::MatrixXd> VanderGBigOPluskm2;

    Eigen::MatrixXd Hmatrix;
    Eigen::MatrixXd HmatrixKp1;

    Eigen::VectorXd ReferenceEdgeInternalPoints;
    Eigen::VectorXd ReferenceEdgeDofInternalPoints;
    Eigen::VectorXd EdgeBasisCoefficients;

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

struct VEM_DF_PCC_2D_Pressure_LocalSpace_Data final
{
    unsigned int Order;
    unsigned int Dimension;

    unsigned int Nk;
    unsigned int Nkm1;

    unsigned int NumBasisFunctions;

    double Diameter;
    Eigen::Vector3d Centroid;

    Gedim::Quadrature::QuadratureData InternalQuadrature;

    Eigen::MatrixXd VanderInternal;

    Eigen::MatrixXd Hmatrix;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
