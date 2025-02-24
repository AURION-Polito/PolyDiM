#ifndef __VEM_DF_PCC_3D_LocalSpace_Data_HPP
#define __VEM_DF_PCC_3D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
struct VEM_DF_PCC_3D_Polyhedron_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;
    double Tolerance3D;

    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    std::vector<Eigen::MatrixXi> Faces;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::MatrixXd> TetrahedronVertices;

    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
    std::vector<Eigen::Vector3d> FacesNormal;
    std::vector<bool> FacesNormalDirection;
    std::vector<bool> FacesNormalGlobalDirection;
    std::vector<std::array<Eigen::Vector3d, 2>> FacesTangents;
    std::vector<std::array<bool, 2>> FacesTangentsGlobalDirection;

    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
};

struct VEM_DF_PCC_3D_Velocity_LocalSpace_Data final
{
    unsigned int Order;
    unsigned int Dimension;

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
    unsigned int NumBasisFunctions;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC BoundaryQuadrature;
    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC BoundaryQuadratureKL;

    double Diameter;
    double Measure;
    Eigen::Vector3d Centroid;

    Eigen::MatrixXd VanderInternal;
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives;

    Eigen::MatrixXd VanderBoundary;

    Eigen::MatrixXd VanderBoundaryKL;

    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives;

    std::vector<std::vector<Eigen::MatrixXd>> VectorDecompositionMatrices;
    std::vector<Eigen::MatrixXd> VanderGBigOPlus;
    std::vector<Eigen::MatrixXd> VanderGBigOPluskm2;

    Eigen::MatrixXd Hmatrix;
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
