// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __VEM_MCC_3D_LocalSpace_Data_HPP
#define __VEM_MCC_3D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "QuadratureData.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{

struct VEM_MCC_3D_Polyhedron_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;
    double Tolerance3D;

    Eigen::MatrixXd Vertices;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::MatrixXd> TetrahedronVertices;

    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
    std::vector<Eigen::Vector3d> FacesNormal;
    std::vector<bool> FacesNormalDirection;
    std::vector<bool> FacesGlobalNormalDirection;
    std::vector<double> FacesMeasure;
    std::vector<Eigen::Vector3d> FacesCentroid2D;
    std::vector<double> FacesDiameter;
    std::vector<std::vector<Eigen::Matrix3d>> FacesTriangulationVertices2D;
};

struct VEM_MCC_3D_Velocity_LocalSpace_Data final
{
    unsigned int Order;

    unsigned int Dimension;

    unsigned int NumBoundaryBasisFunctions;

    unsigned int NumNablaInternalBasisFunctions;

    unsigned int NumBigOPlusInternalBasisFunctions;
    unsigned int NumInternalBasisFunctions;
    unsigned int NumBasisFunctions;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_MCC BoundaryQuadrature;

    Eigen::MatrixXd VanderInternal;
    Eigen::MatrixXd VanderInternalKp1;

    std::vector<Eigen::MatrixXd> FacesVanderInternal;

    Eigen::VectorXd EdgeBasisCoefficients;

    Eigen::MatrixXd Pi0k;

    Eigen::MatrixXd Wmatrix;

    Eigen::MatrixXd Hmatrix;
    Eigen::LLT<Eigen::MatrixXd> H_km1_LLT;

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

    double Diameter;
    double Measure;
    Eigen::Vector3d Centroid;

    Eigen::MatrixXd TkNablaDof;

    Eigen::MatrixXd Bmatrix;
    Eigen::MatrixXd Dmatrix;
    Eigen::MatrixXd Gmatrix;

    Eigen::MatrixXd Vmatrix;

    std::vector<Eigen::MatrixXd> VanderBasisFunctionValuesOnFace;
};

struct VEM_MCC_3D_Pressure_LocalSpace_Data final
{
    unsigned int Order;
    unsigned int Dimension;

    unsigned int Nk;
    unsigned int Nkm1;

    unsigned int NumBasisFunctions;

    double Diameter;
    Eigen::Vector3d Centroid;

    Gedim::Quadrature::QuadratureData InternalQuadrature;

    Eigen::MatrixXd Qmatrix;
    Eigen::MatrixXd VanderInternal;

    Eigen::MatrixXd Hmatrix;
};

} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
