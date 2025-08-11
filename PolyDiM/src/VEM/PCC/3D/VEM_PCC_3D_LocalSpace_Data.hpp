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

    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
};

struct VEM_PCC_3D_LocalSpace_Data final
{
    Gedim::Quadrature::QuadratureData InternalQuadrature;
    Quadrature::VEM_Quadrature_3D::Faces_QuadratureData_PCC BoundaryQuadrature;

    std::vector<VEM_PCC_2D_LocalSpace_Data> facesLocalSpace;

    unsigned int Dimension;
    unsigned int Order;

    unsigned int NumVertexBasisFunctions;
    unsigned int NumEdgeBasisFunctions;
    unsigned int NumFaceBasisFunctions;
    unsigned int NumInternalBasisFunctions;
    unsigned int NumBasisFunctions;
    unsigned int NumEdgeDofs;
    unsigned int NumFaceDofs;
    unsigned int NumBoundaryBasisFunctions;
    unsigned int NumProjectorBasisFunctions;

    unsigned int Nkm1;
    unsigned int Nkm2;

    Eigen::MatrixXd VanderInternal;
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives;

    double Diameter;
    Eigen::Vector3d Centroid;
    double Measure;

    Eigen::VectorXd EdgeBasisCoefficients;
    Eigen::MatrixXd PiNabla;
    Eigen::MatrixXd Pi0km1;
    Eigen::MatrixXd Pi0k;

    std::vector<Eigen::MatrixXd> Pi0km1Der;
    Eigen::MatrixXd StabMatrix;
    Eigen::MatrixXd StabMatrixPi0k;
    Eigen::MatrixXd Hmatrix;
    Eigen::LLT<Eigen::MatrixXd> H_km1_LLT;
    Eigen::MatrixXd Cmatrix;
    Eigen::MatrixXd Bmatrix;
    Eigen::MatrixXd Gmatrix;
    Eigen::MatrixXd Dmatrix;
    std::vector<Eigen::MatrixXd> Ematrix;

    Eigen::MatrixXd VanderBoundary;

    Eigen::MatrixXd VanderEdgeDofs;
    Eigen::MatrixXd ScaledHmatrixOnBoundary;

    Eigen::MatrixXd VanderFaceProjections;
    std::vector<Eigen::MatrixXd> FaceScaledMomentsBasis;
    std::vector<Eigen::MatrixXd> FaceProjectedBasisFunctionsValues;
    Eigen::MatrixXd PointEdgeDofsCoordinates;

    std::vector<Eigen::MatrixXd> VanderBoundaryDerivatives;

    Eigen::MatrixXd Qmatrix; // change of basis matrix: pV = mV*Q'
    Eigen::MatrixXd QmatrixInv;
    Eigen::MatrixXd Qmatrixkm1;

    Utilities::VEM_Inertia_Utilities::Inertia_Data inertia_data;
    VEM_PCC_3D_Polyhedron_Geometry inertia_polyhedron;
    double constantStiff;
    double constantMass;

    Eigen::RowVectorXd EdgeInternalPoints;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
