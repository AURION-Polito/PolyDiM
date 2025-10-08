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

#ifndef __VEM_PCC_2D_LocalSpace_Data_HPP
#define __VEM_PCC_2D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "Inertia_Utilities.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

class VEM_PCC_2D_Polygon_Geometry final
{
  public:
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

class VEM_PCC_2D_LocalSpace_Data final
{
  public:
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
    Polydim::VEM::Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadrature;

    Gedim::Quadrature::QuadratureData InternalQuadratureKL;
    Polydim::VEM::Quadrature::VEM_Quadrature_2D::Edges_QuadratureData BoundaryQuadratureKL;

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
    Eigen::MatrixXd Cmatrix;
    Eigen::MatrixXd Bmatrix;
    Eigen::MatrixXd Gmatrix;
    Eigen::MatrixXd Dmatrix;
    std::vector<Eigen::MatrixXd> Ematrix;

    Eigen::MatrixXd Qmatrix;
    Eigen::MatrixXd QmatrixInv;
    Eigen::MatrixXd Qmatrixkm1;

    Polydim::Utilities::Inertia_Data inertia_data;
    Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry inertia_polygon;
    double constantStiff;
    double constantMass;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
