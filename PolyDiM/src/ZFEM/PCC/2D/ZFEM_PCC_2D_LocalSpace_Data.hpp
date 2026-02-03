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

#ifndef __ZFEM_PCC_2D_LocalSpace_Data_HPP
#define __ZFEM_PCC_2D_LocalSpace_Data_HPP

#include "Eigen/Eigen"
#include "FEM_PCC_2D_LocalSpace_Data.hpp"

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{

struct ZFEM_PCC_2D_Polygon_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;

    Eigen::MatrixXd Vertices;
    double Measure;
    double Diameter;
    Eigen::VectorXd EdgesLength;
    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
    Eigen::MatrixXd EdgesNormal;
    Eigen::MatrixXd ChebishevCenter;
    double InRadius;
    std::vector<Eigen::Matrix3d> Traingulations;
};

struct ZFEM_PCC_2D_LocalSpace_Data final
{
    unsigned int Dimension;
    unsigned int Order;

    Eigen::MatrixXi local_to_total;

    Gedim::Quadrature::QuadratureData InternalQuadrature;

    std::vector<Polydim::FEM::PCC::FEM_PCC_2D_Polygon_Geometry> fem_geometry;
    std::vector<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data> fem_local_space_data;
    std::vector<double> triangles_area;

    Eigen::MatrixXd fem_basis_functions_values;
    std::vector<Eigen::MatrixXd> fem_basis_functions_derivative_values;

    Eigen::MatrixXd ZFEM_basis_functions_values;
    std::vector<Eigen::MatrixXd> ZFEM_basis_functions_derivative_values;

    unsigned int NumVertexBasisFunctions;
    unsigned int NumEdgeBasisFunctions;
    unsigned int NumBoundaryBasisFunctions;
    unsigned int NumInternalBasisFunctions;
    unsigned int NumBasisFunctions;

    unsigned int NumVirtualBasisFunctions;
    unsigned int NumTotalBasisFunctions;

    unsigned int Nk;
    Eigen::MatrixXd VanderInternal;
    Eigen::MatrixXd VanderDOFs;
    Eigen::MatrixXd VanderVirtuals;
    std::vector<Eigen::MatrixXd> VanderInternalDerivatives;
    Eigen::MatrixXd Dmatrix;

    Eigen::RowVectorXd ReferenceEdgeDOFsInternalPoints;
    Eigen::VectorXd EdgeBasisCoefficients;

    Eigen::MatrixXd VirtualNodes;
    Eigen::MatrixXd VirtualWeights;
    Eigen::MatrixXd DOFsCoordinates;
    Eigen::MatrixXd FinerNodes;

    Eigen::MatrixXd KernelVertices;
    Eigen::MatrixXd KernelEdgesNormal;
    Eigen::Vector3d KernelIncenter;
    double KernelInRadius;
};

} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
