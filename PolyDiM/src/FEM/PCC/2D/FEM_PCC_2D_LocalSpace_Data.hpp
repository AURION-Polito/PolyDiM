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

#ifndef __FEM_PCC_2D_LocalSpace_Data_HPP
#define __FEM_PCC_2D_LocalSpace_Data_HPP

#include "MapParallelogram.hpp"
#include "MapTriangle.hpp"
#include "QuadratureData.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
enum class FEM_PCC_2D_Types
{
    Triangle = 0,
    Quadrilateral = 1
};

enum class QuadrilateralType
{
    Parallelogram = 0,
    Generic = 1
};

struct FEM_PCC_2D_Polygon_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;

    Eigen::MatrixXd Vertices;
    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
    Eigen::VectorXd EdgesLength;
};

struct FEM_Triangle_PCC_2D_LocalSpace_Data final
{
    Gedim::MapTriangle::MapTriangleData MapData;
    Eigen::Matrix3d B_lap;
    unsigned int Order;
    unsigned int NumberOfBasisFunctions;
    unsigned int NumBoundaryBasisFunctions;
    unsigned int NumInternalBasisFunctions;
    Eigen::MatrixXd Dofs;
    std::vector<unsigned int> DofsMeshOrder;
    std::array<unsigned int, 4> Dof0DsIndex;
    std::array<unsigned int, 4> Dof1DsIndex;
    std::array<unsigned int, 2> Dof2DsIndex;
    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature;
};

struct FEM_Quadrilateral_PCC_2D_LocalSpace_Data final
{
    Eigen::MatrixXd Vertices;
    Gedim::MapParallelogram::MapParallelogramData MapData;
    Eigen::Matrix3d B_lap;
    unsigned int Order;
    unsigned int NumberOfBasisFunctions;
    Eigen::MatrixXd Dofs;
    std::vector<unsigned int> DofsMeshOrder;
    std::array<unsigned int, 5> Dof0DsIndex;
    std::array<unsigned int, 5> Dof1DsIndex;
    std::array<unsigned int, 2> Dof2DsIndex;
    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature;
    Polydim::FEM::PCC::QuadrilateralType quadrilateral_type;
};

struct FEM_PCC_2D_LocalSpace_Data final
{
    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data triangle_local_space_data;
    Polydim::FEM::PCC::FEM_Quadrilateral_PCC_2D_LocalSpace_Data quadrilateral_local_space_data;
    Polydim::FEM::PCC::FEM_PCC_2D_Types fem_type;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature;
    unsigned int NumberOfBasisFunctions;
};

} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
