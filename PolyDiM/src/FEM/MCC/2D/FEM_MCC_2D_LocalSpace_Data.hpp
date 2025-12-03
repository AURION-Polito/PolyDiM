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

#ifndef __FEM_MCC_2D_LocalSpace_Data_HPP
#define __FEM_MCC_2D_LocalSpace_Data_HPP

#include "MapTriangle.hpp"
#include "QuadratureData.hpp"

namespace Polydim
{
namespace FEM
{
namespace MCC
{
enum class FEM_MCC_2D_Types
{
    RT_Triangle = 0
};

struct FEM_MCC_2D_Polygon_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;

    Eigen::MatrixXd Vertices;
    Eigen::VectorXd EdgesLength;
    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
    Eigen::MatrixXd EdgesNormal;
};

struct FEM_Triangle_RT_MCC_2D_LocalSpace_Data final
{
    Gedim::MapTriangle::MapTriangleData MapData;

    unsigned int Order;
    unsigned int NumVelocityBasisFunctions;
    unsigned int NumPressureBasisFunctions;

    Eigen::MatrixXd CompatibilityMatrix; // to ensure global comabitibility for edge degrees of freedom

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature;
};

struct FEM_MCC_2D_LocalSpace_Data final
{
    Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data rt_triangle_local_space_data;
    Polydim::FEM::MCC::FEM_MCC_2D_Types fem_type;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature;
    unsigned int NumVelocityBasisFunctions;
    unsigned int NumPressureBasisFunctions;
};

} // namespace MCC
} // namespace FEM
} // namespace Polydim

#endif
