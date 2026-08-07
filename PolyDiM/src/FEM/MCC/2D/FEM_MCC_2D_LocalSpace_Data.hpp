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

/// @brief Enumeration of the available 2D mixed (velocity/pressure) finite element types.
///
/// Each value identifies a concrete local space implementation that can be built on a
/// two-dimensional element for the mixed MCC formulation.
enum class FEM_MCC_2D_Types
{
    RT_Triangle = 0 ///< Raviart-Thomas element defined on a triangle.
};

/// @brief Geometric description of a polygon on which a 2D MCC local space is constructed.
///
/// Collects the vertices together with the per-edge metric quantities (lengths, tangents,
/// normals and orientations) required to assemble the local mixed finite element space and
/// to enforce the H(div) degrees of freedom on the polygon boundary.
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

/// @brief Local space data for a Raviart-Thomas (RT) triangular element in the 2D MCC formulation.
///
/// Holds everything needed to evaluate the RT velocity/pressure basis on a single triangle:
/// the reference-to-physical mapping, the polynomial order, the basis function counts and the
/// quadrature rules used on the element interior and on its boundary.
struct FEM_Triangle_RT_MCC_2D_LocalSpace_Data final
{
    Gedim::MapTriangle::MapTriangleData MapData;

    unsigned int Order;
    unsigned int NumVelocityBasisFunctions; ///< Number of velocity (H(div)) basis functions on the element.
    unsigned int NumPressureBasisFunctions; ///< Number of pressure (L2) basis functions on the element.
    std::array<bool, 3> EdgesDirection;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature;
};

/// @brief Aggregated local space data for a 2D MCC element.
///
/// Wraps the concrete element-specific data (currently the RT triangle local space) together
/// with the selected @ref FEM_MCC_2D_Types and the quadrature rules and basis function counts
/// exposed to the assembler, providing a uniform interface regardless of the underlying element type.
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
