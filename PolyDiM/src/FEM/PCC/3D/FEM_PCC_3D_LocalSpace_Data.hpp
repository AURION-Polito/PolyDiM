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

#ifndef __FEM_Hexahedron_PCC_3D_LocalSpace_Data_HPP
#define __FEM_Hexahedron_PCC_3D_LocalSpace_Data_HPP

#include "FEM_PCC_2D_LocalSpace_Data.hpp"
#include "MapHexahedron.hpp"
#include "MapTetrahedron.hpp"
#include "QuadratureData.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

enum class FEM_PCC_3D_Types
{
    Tetrahedron = 0,
    Hexahedron = 1
};

struct FEM_PCC_3D_Polyhedron_Geometry final
{
    double Tolerance1D;
    double Tolerance2D;
    double Tolerance3D;

    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    std::vector<Eigen::MatrixXi> Faces;
    std::vector<FEM_PCC_2D_Polygon_Geometry> Faces_2D_Geometry;
    std::vector<bool> EdgesDirection;
    std::vector<bool> FacesDirection;
    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
};

struct FEM_Hexahedron_PCC_3D_LocalSpace_Data final
{
    Gedim::MapHexahedron::MapHexahedronData MapData;
    unsigned int Order;
    unsigned int NumberOfBasisFunctions;
    Eigen::MatrixXd Dofs;
    std::array<unsigned int, 12> polyhedron_to_reference_edge_index;
    std::array<bool, 12> polyhedron_to_reference_edge_direction;
    std::array<unsigned int, 6> polyhedron_to_reference_face_index;
    std::vector<unsigned int> DofsMeshOrder;
    std::array<unsigned int, 9> Dof0DsIndex;
    std::array<unsigned int, 13> Dof1DsIndex;
    std::array<unsigned int, 7> Dof2DsIndex;
    std::array<unsigned int, 2> Dof3DsIndex;
    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::array<FEM_Quadrilateral_PCC_2D_LocalSpace_Data, 6> Boundary_LocalSpace_Data;
    std::array<Gedim::Quadrature::QuadratureData, 6> BoundaryQuadrature;
};

struct FEM_Tetrahedron_PCC_3D_LocalSpace_Data final
{
    Gedim::MapTetrahedron::MapTetrahedronData MapData;
    unsigned int Order;
    unsigned int NumberOfBasisFunctions;
    Eigen::MatrixXd Dofs;
    std::array<unsigned int, 6> polyhedron_to_reference_edge_index;
    std::array<bool, 6> polyhedron_to_reference_edge_direction;
    std::array<unsigned int, 4> polyhedron_to_reference_face_index;
    std::vector<unsigned int> DofsMeshOrder;
    std::array<unsigned int, 5> Dof0DsIndex;
    std::array<unsigned int, 7> Dof1DsIndex;
    std::array<unsigned int, 5> Dof2DsIndex;
    std::array<unsigned int, 2> Dof3DsIndex;
    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::array<FEM_Triangle_PCC_2D_LocalSpace_Data, 4> Boundary_LocalSpace_Data;
    std::array<Gedim::Quadrature::QuadratureData, 4> BoundaryQuadrature;
};

struct FEM_PCC_3D_LocalSpace_Data final
{
    FEM_Hexahedron_PCC_3D_LocalSpace_Data hexahedron_local_space_data;
    FEM_Tetrahedron_PCC_3D_LocalSpace_Data tetrahedron_local_space_data;
    FEM_PCC_3D_Types fem_type;

    Gedim::Quadrature::QuadratureData InternalQuadrature;
    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature;
    unsigned int NumberOfBasisFunctions;
};

} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
