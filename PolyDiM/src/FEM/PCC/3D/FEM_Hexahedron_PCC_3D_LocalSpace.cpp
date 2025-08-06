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

#include "FEM_Hexahedron_PCC_3D_LocalSpace.hpp"

using namespace Eigen;

namespace Polydim
{
namespace FEM
{
namespace PCC
{
// ***************************************************************************
FEM_Hexahedron_PCC_3D_LocalSpace_Data FEM_Hexahedron_PCC_3D_LocalSpace::CreateLocalSpace(
    const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data,
    const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    FEM_Hexahedron_PCC_3D_LocalSpace_Data localSpace;

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = polyhedron.Tolerance1D;
    geometry_utilities_config.Tolerance2D = polyhedron.Tolerance2D;
    geometry_utilities_config.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MapHexahedron mapHexahedron(geometry_utilities);
    const std::vector<unsigned int> polyhedronCoordinateSystem =
        geometry_utilities.PolyhedronCoordinateSystem(polyhedron.Vertices, polyhedron.Edges);
    localSpace.MapData = mapHexahedron.Compute(polyhedron.Vertices, polyhedronCoordinateSystem);

    localSpace.Order = reference_element_data.Order;
    localSpace.NumberOfBasisFunctions = reference_element_data.NumBasisFunctions;

    for (unsigned int e = 0; e < 12; ++e)
    {
        const unsigned int edge_origin_index = polyhedron.Edges(0, e);
        const unsigned int edge_end_index = polyhedron.Edges(1, e);

        const auto &reference_edge = reference_element_data.Edges_by_vertices.at({edge_origin_index, edge_end_index});

        localSpace.polyhedron_to_reference_edge_index[e] = reference_edge.first;
        localSpace.polyhedron_to_reference_edge_direction[e] = reference_edge.second;
    }

    for (unsigned int f = 0; f < 6; ++f)
    {
        const unsigned int face_first_edge_index = polyhedron.Faces[f](1, 0);
        const unsigned int face_reference_first_edge_index = localSpace.polyhedron_to_reference_edge_index[face_first_edge_index];
        const unsigned int face_first_second_index = polyhedron.Faces[f](1, 2);
        const unsigned int face_reference_second_edge_index = localSpace.polyhedron_to_reference_edge_index[face_first_second_index];

        localSpace.polyhedron_to_reference_face_index[f] =
            reference_element_data.Faces_by_edges.at({face_reference_first_edge_index, face_reference_second_edge_index});
    }

    localSpace.DofsMeshOrder.resize(localSpace.NumberOfBasisFunctions, 0);

    localSpace.Dof0DsIndex.fill(0);
    for (unsigned int v = 0; v < 8; v++)
    {
        localSpace.Dof0DsIndex[v + 1] = localSpace.Dof0DsIndex[v] + reference_element_data.NumDofs0D;

        unsigned int vertex_dofCounter = reference_element_data.NumDofs0D * v;
        for (unsigned int d = localSpace.Dof0DsIndex[v]; d < localSpace.Dof0DsIndex[v + 1]; d++)
        {
            localSpace.DofsMeshOrder[vertex_dofCounter] = d;
            vertex_dofCounter++;
        }
    }

    localSpace.Dof1DsIndex.fill(localSpace.Dof0DsIndex[8]);
    for (unsigned int e = 0; e < 12; ++e)
    {
        localSpace.Dof1DsIndex[e + 1] = localSpace.Dof1DsIndex[e] + reference_element_data.NumDofs1D;

        const unsigned int ref_e = localSpace.polyhedron_to_reference_edge_index[e];
        const bool ref_e_direction = localSpace.polyhedron_to_reference_edge_direction[e];
        unsigned int edge_dof_counter = reference_element_data.NumDofs0D * 8 + reference_element_data.NumDofs1D * ref_e;

        if (polyhedron.EdgesDirection.at(e) == ref_e_direction)
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e]; d < localSpace.Dof1DsIndex[e + 1]; d++)
            {
                localSpace.DofsMeshOrder[edge_dof_counter] = d;
                edge_dof_counter++;
            }
        }
        else
        {
            for (unsigned int d = localSpace.Dof1DsIndex[e + 1] - 1; d < UINT_MAX && d >= localSpace.Dof1DsIndex[e]; d--)
            {
                localSpace.DofsMeshOrder[edge_dof_counter] = d;
                edge_dof_counter++;
            }
        }
    }

    localSpace.Dof2DsIndex.fill(localSpace.Dof1DsIndex[12]);
    for (unsigned int f = 0; f < 6; ++f)
    {
        localSpace.Dof2DsIndex[f + 1] = localSpace.Dof2DsIndex[f] + reference_element_data.NumDofs2D;

        const unsigned int ref_f = localSpace.polyhedron_to_reference_face_index[f];
        unsigned int face_dof_counter = reference_element_data.NumDofs0D * 8 + reference_element_data.NumDofs1D * 12 +
                                        reference_element_data.NumDofs2D * ref_f;
        if (polyhedron.FacesDirection.at(f))
        {
            for (unsigned int d = localSpace.Dof2DsIndex[f]; d < localSpace.Dof2DsIndex[f + 1]; d++)
            {
                localSpace.DofsMeshOrder[face_dof_counter] = d;
                face_dof_counter++;
            }
        }
        else
        {
            for (unsigned int d = localSpace.Dof2DsIndex[f + 1] - 1; d < UINT_MAX && d >= localSpace.Dof2DsIndex[f]; d--)
            {
                localSpace.DofsMeshOrder[face_dof_counter] = d;
                face_dof_counter++;
            }
        }
    }

    localSpace.Dof3DsIndex.fill(localSpace.Dof2DsIndex[6]);
    localSpace.Dof3DsIndex[1] = localSpace.Dof3DsIndex[0] + reference_element_data.NumDofs3D;

    unsigned int cell_dof_counter =
        reference_element_data.NumDofs0D * 8 + reference_element_data.NumDofs1D * 12 + reference_element_data.NumDofs2D * 6;
    for (unsigned int d = localSpace.Dof3DsIndex[0]; d < localSpace.Dof3DsIndex[1]; d++)
    {
        localSpace.DofsMeshOrder[cell_dof_counter] = d;
        cell_dof_counter++;
    }

    localSpace.Dofs = MapValues(localSpace, Gedim::MapHexahedron::F(localSpace.MapData, reference_element_data.DofPositions));

    FEM_Quadrilateral_PCC_2D_LocalSpace face_local_space;

    for (unsigned int f = 0; f < 6; ++f)
    {
        const auto &face_geometry = polyhedron.Faces_2D_Geometry[f];

        FEM_PCC_2D_Polygon_Geometry fem_face_geometry = {polyhedron.Tolerance1D,
                                                         polyhedron.Tolerance2D,
                                                         face_geometry.Vertices,
                                                         face_geometry.EdgesDirection,
                                                         face_geometry.EdgesTangent,
                                                         face_geometry.EdgesLength};

        localSpace.Boundary_LocalSpace_Data[f] =
            face_local_space.CreateLocalSpace(reference_element_data.BoundaryReferenceElement_Data, fem_face_geometry);
    }

    localSpace.InternalQuadrature = InternalQuadrature(reference_element_data.ReferenceHexahedronQuadrature, localSpace.MapData);
    localSpace.BoundaryQuadrature = BoundaryQuadrature(localSpace.Boundary_LocalSpace_Data, polyhedron);

    return localSpace;
}
// ***************************************************************************
MatrixXd FEM_Hexahedron_PCC_3D_LocalSpace::MapValues(const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
                                                     const Eigen::MatrixXd &referenceValues) const
{
    Eigen::MatrixXd basisFunctionValuesOrdered(referenceValues.rows(), local_space.NumberOfBasisFunctions);

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
        basisFunctionValuesOrdered.col(local_space.DofsMeshOrder.at(d)) << referenceValues.col(d);

    return basisFunctionValuesOrdered;
}
// ***************************************************************************
std::vector<MatrixXd> FEM_Hexahedron_PCC_3D_LocalSpace::MapDerivativeValues(const FEM_Hexahedron_PCC_3D_LocalSpace_Data &local_space,
                                                                            const std::vector<Eigen::MatrixXd> &referenceDerivateValues) const
{
    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(
        3,
        Eigen::MatrixXd::Zero(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int i = 0; i < 3; i++)
    {
        basisFunctionsDerivativeValues[i] = local_space.MapData.QInv(i, i) * referenceDerivateValues[i];
        for (unsigned int j = 0; j < i; j++)
        {
            basisFunctionsDerivativeValues[i] += local_space.MapData.QInv(j, i) * referenceDerivateValues[j];
            basisFunctionsDerivativeValues[j] += local_space.MapData.QInv(i, j) * referenceDerivateValues[i];
        }
    }

    std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValuesOrdered(
        3,
        Eigen::MatrixXd(referenceDerivateValues.at(0).rows(), local_space.NumberOfBasisFunctions));

    for (unsigned int d = 0; d < local_space.NumberOfBasisFunctions; d++)
    {
        const unsigned int &dofOrder = local_space.DofsMeshOrder.at(d);
        basisFunctionsDerivativeValuesOrdered.at(0).col(dofOrder) << basisFunctionsDerivativeValues.at(0).col(d);
        basisFunctionsDerivativeValuesOrdered.at(1).col(dofOrder) << basisFunctionsDerivativeValues.at(1).col(d);
        basisFunctionsDerivativeValuesOrdered.at(2).col(dofOrder) << basisFunctionsDerivativeValues.at(2).col(d);
    }

    return basisFunctionsDerivativeValuesOrdered;
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData FEM_Hexahedron_PCC_3D_LocalSpace::InternalQuadrature(
    const Gedim::Quadrature::QuadratureData &reference_quadrature,
    const Gedim::MapHexahedron::MapHexahedronData &mapData) const
{
    Gedim::Quadrature::QuadratureData quadrature;

    quadrature.Points = Gedim::MapHexahedron::F(mapData, reference_quadrature.Points);
    quadrature.Weights = reference_quadrature.Weights.array() *
                         Gedim::MapHexahedron::DetJ(mapData, reference_quadrature.Points).array().abs();

    return quadrature;
}
// ***************************************************************************
std::array<Gedim::Quadrature::QuadratureData, 6> FEM_Hexahedron_PCC_3D_LocalSpace::BoundaryQuadrature(
    const std::array<FEM_Quadrilateral_PCC_2D_LocalSpace_Data, 6> &faces_local_space_data,
    const FEM_PCC_3D_Polyhedron_Geometry &polyhedron) const
{
    std::array<Gedim::Quadrature::QuadratureData, 6> faces_quadrature;

    for (unsigned int f = 0; f < 6; ++f)
    {
        auto &face_quadrature = faces_quadrature.at(f);
        face_quadrature = faces_local_space_data[f].InternalQuadrature;

        const Eigen::Vector3d &face_translation = polyhedron.FacesTranslation[f];
        const Eigen::Matrix3d &face_rotation = polyhedron.FacesRotationMatrix[f];

        face_quadrature.Points = (face_rotation * face_quadrature.Points).colwise() + face_translation;
    }

    return faces_quadrature;
}
// ***************************************************************************
} // namespace PCC

} // namespace FEM
} // namespace Polydim
