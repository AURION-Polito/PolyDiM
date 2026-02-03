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

#include "LocalSpace_PCC_2D.hpp"
#include "CommonUtilities.hpp"
#include "VTKUtilities.hpp"
#include "ZFEM_PCC_2D_Creator.hpp"
#include <memory>

namespace Polydim
{
namespace PDETools
{
namespace LocalSpace_PCC_2D
{
//***************************************************************************
ReferenceElement_Data CreateReferenceElement(const MethodTypes &method_type, const unsigned int method_order)
{
    ReferenceElement_Data reference_element_data;
    reference_element_data.Method_Type = method_type;
    reference_element_data.Order = method_order;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        reference_element_data.FEM_ReferenceElement = std::make_unique<FEM::PCC::FEM_PCC_2D_ReferenceElement>();
        reference_element_data.FEM_ReferenceElement_Data = reference_element_data.FEM_ReferenceElement->Create(method_order);
        reference_element_data.FEM_LocalSpace = std::make_unique<FEM::PCC::FEM_PCC_2D_LocalSpace>();
    }
    break;
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::VEM_PCC:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_LocalSpace;
            break;
        case MethodTypes::VEM_PCC_Inertia:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Inertia_LocalSpace;
            break;
        case MethodTypes::VEM_PCC_Ortho:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Ortho_LocalSpace;
            break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }

        reference_element_data.VEM_ReferenceElement =
            Polydim::VEM::PCC::create_VEM_PCC_2D_reference_element(reference_element_data.VEM_Type);
        reference_element_data.VEM_ReferenceElement_Data = reference_element_data.VEM_ReferenceElement->Create(method_order);
        reference_element_data.VEM_LocalSpace = Polydim::VEM::PCC::create_VEM_PCC_2D_local_space(reference_element_data.VEM_Type);
    }
    break;
    case MethodTypes::ZFEM_PCC: {
        reference_element_data.ZFEM_ReferenceElement = Polydim::ZFEM::PCC::create_ZFEM_PCC_2D_reference_element(
            ZFEM::PCC::ZFEM_PCC_2D_LocalSpace_Types::ZFEM_PCC_2D_LocalSpace);
        reference_element_data.ZFEM_ReferenceElement_Data = reference_element_data.ZFEM_ReferenceElement->Create(method_order);
        reference_element_data.ZFEM_LocalSpace =
            Polydim::ZFEM::PCC::create_ZFEM_PCC_2D_local_space(ZFEM::PCC::ZFEM_PCC_2D_LocalSpace_Types::ZFEM_PCC_2D_LocalSpace);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return reference_element_data;
}
//***************************************************************************
Gedim::MeshUtilities::MeshGeometricData2DConfig MeshGeometricDataConfigiguration(const ReferenceElement_Data &reference_element_data)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_ReferenceElement_Data.mesh_geometric_data_config;
    }
    break;
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_ReferenceElement_Data.mesh_geometric_data_config;
    }
    break;
    case MethodTypes::ZFEM_PCC: {
        return reference_element_data.ZFEM_ReferenceElement_Data.mesh_geometric_data_config;
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
LocalSpace_Data CreateLocalSpace(const double &geometric_tolerance_1D,
                                 const double &geometric_tolerance_2D,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const ReferenceElement_Data &reference_element_data)
{
    LocalSpace_Data local_space_data;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        local_space_data.FEM_Geometry = {geometric_tolerance_1D,
                                         geometric_tolerance_2D,
                                         mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index)};

        local_space_data.FEM_LocalSpace_Data =
            reference_element_data.FEM_LocalSpace->CreateLocalSpace(reference_element_data.FEM_ReferenceElement_Data,
                                                                    local_space_data.FEM_Geometry);
    }
    break;
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        local_space_data.VEM_Geometry = {geometric_tolerance_1D,
                                         geometric_tolerance_2D,
                                         mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsCentroids.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsAreas.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsDiameters.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsTriangulations.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeNormals.at(cell2D_index)};
        local_space_data.VEM_LocalSpace_Data =
            reference_element_data.VEM_LocalSpace->CreateLocalSpace(reference_element_data.VEM_ReferenceElement_Data,
                                                                    local_space_data.VEM_Geometry);
    }
    break;
    case MethodTypes::ZFEM_PCC: {
        local_space_data.ZFEM_Geometry = {geometric_tolerance_1D,
                                          geometric_tolerance_2D,
                                          mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsAreas.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsDiameters.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsEdgeNormals.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsChebyshevCenter.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsInRadius.at(cell2D_index),
                                          mesh_geometric_data.Cell2DsTriangulationsByChebyshevCenter.at(cell2D_index)};

        local_space_data.ZFEM_LocalSpace_Data =
            reference_element_data.ZFEM_LocalSpace->CreateLocalSpace(reference_element_data.ZFEM_ReferenceElement_Data,
                                                                     local_space_data.ZFEM_Geometry);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return local_space_data;
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                     const LocalSpace_Data &local_space_data,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                  local_space_data.FEM_LocalSpace_Data);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsValues(local_space_data.VEM_LocalSpace_Data, projectionType);
    }
    case MethodTypes::ZFEM_PCC: {
        return reference_element_data.ZFEM_LocalSpace->ComputeBasisFunctionsValues(local_space_data.ZFEM_LocalSpace_Data);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                     const LocalSpace_Data &local_space_data,
                                     const Eigen::MatrixXd &points,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                  local_space_data.FEM_LocalSpace_Data,
                                                                                  points);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.VEM_ReferenceElement_Data,
                                                                                  local_space_data.VEM_LocalSpace_Data,
                                                                                  projectionType,
                                                                                  points);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsLaplacianValues(const ReferenceElement_Data &reference_element_data,
                                              const LocalSpace_Data &local_space_data,
                                              const Polydim::VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsLaplacianValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                           local_space_data.FEM_LocalSpace_Data);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsLaplacianValues(local_space_data.VEM_LocalSpace_Data,
                                                                                           projectionType);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd BasisFunctionsValuesOnEdge(const unsigned int &edge_local_index,
                                           const ReferenceElement_Data &reference_element_data,
                                           const LocalSpace_Data &local_space_data,
                                           const Eigen::MatrixXd &pointsCurvilinearCoordinates)
{
    Gedim::Utilities::Unused(edge_local_index);

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValuesOnEdge(reference_element_data.FEM_ReferenceElement_Data,
                                                                                        local_space_data.FEM_LocalSpace_Data,
                                                                                        pointsCurvilinearCoordinates);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeValuesOnEdge(reference_element_data.VEM_ReferenceElement_Data,
                                                                          pointsCurvilinearCoordinates);
    }
    case MethodTypes::ZFEM_PCC: {
        return reference_element_data.ZFEM_LocalSpace->ComputeValuesOnEdge(reference_element_data.ZFEM_ReferenceElement_Data,
                                                                           local_space_data.ZFEM_LocalSpace_Data,
                                                                           pointsCurvilinearCoordinates);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(const ReferenceElement_Data &reference_element_data,
                                                            const LocalSpace_Data &local_space_data,
                                                            const VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                            local_space_data.FEM_LocalSpace_Data);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(local_space_data.VEM_LocalSpace_Data,
                                                                                            projectionType);
    }
    case MethodTypes::ZFEM_PCC: {
        return reference_element_data.ZFEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(local_space_data.ZFEM_LocalSpace_Data);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(const ReferenceElement_Data &reference_element_data,
                                                            const LocalSpace_Data &local_space_data,
                                                            const Eigen::MatrixXd &points,
                                                            const VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                            local_space_data.FEM_LocalSpace_Data,
                                                                                            points);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(reference_element_data.VEM_ReferenceElement_Data,
                                                                                            local_space_data.VEM_LocalSpace_Data,
                                                                                            projectionType,
                                                                                            points);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Gedim::Quadrature::QuadratureData InternalQuadrature(const ReferenceElement_Data &reference_element_data,
                                                     const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return local_space_data.FEM_LocalSpace_Data.InternalQuadrature;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.InternalQuadrature;
    }
    case MethodTypes::ZFEM_PCC: {
        return local_space_data.ZFEM_LocalSpace_Data.InternalQuadrature;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
unsigned int Size(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.NumBasisFunctions;
    }
    case MethodTypes::ZFEM_PCC: {
        return local_space_data.ZFEM_LocalSpace_Data.NumBasisFunctions;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd StabilizationMatrix(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const VEM::PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return Eigen::MatrixXd::Zero(local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions,
                                     local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeDofiDofiStabilizationMatrix(local_space_data.VEM_LocalSpace_Data,
                                                                                         projectionType);
    }
    case MethodTypes::ZFEM_PCC: {
        return Eigen::MatrixXd::Zero(local_space_data.ZFEM_LocalSpace_Data.NumBasisFunctions,
                                     local_space_data.ZFEM_LocalSpace_Data.NumBasisFunctions);
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd EdgeDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                    const LocalSpace_Data &local_space_data,
                                    const unsigned int edge_local_index)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->EdgeDOFsCoordinates(reference_element_data.FEM_ReferenceElement_Data,
                                                                          local_space_data.FEM_LocalSpace_Data,
                                                                          edge_local_index);
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        const auto &referenceEdgeDOFsPoint = reference_element_data.VEM_ReferenceElement_Data.Quadrature.ReferenceEdgeDOFsInternalPoints;
        const unsigned int num_edge_dofs = referenceEdgeDOFsPoint.cols();

        if (num_edge_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const unsigned int num_edges = local_space_data.VEM_Geometry.Vertices.cols();
        const Eigen::Vector3d edge_origin = local_space_data.VEM_Geometry.EdgesDirection.at(edge_local_index)
                                                ? local_space_data.VEM_Geometry.Vertices.col(edge_local_index)
                                                : local_space_data.VEM_Geometry.Vertices.col((edge_local_index + 1) % num_edges);

        const Eigen::Vector3d edge_tangent = local_space_data.VEM_Geometry.EdgesTangent.col(edge_local_index);
        const double edge_direction = local_space_data.VEM_Geometry.EdgesDirection[edge_local_index] ? 1.0 : -1.0;

        Eigen::MatrixXd edge_dofs_coordinates = Eigen::MatrixXd::Zero(3, num_edge_dofs);
        for (unsigned int r = 0; r < num_edge_dofs; r++)
        {
            edge_dofs_coordinates.col(r) << edge_origin + edge_direction * referenceEdgeDOFsPoint(0, r) * edge_tangent;
        }

        return edge_dofs_coordinates;
    }
    case MethodTypes::ZFEM_PCC: {
        const auto &referenceEdgeDOFsPoint = local_space_data.ZFEM_LocalSpace_Data.ReferenceEdgeDOFsInternalPoints;
        const unsigned int num_edge_dofs = referenceEdgeDOFsPoint.cols();

        if (num_edge_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const unsigned int num_edges = local_space_data.ZFEM_Geometry.Vertices.cols();
        const Eigen::Vector3d edge_origin = local_space_data.ZFEM_Geometry.EdgesDirection.at(edge_local_index)
                                                ? local_space_data.ZFEM_Geometry.Vertices.col(edge_local_index)
                                                : local_space_data.ZFEM_Geometry.Vertices.col((edge_local_index + 1) % num_edges);

        const Eigen::Vector3d edge_tangent = local_space_data.ZFEM_Geometry.EdgesTangent.col(edge_local_index);
        const double edge_direction = local_space_data.ZFEM_Geometry.EdgesDirection[edge_local_index] ? 1.0 : -1.0;

        Eigen::MatrixXd edge_dofs_coordinates = Eigen::MatrixXd::Zero(3, num_edge_dofs);
        for (unsigned int r = 0; r < num_edge_dofs; r++)
        {
            edge_dofs_coordinates.col(r) << edge_origin + edge_direction * referenceEdgeDOFsPoint(0, r) * edge_tangent;
        }

        return edge_dofs_coordinates;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
// ***************************************************************************
Gedim::Quadrature::QuadratureData InternalDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                                          const LocalSpace_Data &local_space_data)
{
    Gedim::Quadrature::QuadratureData face_dofs_coordinates;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {

        face_dofs_coordinates.Points =
            reference_element_data.FEM_LocalSpace->InternalDOFsCoordinates(reference_element_data.FEM_ReferenceElement_Data,
                                                                           local_space_data.FEM_LocalSpace_Data);
        return face_dofs_coordinates;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {

        face_dofs_coordinates.Points = local_space_data.VEM_LocalSpace_Data.InternalQuadrature.Points;

        face_dofs_coordinates.Weights = local_space_data.VEM_LocalSpace_Data.InternalQuadrature.Weights;

        return face_dofs_coordinates;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
// ***************************************************************************
Eigen::VectorXd InternalDofs(const ReferenceElement_Data &reference_element_data,
                             const LocalSpace_Data &local_space_data,
                             const Eigen::VectorXd &values_at_dofs,
                             const Gedim::Quadrature::QuadratureData &internal_dofs_coordinates)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        return values_at_dofs;
    }
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {

        const auto scaled_polynomial =
            reference_element_data.VEM_LocalSpace->ComputeScaledPolynomialsValues(local_space_data.VEM_LocalSpace_Data);

        return scaled_polynomial.transpose() * internal_dofs_coordinates.Weights.asDiagonal() * values_at_dofs;
    }
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
PDETools::DOFs::DOFsManager::MeshDOFsInfo SetMeshDOFsInfo(
    const ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> &boundary_info)
{
    PDETools::DOFs::DOFsManager::MeshDOFsInfo mesh_dof_info;

    const unsigned int numCell0Ds = mesh.Cell0DTotalNumber();

    mesh_dof_info.CellsNumDOFs[0].resize(numCell0Ds);
    mesh_dof_info.CellsBoundaryInfo[0].resize(numCell0Ds);

    for (unsigned int c = 0; c < numCell0Ds; ++c)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::FEM_PCC: {
            mesh_dof_info.CellsNumDOFs[0][c] = reference_element_data.FEM_ReferenceElement_Data.NumDofs0D;
        }
        break;
        case MethodTypes::VEM_PCC:
        case MethodTypes::VEM_PCC_Inertia:
        case MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[0][c] = reference_element_data.VEM_ReferenceElement_Data.NumDofs0D;
        }
        break;
        case MethodTypes::ZFEM_PCC: {
            mesh_dof_info.CellsNumDOFs[0][c] = reference_element_data.ZFEM_ReferenceElement_Data.NumDofs0D;
        }
        break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }

        mesh_dof_info.CellsBoundaryInfo[0][c] = boundary_info.at(mesh.Cell0DMarker(c));
    }

    const unsigned int numCell1Ds = mesh.Cell1DTotalNumber();

    mesh_dof_info.CellsNumDOFs[1].resize(numCell1Ds);
    mesh_dof_info.CellsBoundaryInfo[1].resize(numCell1Ds);

    for (unsigned int c = 0; c < numCell1Ds; c++)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::FEM_PCC: {
            mesh_dof_info.CellsNumDOFs[1][c] = reference_element_data.FEM_ReferenceElement_Data.NumDofs1D;
        }
        break;
        case MethodTypes::VEM_PCC:
        case MethodTypes::VEM_PCC_Inertia:
        case MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[1][c] = reference_element_data.VEM_ReferenceElement_Data.NumDofs1D;
        }
        break;
        case MethodTypes::ZFEM_PCC: {
            mesh_dof_info.CellsNumDOFs[1][c] = reference_element_data.ZFEM_ReferenceElement_Data.NumDofs1D;
        }
        break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
        mesh_dof_info.CellsBoundaryInfo[1][c] = boundary_info.at(mesh.Cell1DMarker(c));
    }

    const unsigned int numCell2Ds = mesh.Cell2DTotalNumber();

    mesh_dof_info.CellsNumDOFs[2].resize(numCell2Ds);
    mesh_dof_info.CellsBoundaryInfo[2].resize(numCell2Ds);

    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::FEM_PCC: {
            if (mesh.Cell2DNumberVertices(c) == 3)
                mesh_dof_info.CellsNumDOFs[2][c] =
                    reference_element_data.FEM_ReferenceElement_Data.triangle_reference_element_data.NumDofs2D;
            else if (mesh.Cell2DNumberVertices(c) == 4)
                mesh_dof_info.CellsNumDOFs[2][c] =
                    reference_element_data.FEM_ReferenceElement_Data.quadrilateral_reference_element_data.NumDofs2D;
            else
                throw std::runtime_error("not valid element");
        }
        break;
        case MethodTypes::VEM_PCC:
        case MethodTypes::VEM_PCC_Inertia:
        case MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[2][c] = reference_element_data.VEM_ReferenceElement_Data.NumDofs2D;
        }
        break;
        case MethodTypes::ZFEM_PCC: {
            mesh_dof_info.CellsNumDOFs[2][c] = reference_element_data.ZFEM_ReferenceElement_Data.NumDofs2D;
        }
        break;
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
        mesh_dof_info.CellsBoundaryInfo[2][c] = boundary_info.at(mesh.Cell2DMarker(c));
    }

    return mesh_dof_info;
}
//***************************************************************************
Performance_Data ComputePerformance(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    Performance_Data performance;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::FEM_PCC: {
        performance.performance_data.NumInternalQuadraturePoints =
            local_space_data.FEM_LocalSpace_Data.InternalQuadrature.Weights.size();
    }
    break;
    case MethodTypes::VEM_PCC:
    case MethodTypes::VEM_PCC_Inertia:
    case MethodTypes::VEM_PCC_Ortho: {
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        performance.performance_data.vem_analysis_data =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(),
                                        reference_element_data.VEM_ReferenceElement_Data.Monomials,
                                        *reference_element_data.VEM_LocalSpace,
                                        local_space_data.VEM_LocalSpace_Data);

        performance.performance_data.NumInternalQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data.InternalQuadrature.Weights.size();
        performance.performance_data.NumBoundaryQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data.BoundaryQuadrature.Quadrature.Weights.size();
    }
    break;
    case MethodTypes::ZFEM_PCC: {

        Polydim::ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis performanceAnalysis;

        performance.performance_data.NumInternalQuadraturePoints =
            local_space_data.ZFEM_LocalSpace_Data.InternalQuadrature.Weights.size();

        performance.performance_data.zfem_analysis_data =
            performanceAnalysis.Compute2D(Polydim::Utilities::Monomials_2D(),
                                          reference_element_data.ZFEM_ReferenceElement_Data.monomials_data,
                                          *reference_element_data.ZFEM_LocalSpace,
                                          local_space_data.ZFEM_LocalSpace_Data);
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return performance;
}
//***************************************************************************
void export_dofs(const Gedim::GeometryUtilities &geometry_utilities,
                 const Gedim::MeshMatricesDAO &mesh,
                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                 const Gedim::IArray &right_hand_side,
                 const Gedim::IArray &solution,
                 const Gedim::IArray &solution_strongs,
                 const std::string &file_path)
{
    std::list<Eigen::Vector3d> dofs_coordinate;
    std::list<double> solution_values;
    std::list<double> rhs_values;
    std::list<double> dof_global_index_values;
    std::list<double> dof_type_values;
    std::list<double> dof_cell_index_values;
    std::list<double> dof_dimension_values;
    std::list<double> dof_boundary_type_values;
    std::list<double> dof_boundary_marker_values;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[0].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(0);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            dofs_coordinate.push_back(mesh.Cell0DCoordinates(c));
            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(solution_strongs.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(right_hand_side.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[1].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        const std::vector<double> local_edge_coordinates = geometry_utilities.EquispaceCoordinates(num_loc_dofs, 0.0, 1.0, false);
        const Eigen::Vector3d edge_origin = mesh.Cell1DOriginCoordinates(c);
        const Eigen::Vector3d edge_tangent = mesh.Cell1DEndCoordinates(c) - edge_origin;

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(1);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            dofs_coordinate.push_back(edge_origin + local_edge_coordinates[loc_i] * edge_tangent);
            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(solution_strongs.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(right_hand_side.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[2].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        const auto local_polygon_coordinates = geometry_utilities.EquispaceCoordinates(num_loc_dofs + 1, 0.0, 1.0, true);
        const Eigen::Vector3d polygon_centroid = mesh_geometric_data.Cell2DsCentroids.at(c);
        const auto polygonCentroidEdgesDistance =
            geometry_utilities.PolygonCentroidEdgesDistance(mesh_geometric_data.Cell2DsVertices.at(c),
                                                            mesh_geometric_data.Cell2DsCentroids.at(c),
                                                            mesh_geometric_data.Cell2DsEdgeNormals.at(c));
        const double circle_diameter = 0.5 * geometry_utilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(2);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);

            dofs_coordinate.push_back(
                polygon_centroid +
                circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  0.0));

            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(solution_strongs.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(right_hand_side.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    Eigen::MatrixXd coordinates(3, dofs_coordinate.size());
    unsigned int c = 0;
    for (const auto &dof_coordinate : dofs_coordinate)
        coordinates.col(c++) << dof_coordinate;
    const auto rhs_values_data = std::vector<double>(rhs_values.begin(), rhs_values.end());
    const auto solution_values_data = std::vector<double>(solution_values.begin(), solution_values.end());
    const auto dof_global_index_values_data =
        std::vector<double>(dof_global_index_values.begin(), dof_global_index_values.end());
    const auto dof_type_values_data = std::vector<double>(dof_type_values.begin(), dof_type_values.end());
    const auto dof_cell_index_values_data = std::vector<double>(dof_cell_index_values.begin(), dof_cell_index_values.end());
    const auto dof_dimension_values_data = std::vector<double>(dof_dimension_values.begin(), dof_dimension_values.end());
    const auto dof_boundary_type_values_data =
        std::vector<double>(dof_boundary_type_values.begin(), dof_boundary_type_values.end());
    const auto dof_boundary_marker_values_data =
        std::vector<double>(dof_boundary_marker_values.begin(), dof_boundary_marker_values.end());

    Gedim::VTKUtilities exporter;
    exporter.AddPoints(
        coordinates,
        {{"cell_dimension",
          Gedim::VTPProperty::Formats::Points,
          static_cast<unsigned int>(dof_dimension_values_data.size()),
          dof_dimension_values_data.data()},
         {"cell_index",
          Gedim::VTPProperty::Formats::Points,
          static_cast<unsigned int>(dof_cell_index_values_data.size()),
          dof_cell_index_values_data.data()},
         {"boundary_type",
          Gedim::VTPProperty::Formats::Points,
          static_cast<unsigned int>(dof_boundary_type_values_data.size()),
          dof_boundary_type_values_data.data()},
         {"boundary_marker",
          Gedim::VTPProperty::Formats::Points,
          static_cast<unsigned int>(dof_boundary_marker_values_data.size()),
          dof_boundary_marker_values_data.data()},
         {"dof_global_index",
          Gedim::VTPProperty::Formats::Points,
          static_cast<unsigned int>(dof_global_index_values_data.size()),
          dof_global_index_values_data.data()},
         {"dof_type",
          Gedim::VTPProperty::Formats::Points,
          static_cast<unsigned int>(dof_type_values_data.size()),
          dof_type_values_data.data()},
         {"rhs", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(rhs_values_data.size()), rhs_values_data.data()},
         {"solution",
          Gedim::VTPProperty::Formats::Points,
          static_cast<unsigned int>(solution_values_data.size()),
          solution_values_data.data()}});

    exporter.Export(file_path);
}
//***************************************************************************
} // namespace LocalSpace_PCC_2D
} // namespace PDETools
} // namespace Polydim
