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

#include "local_space.hpp"
#include <memory>

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_2D
{
namespace local_space
{
//***************************************************************************
ReferenceElement_Data CreateReferenceElement(const Program_configuration::MethodTypes &method_type, const unsigned int method_order)
{
    ReferenceElement_Data reference_element_data;
    reference_element_data.Method_Type = method_type;
    reference_element_data.Order = method_order;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_PCC: {
        reference_element_data.FEM_ReferenceElement = std::make_unique<FEM::PCC::FEM_PCC_2D_ReferenceElement>();
        reference_element_data.FEM_ReferenceElement_Data = reference_element_data.FEM_ReferenceElement->Create(method_order);
        reference_element_data.FEM_LocalSpace = std::make_unique<FEM::PCC::FEM_PCC_2D_LocalSpace>();
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        switch (reference_element_data.Method_Type)
        {
        case Program_configuration::MethodTypes::VEM_PCC:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_PCC_Inertia:
            reference_element_data.VEM_Type = VEM::PCC::VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Inertia_LocalSpace;
            break;
        case Program_configuration::MethodTypes::VEM_PCC_Ortho:
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
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return reference_element_data;
}
//***************************************************************************
LocalSpace_Data CreateLocalSpace(const Polydim::examples::Elliptic_PCC_2D::Program_configuration &config,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const ReferenceElement_Data &reference_element_data)
{
    LocalSpace_Data local_space_data;

    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_PCC: {
        local_space_data.FEM_Geometry = {config.GeometricTolerance1D(),
                                         config.GeometricTolerance2D(),
                                         mesh_geometric_data.Cell2DsVertices.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index),
                                         mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index)};

        local_space_data.FEM_LocalSpace_Data =
            reference_element_data.FEM_LocalSpace->CreateLocalSpace(reference_element_data.FEM_ReferenceElement_Data,
                                                                    local_space_data.FEM_Geometry);
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        local_space_data.VEM_Geometry = {config.GeometricTolerance1D(),
                                         config.GeometricTolerance2D(),
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                  local_space_data.FEM_LocalSpace_Data);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsValues(local_space_data.VEM_LocalSpace_Data, projectionType);
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsLaplacianValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                           local_space_data.FEM_LocalSpace_Data);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
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
    switch (reference_element_data.Method_Type)
    {
    case Program_configuration::MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsValuesOnEdge(reference_element_data.FEM_ReferenceElement_Data,
                                                                                        local_space_data.FEM_LocalSpace_Data,
                                                                                        pointsCurvilinearCoordinates);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeValuesOnEdge(reference_element_data.VEM_ReferenceElement_Data,
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(reference_element_data.FEM_ReferenceElement_Data,
                                                                                            local_space_data.FEM_LocalSpace_Data);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeBasisFunctionsDerivativeValues(local_space_data.VEM_LocalSpace_Data,
                                                                                            projectionType);
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        return local_space_data.FEM_LocalSpace_Data.InternalQuadrature;
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.InternalQuadrature;
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        return local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions;
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return local_space_data.VEM_LocalSpace_Data.NumBasisFunctions;
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        return Eigen::MatrixXd::Zero(local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions,
                                     local_space_data.FEM_LocalSpace_Data.NumberOfBasisFunctions);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        return reference_element_data.VEM_LocalSpace->ComputeDofiDofiStabilizationMatrix(local_space_data.VEM_LocalSpace_Data,
                                                                                         projectionType);
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        return reference_element_data.FEM_LocalSpace->EdgeDOFsCoordinates(reference_element_data.FEM_ReferenceElement_Data,
                                                                          local_space_data.FEM_LocalSpace_Data,
                                                                          edge_local_index);
    }
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
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
        case Program_configuration::MethodTypes::FEM_PCC: {
            mesh_dof_info.CellsNumDOFs[0][c] = reference_element_data.FEM_ReferenceElement_Data.NumDofs0D;
        }
        break;
        case Program_configuration::MethodTypes::VEM_PCC:
        case Program_configuration::MethodTypes::VEM_PCC_Inertia:
        case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[0][c] = reference_element_data.VEM_ReferenceElement_Data.NumDofs0D;
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
        case Program_configuration::MethodTypes::FEM_PCC: {
            mesh_dof_info.CellsNumDOFs[1][c] = reference_element_data.FEM_ReferenceElement_Data.NumDofs1D;
        }
        break;
        case Program_configuration::MethodTypes::VEM_PCC:
        case Program_configuration::MethodTypes::VEM_PCC_Inertia:
        case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[1][c] = reference_element_data.VEM_ReferenceElement_Data.NumDofs1D;
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
        case Program_configuration::MethodTypes::FEM_PCC: {
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
        case Program_configuration::MethodTypes::VEM_PCC:
        case Program_configuration::MethodTypes::VEM_PCC_Inertia:
        case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
            mesh_dof_info.CellsNumDOFs[2][c] = reference_element_data.VEM_ReferenceElement_Data.NumDofs2D;
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
    case Program_configuration::MethodTypes::FEM_PCC: {
        performance.VEM_Performance_Data.NumInternalQuadraturePoints =
            local_space_data.FEM_LocalSpace_Data.InternalQuadrature.Weights.size();
    }
    break;
    case Program_configuration::MethodTypes::VEM_PCC:
    case Program_configuration::MethodTypes::VEM_PCC_Inertia:
    case Program_configuration::MethodTypes::VEM_PCC_Ortho: {
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        performance.VEM_Performance_Data.Analysis =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(),
                                        reference_element_data.VEM_ReferenceElement_Data.Monomials,
                                        *reference_element_data.VEM_LocalSpace,
                                        local_space_data.VEM_LocalSpace_Data);

        performance.VEM_Performance_Data.NumInternalQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data.InternalQuadrature.Weights.size();
        performance.VEM_Performance_Data.NumBoundaryQuadraturePoints =
            local_space_data.VEM_LocalSpace_Data.BoundaryQuadrature.Quadrature.Weights.size();
    }
    break;
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return performance;
}
//***************************************************************************
} // namespace local_space
} // namespace Elliptic_PCC_2D
} // namespace examples
} // namespace Polydim
