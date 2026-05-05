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

#include "LocalSpace_DF_PCC_2D.hpp"
#include "CommonUtilities.hpp"
#include "VEM_DF_PCC_2D_Creator.hpp"
#include "VTKUtilities.hpp"

namespace Polydim
{
namespace PDETools
{
namespace LocalSpace_DF_PCC_2D
{
//***************************************************************************
ReferenceElement_Data CreateReferenceElement(const MethodTypes &method_type, const unsigned int method_order)
{
    ReferenceElement_Data reference_element_data;
    reference_element_data.Method_Type = method_type;
    reference_element_data.Order = method_order;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {

        switch (method_type)
        {
        case MethodTypes::VEM_DF_PCC_REDUCED:
            reference_element_data.VEM_Type = Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types::VEM_DF_PCC_2D_Reduced_LocalSpace;
            break;
        case MethodTypes::VEM_DF_PCC_FULL:
            reference_element_data.VEM_Type = Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types::VEM_DF_PCC_2D_LocalSpace;
            break;
        default:
            break;
        }

        reference_element_data.VEM_Pressure_ReferenceElement =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_pressure_reference_element(reference_element_data.VEM_Type);
        reference_element_data.VEM_Pressure_ReferenceElement_Data =
            reference_element_data.VEM_Pressure_ReferenceElement->Create(method_order);
        reference_element_data.VEM_Velocity_ReferenceElement =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_velocity_reference_element(reference_element_data.VEM_Type);
        reference_element_data.VEM_Velocity_ReferenceElement_Data =
            reference_element_data.VEM_Velocity_ReferenceElement->Create(method_order);

        reference_element_data.VEM_Pressure_LocalSpace =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_pressure_local_space(reference_element_data.VEM_Type);
        reference_element_data.VEM_Velocity_LocalSpace =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_velocity_local_space(reference_element_data.VEM_Type);
    }
    break;
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return reference_element_data;
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
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
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

        local_space_data.VEM_Pressure_LocalSpace_Data =
            reference_element_data.VEM_Pressure_LocalSpace->CreateLocalSpace(reference_element_data.VEM_Pressure_ReferenceElement_Data,
                                                                             local_space_data.VEM_Geometry);

        local_space_data.VEM_Velocity_LocalSpace_Data =
            reference_element_data.VEM_Velocity_LocalSpace->CreateLocalSpace(reference_element_data.VEM_Velocity_ReferenceElement_Data,
                                                                             local_space_data.VEM_Geometry);
    }
    break;
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return local_space_data;
}
//***************************************************************************
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                                          const LocalSpace_Data &local_space_data,
                                                          const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Velocity_LocalSpace->ComputeBasisFunctionsValues(local_space_data.VEM_Velocity_LocalSpace_Data,
                                                                                           projectionType);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                                          const LocalSpace_Data &local_space_data,
                                                          const Eigen::MatrixXd &points,
                                                          const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Velocity_LocalSpace->ComputeBasisFunctionsValues(
            reference_element_data.VEM_Velocity_ReferenceElement_Data,
            local_space_data.VEM_Geometry,
            local_space_data.VEM_Velocity_LocalSpace_Data,
            projectionType,
            points);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd PressureBasisFunctionsValues(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Pressure_LocalSpace->ComputeBasisFunctionsValues(local_space_data.VEM_Pressure_LocalSpace_Data);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd PressureBasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                             const LocalSpace_Data &local_space_data,
                                             const Eigen::MatrixXd &points)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Pressure_LocalSpace->ComputeBasisFunctionsValues(
            reference_element_data.VEM_Pressure_ReferenceElement_Data,
            local_space_data.VEM_Pressure_LocalSpace_Data,
            points);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd VelocityBasisFunctionsDivergenceValues(const ReferenceElement_Data &reference_element_data,
                                                       const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Velocity_LocalSpace->ComputeBasisFunctionsDivergenceValues(
            local_space_data.VEM_Velocity_LocalSpace_Data);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd VelocityBasisFunctionsValuesOnEdge(const unsigned int &edge_local_index,
                                                   const ReferenceElement_Data &reference_element_data,
                                                   const LocalSpace_Data &local_space_data,
                                                   const Eigen::MatrixXd &pointsCurvilinearCoordinates)
{
    Gedim::Utilities::Unused(edge_local_index);
    Gedim::Utilities::Unused(local_space_data);

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Velocity_LocalSpace->ComputeValuesOnEdge(reference_element_data.VEM_Velocity_ReferenceElement_Data,
                                                                                   pointsCurvilinearCoordinates);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsDerivativeValues(const ReferenceElement_Data &reference_element_data,
                                                                    const LocalSpace_Data &local_space_data,
                                                                    const VEM::DF_PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Velocity_LocalSpace->ComputeBasisFunctionsDerivativeValues(local_space_data.VEM_Velocity_LocalSpace_Data,
                                                                                                     projectionType);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::vector<Eigen::MatrixXd> VelocityBasisFunctionsDerivativeValues(const ReferenceElement_Data &reference_element_data,
                                                                    const LocalSpace_Data &local_space_data,
                                                                    const Eigen::MatrixXd &points,
                                                                    const VEM::DF_PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Velocity_LocalSpace->ComputeBasisFunctionsDerivativeValues(
            reference_element_data.VEM_Velocity_ReferenceElement_Data,
            local_space_data.VEM_Geometry,
            local_space_data.VEM_Velocity_LocalSpace_Data,
            projectionType,
            points);
    }
    case MethodTypes::TAYLOR_HOOD:
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
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return local_space_data.VEM_Velocity_LocalSpace_Data.InternalQuadrature;
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
unsigned int VelocitySize(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return local_space_data.VEM_Velocity_LocalSpace_Data.NumBasisFunctions;
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd VelocityStabilizationMatrix(const ReferenceElement_Data &reference_element_data,
                                            const LocalSpace_Data &local_space_data,
                                            const VEM::DF_PCC::ProjectionTypes &projectionType)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        return reference_element_data.VEM_Velocity_LocalSpace->ComputeDofiDofiStabilizationMatrix(local_space_data.VEM_Velocity_LocalSpace_Data,
                                                                                                  projectionType);
    }
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
Eigen::MatrixXd VelocityEdgeDofsCoordinates(const ReferenceElement_Data &reference_element_data,
                                            const LocalSpace_Data &local_space_data,
                                            const unsigned int edge_local_index)
{
    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        const auto &referenceEdgeDOFsPoint =
            reference_element_data.VEM_Velocity_ReferenceElement_Data.Quadrature.ReferenceEdgeDOFsInternalPoints;
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
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }
}
//***************************************************************************
std::vector<PDETools::DOFs::DOFsManager::MeshDOFsInfo> SetMeshDOFsInfo(
    const ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 2> &boundary_info)
{
    std::vector<PDETools::DOFs::DOFsManager::MeshDOFsInfo> mesh_dofs_info;

    const unsigned int numCell0Ds = mesh.Cell0DTotalNumber();
    const unsigned int numCell1Ds = mesh.Cell1DTotalNumber();
    const unsigned int numCell2Ds = mesh.Cell2DTotalNumber();

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        mesh_dofs_info.resize(4);
    }
    break;
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    const unsigned int num_mesh_do_info = mesh_dofs_info.size();

    for (unsigned int n = 0; n < num_mesh_do_info; n++)
    {
        mesh_dofs_info[n].CellsNumDOFs[0].resize(numCell0Ds);
        mesh_dofs_info[n].CellsBoundaryInfo[0].resize(numCell0Ds);

        mesh_dofs_info[n].CellsNumDOFs[1].resize(numCell1Ds);
        mesh_dofs_info[n].CellsBoundaryInfo[1].resize(numCell1Ds);

        mesh_dofs_info[n].CellsNumDOFs[2].resize(numCell2Ds);
        mesh_dofs_info[n].CellsBoundaryInfo[2].resize(numCell2Ds);
    }

    for (unsigned int c = 0; c < numCell0Ds; ++c)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::VEM_DF_PCC_REDUCED:
        case MethodTypes::VEM_DF_PCC_FULL: {

            for (unsigned int n = 0; n < 2; n++)
            {
                mesh_dofs_info[n].CellsNumDOFs[0][c] = reference_element_data.VEM_Velocity_ReferenceElement_Data.NumDofs0D;
                mesh_dofs_info[n].CellsBoundaryInfo[0][c] = boundary_info[n].at(mesh.Cell0DMarker(c));
            }

            mesh_dofs_info[2].CellsNumDOFs[0][c] = 0;
            mesh_dofs_info[2].CellsBoundaryInfo[0][c] = {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0};

            mesh_dofs_info[3].CellsNumDOFs[0][c] = reference_element_data.VEM_Pressure_ReferenceElement_Data.NumDofs0D;
            mesh_dofs_info[3].CellsBoundaryInfo[0][c] = {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0};
        }
        break;
        case MethodTypes::TAYLOR_HOOD:
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
    }

    for (unsigned int c = 0; c < numCell1Ds; c++)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::VEM_DF_PCC_REDUCED:
        case MethodTypes::VEM_DF_PCC_FULL: {

            for (unsigned int n = 0; n < 2; n++)
            {
                mesh_dofs_info[n].CellsNumDOFs[1][c] = reference_element_data.VEM_Velocity_ReferenceElement_Data.NumDofs1D;
                mesh_dofs_info[n].CellsBoundaryInfo[1][c] = boundary_info[n].at(mesh.Cell1DMarker(c));
            }

            mesh_dofs_info[2].CellsNumDOFs[1][c] = 0;
            mesh_dofs_info[2].CellsBoundaryInfo[1][c] = {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0};

            mesh_dofs_info[3].CellsNumDOFs[1][c] = reference_element_data.VEM_Pressure_ReferenceElement_Data.NumDofs1D;
            mesh_dofs_info[3].CellsBoundaryInfo[1][c] = {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0};
        }
        break;
        case MethodTypes::TAYLOR_HOOD:
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
    }

    for (unsigned int c = 0; c < numCell2Ds; c++)
    {
        switch (reference_element_data.Method_Type)
        {
        case MethodTypes::VEM_DF_PCC_REDUCED:
        case MethodTypes::VEM_DF_PCC_FULL: {

            for (unsigned int n = 0; n < 2; n++)
            {
                mesh_dofs_info[n].CellsNumDOFs[2][c] = 0;
                mesh_dofs_info[n].CellsBoundaryInfo[2][c] = boundary_info[n].at(mesh.Cell2DMarker(c));
            }

            mesh_dofs_info[2].CellsNumDOFs[2][c] = reference_element_data.VEM_Velocity_ReferenceElement_Data.NumDofs2D_BigOPlus +
                                                   reference_element_data.VEM_Velocity_ReferenceElement_Data.NumDofs2D_Divergence;
            mesh_dofs_info[2].CellsBoundaryInfo[2][c] = {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0};

            mesh_dofs_info[3].CellsNumDOFs[2][c] = reference_element_data.VEM_Pressure_ReferenceElement_Data.NumDofs2D;
            mesh_dofs_info[3].CellsBoundaryInfo[2][c] = {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0};
        }
        break;
        case MethodTypes::TAYLOR_HOOD:
        default:
            throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
        }
    }

    return mesh_dofs_info;
}
// ***************************************************************************
Performance_Data ComputePerformance(const ReferenceElement_Data &reference_element_data, const LocalSpace_Data &local_space_data)
{
    Performance_Data performance;

    switch (reference_element_data.Method_Type)
    {
    case MethodTypes::VEM_DF_PCC_REDUCED:
    case MethodTypes::VEM_DF_PCC_FULL: {
        Polydim::VEM::DF_PCC::VEM_DF_PCC_PerformanceAnalysis performanceAnalysis;

        performance.performance_data.vem_analysis_data =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(),
                                        reference_element_data.VEM_Velocity_ReferenceElement_Data.Monomials,
                                        *reference_element_data.VEM_Velocity_LocalSpace,
                                        local_space_data.VEM_Velocity_LocalSpace_Data);

        performance.performance_data.NumInternalQuadraturePoints =
            local_space_data.VEM_Velocity_LocalSpace_Data.InternalQuadrature.Weights.size();
        performance.performance_data.NumBoundaryQuadraturePoints =
            local_space_data.VEM_Velocity_LocalSpace_Data.BoundaryQuadrature.Quadrature.Weights.size();
    }
    break;
    case MethodTypes::TAYLOR_HOOD:
    default:
        throw std::runtime_error("method type " + std::to_string((unsigned int)reference_element_data.Method_Type) + " not supported");
    }

    return performance;
}
// ***************************************************************************
void export_velocity_dofs(const Gedim::GeometryUtilities &geometry_utilities,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Gedim::IArray &right_hand_side,
                          const Gedim::IArray &solution,
                          const Gedim::IArray &solution_strongs,
                          const std::string &file_path)
{

    std::list<Eigen::Vector3d> dofs_coordinate;
    std::list<double> dof_dimension_values;
    std::list<double> dof_cell_index_values;
    std::list<std::array<double, 2>> solution_values;
    std::list<std::array<double, 2>> rhs_values;
    std::list<std::array<double, 2>> dof_global_index_values;
    std::list<std::array<double, 2>> dof_type_values;
    std::list<std::array<double, 2>> dof_boundary_type_values;
    std::list<std::array<double, 2>> dof_boundary_marker_values;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
    {
        for (unsigned int loc_i = 0; loc_i < dofs_data[0].CellsDOFs[0].at(c).size(); ++loc_i)
        {
            dofs_coordinate.push_back(mesh.Cell0DCoordinates(c));
            dof_dimension_values.push_back(0);
            dof_cell_index_values.push_back(c);

            std::array<double, 2> sol;
            std::array<double, 2> rhs;
            std::array<double, 2> dof_global;
            std::array<double, 2> dof_type;
            std::array<double, 2> dof_boundary_type;
            std::array<double, 2> dof_boundary_marker;

            for (unsigned int h = 0; h < 2; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(0).at(c);

                const auto &local_dofs = dofs_data[h].CellsDOFs[0].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h];
                    sol[h] = solution_strongs.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h];
                    sol[h] = solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    rhs[h] = right_hand_side.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const std::vector<double> local_edge_coordinates =
            geometry_utilities.EquispaceCoordinates(dofs_data[0].CellsDOFs[1].at(c).size(), 0.0, 1.0, false);
        const Eigen::Vector3d edge_origin = mesh.Cell1DOriginCoordinates(c);
        const Eigen::Vector3d edge_tangent = mesh.Cell1DEndCoordinates(c) - edge_origin;

        for (unsigned int loc_i = 0; loc_i < dofs_data[0].CellsDOFs[1].at(c).size(); ++loc_i)
        {
            std::array<double, 2> sol;
            std::array<double, 2> rhs;
            std::array<double, 2> dof_global;
            std::array<double, 2> dof_type;
            std::array<double, 2> dof_boundary_type;
            std::array<double, 2> dof_boundary_marker;

            dofs_coordinate.push_back(edge_origin + local_edge_coordinates[loc_i] * edge_tangent);
            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(1);

            for (unsigned int h = 0; h < 2; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(c);

                const auto &local_dofs = dofs_data[h].CellsDOFs[1].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h];
                    sol[h] = solution_strongs.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h];
                    sol[h] = solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    rhs[h] = right_hand_side.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const unsigned int num_loc_dofs = dofs_data[2].CellsDOFs[2].at(c).size();

        if (num_loc_dofs == 0)
            break;

        const auto local_polygon_coordinates = geometry_utilities.EquispaceCoordinates(num_loc_dofs + 1, 0.0, 1.0, true);
        const Eigen::Vector3d polygon_centroid = mesh_geometric_data.Cell2DsCentroids.at(c);
        const auto polygonCentroidEdgesDistance =
            geometry_utilities.PolygonCentroidEdgesDistance(mesh_geometric_data.Cell2DsVertices.at(c),
                                                            mesh_geometric_data.Cell2DsCentroids.at(c),
                                                            mesh_geometric_data.Cell2DsEdgeNormals.at(c));
        const double circle_diameter = 0.5 * geometry_utilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(2);
            dofs_coordinate.push_back(
                polygon_centroid + circle_diameter * Eigen::Vector3d(cos(2.0 * M_PI * local_polygon_coordinates.at(loc_i)),
                                                                     sin(2.0 * M_PI * local_polygon_coordinates.at(loc_i)),
                                                                     0.0));

            std::array<double, 2> sol = {std::nan(""), std::nan("")};
            std::array<double, 2> rhs = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_global = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_type = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_boundary_type = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_boundary_marker = {std::nan(""), std::nan("")};

            for (unsigned int h = 0; h < 1; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h + 2].CellsBoundaryInfo.at(2).at(c);
                const auto &local_dofs = dofs_data[h + 2].CellsDOFs[2].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h + 2];
                    sol[h] = solution_strongs.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h + 2]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h + 2];
                    sol[h] = solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 2]);
                    rhs[h] = right_hand_side.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 2]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    {
        const unsigned int n = dofs_coordinate.size();
        Eigen::MatrixXd coordinates(3, dofs_coordinate.size());
        unsigned int c = 0;
        for (const auto &dof_coordinate : dofs_coordinate)
            coordinates.col(c++) << dof_coordinate;
        const auto dof_dimension_values_data = std::vector<double>(dof_dimension_values.begin(), dof_dimension_values.end());
        const auto dof_cell_index_values_data = std::vector<double>(dof_cell_index_values.begin(), dof_cell_index_values.end());
        Eigen::VectorXd rhs_values_data(2 * n);
        c = 0;
        for (const auto &v : rhs_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                rhs_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd solution_values_data(2 * n);
        c = 0;
        for (const auto &v : solution_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                solution_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_global_index_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_global_index_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_global_index_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_type_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_type_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_type_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_boundary_type_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_boundary_type_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_boundary_type_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_boundary_marker_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_boundary_marker_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_boundary_marker_values_data[2 * c + d] = v[d];
            c++;
        }

        Gedim::VTKUtilities exporter;
        exporter.AddPoints(coordinates,
                           {{"cell_dimension",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_dimension_values_data.size()),
                             dof_dimension_values_data.data()},
                            {"cell_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_cell_index_values_data.size()),
                             dof_cell_index_values_data.data()},
                            {"boundary_type",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_boundary_type_values_data.size()),
                             dof_boundary_type_values_data.data()},
                            {"boundary_marker",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_boundary_marker_values_data.size()),
                             dof_boundary_marker_values_data.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_global_index_values_data.size()),
                             dof_global_index_values_data.data()},
                            {"dof_type",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_type_values_data.size()),
                             dof_type_values_data.data()},
                            {"rhs",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(rhs_values_data.size()),
                             rhs_values_data.data()},
                            {"solution",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(solution_values_data.size()),
                             solution_values_data.data()}});

        exporter.Export(file_path);
    }
}
//***************************************************************************
} // namespace LocalSpace_DF_PCC_2D
} // namespace PDETools
} // namespace Polydim
