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

#include "assembler.hpp"

#include "EllipticEquation.hpp"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_2D
{
// ***************************************************************************
void Assembler::ComputeStrongTerm(const unsigned int &cell2DIndex,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                  const local_space::ReferenceElement_Data &reference_element_data,
                                  const local_space::LocalSpace_Data &local_space_data,
                                  const test::I_Test &test,
                                  Elliptic_MCC_2D_Problem_Data &assembler_data) const
{
    for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(cell2DIndex); e++)
    {
        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex, e);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);
        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
            local_dofs.size() == 0)
            continue;

        const auto edge_dofs_coordinates = local_space::EdgeDofsCoordinates(reference_element_data, local_space_data, e);

        const VectorXd neumann_values = test.strong_boundary_condition(boundary_info.Marker, edge_dofs_coordinates.Points);
        const VectorXd strong_boundary_values =
            local_space::EdgeDofs(reference_element_data, local_space_data, e, edge_dofs_coordinates, neumann_values);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                assembler_data.solutionNeumann.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }
}
// ***************************************************************************
void Assembler::ComputeWeakTerm(const unsigned int cell2DIndex,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                const local_space::ReferenceElement_Data &reference_element_data,
                                const local_space::LocalSpace_Data &local_space_data,
                                const test::I_Test &test,
                                Elliptic_MCC_2D_Problem_Data &assembler_data) const
{
    for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(cell2DIndex); e++)
    {
        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex, e);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
            continue;

        // compute values of Neumann condition
        const auto quadrature = local_space::EdgeQuadrature(reference_element_data, local_space_data, e);

        const auto dirichlet_values = test.weak_boundary_condition(boundary_info.Marker, quadrature.Points);
        const auto basis_function_values =
            local_space::VelocityBasisFunctionsValuesOnEdges(e, reference_element_data, local_space_data, quadrature.Points);
        const VectorXd weak_boundary_values = -basis_function_values.transpose() * quadrature.Weights.asDiagonal() * dirichlet_values;

        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); loc_i++)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                continue;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index, weak_boundary_values(loc_i));
            }
            break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }
}
// ***************************************************************************
Assembler::Elliptic_MCC_2D_Problem_Data Assembler::Assemble(
    const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const local_space::ReferenceElement_Data &reference_element_data,
    const Polydim::examples::Elliptic_MCC_2D::test::I_Test &test) const
{
    Elliptic_MCC_2D_Problem_Data result;
    result.globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    result.neumannMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    result.rightHandSide.SetSize(count_dofs.num_total_dofs);
    result.solution.SetSize(count_dofs.num_total_dofs);
    result.solutionNeumann.SetSize(count_dofs.num_total_strong);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {

        const auto local_space_data = local_space::CreateLocalSpace(config, mesh_geometric_data, c, reference_element_data);

        const auto velocity_basis_functions_values =
            local_space::VelocityBasisFunctionsValues(reference_element_data, local_space_data);
        const auto velocity_basis_functions_divergence_values =
            local_space::VelocityBasisFunctionsDivergenceValues(reference_element_data, local_space_data);
        const auto pressure_basis_functions_values =
            local_space::PressureBasisFunctionsValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature = local_space::InternalQuadrature(reference_element_data, local_space_data);

        const auto reaction_term_values = test.reaction_term(cell2D_internal_quadrature.Points);
        const auto advection_term_values = test.mixed_advection_term(cell2D_internal_quadrature.Points);
        const auto diffusion_term_values = test.inverse_diffusion_term(cell2D_internal_quadrature.Points);
        const auto source_term_values = test.source_term(cell2D_internal_quadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                           velocity_basis_functions_values,
                                                           cell2D_internal_quadrature.Weights);

        double kmax = 0.0;
        for (const auto &diffusion_term : diffusion_term_values)
        {
            const double max_k = diffusion_term.cwiseAbs().maxCoeff();
            kmax = kmax < max_k ? max_k : kmax;
        }

        local_A += kmax * local_space::StabilizationMatrix(reference_element_data, local_space_data);

        const auto local_M = equation.ComputeCellReactionMatrix(reaction_term_values,
                                                                pressure_basis_functions_values,
                                                                cell2D_internal_quadrature.Weights);

        const auto local_T = equation.ComputeCellAdvectionMatrix(advection_term_values,
                                                                 pressure_basis_functions_values,
                                                                 velocity_basis_functions_values,
                                                                 cell2D_internal_quadrature.Weights);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        cell2D_internal_quadrature.Weights.asDiagonal() * velocity_basis_functions_divergence_values;

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               pressure_basis_functions_values,
                                                               cell2D_internal_quadrature.Weights);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(c).size();

        Eigen::MatrixXd elemental_matrix = MatrixXd::Zero(local_count_dofs.num_total_dofs, local_count_dofs.num_total_dofs);
        Eigen::VectorXd elemental_rhs = VectorXd::Zero(local_count_dofs.num_total_dofs);
        elemental_matrix << local_A, -(local_B + local_T).transpose(), local_B, local_M;
        elemental_rhs << VectorXd::Zero(local_count_dofs.num_total_dofs - num_local_dofs_pressure), local_rhs;

        assert(local_space::VelocitySize(reference_element_data, local_space_data) ==
               local_count_dofs.num_total_dofs - num_local_dofs_pressure);

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data = {
            {std::cref(dofs_data[0]), std::cref(dofs_data[1])},
            local_count_dofs.offsets_DOFs,
            count_dofs.offsets_DOFs,
            count_dofs.offsets_Strongs};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          elemental_matrix,
                                                                                          elemental_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.neumannMatrixA,
                                                                                          result.rightHandSide);

        ComputeWeakTerm(c, mesh, mesh_dofs_info[0], dofs_data[0], count_dofs, reference_element_data, local_space_data, test, result);

        ComputeStrongTerm(c, mesh, mesh_dofs_info[0], dofs_data[0], count_dofs, reference_element_data, local_space_data, test, result);
    }

    result.rightHandSide.Create();
    result.solutionNeumann.Create();
    result.globalMatrixA.Create();
    result.neumannMatrixA.Create();

    if (count_dofs.num_total_strong > 0)
        result.rightHandSide.SubtractionMultiplication(result.neumannMatrixA, result.solutionNeumann);

    return result;
}
// ***************************************************************************
Assembler::Performance_Data Assembler::ComputePerformance(const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
                                                          const Gedim::MeshMatricesDAO &mesh,
                                                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                          const local_space::ReferenceElement_Data &reference_element_data) const
{
    Assembler::Performance_Data result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = local_space::CreateLocalSpace(config, mesh_geometric_data, c, reference_element_data);

        result.Cell2DsPerformance[c] = local_space::ComputePerformance(reference_element_data, local_space_data);
    }

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Polydim::examples::Elliptic_MCC_2D::Program_configuration &config,
                                                           const Gedim::MeshMatricesDAO &mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                           const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                                           const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                                           const local_space::ReferenceElement_Data &reference_element_data,
                                                           const Elliptic_MCC_2D_Problem_Data &assembler_data,
                                                           const Polydim::examples::Elliptic_MCC_2D::test::I_Test &test) const
{

    PostProcess_Data result;

    result.residual_norm = 0.0;
    if (count_dofs.num_total_dofs > 0)
    {
        Gedim::Eigen_Array<> residual;
        residual.SetSize(count_dofs.num_total_dofs);
        residual.SumMultiplication(assembler_data.globalMatrixA, assembler_data.solution);
        residual -= assembler_data.rightHandSide;

        result.residual_norm = residual.Norm();
    }

    result.cell2Ds_numeric_pressure.resize(mesh.Cell2DTotalNumber());
    result.cell2Ds_exact_pressure.resize(mesh.Cell2DTotalNumber());

    result.cell2Ds_error_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_super_error_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_L2_velocity.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2_velocity.setZero(mesh.Cell2DTotalNumber());
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = local_space::CreateLocalSpace(config, mesh_geometric_data, c, reference_element_data);

        const auto velocity_basis_functions_values =
            local_space::VelocityBasisFunctionsValues(reference_element_data, local_space_data);
        const auto pressure_basis_functions_values =
            local_space::PressureBasisFunctionsValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature = local_space::InternalQuadrature(reference_element_data, local_space_data);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(c).size();

        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                assembler_data.solution,
                                                                                assembler_data.solutionNeumann);

        const Eigen::VectorXd velocity_dofs_values =
            dofs_values.segment(0, local_count_dofs.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd pressure_dofs_values =
            dofs_values.segment(local_count_dofs.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        const auto exact_pressure_values = test.exact_pressure(cell2D_internal_quadrature.Points);
        const auto exact_velocity_values = test.exact_velocity(cell2D_internal_quadrature.Points);

        // Interpolate Exact Solution
        const VectorXd rightHandSide =
            pressure_basis_functions_values.transpose() *
            local_space_data.VEM_LocalSpace_Data_Pressure.InternalQuadrature.Weights.asDiagonal() * exact_pressure_values;
        const MatrixXd Hmatrix = pressure_basis_functions_values.transpose() *
                                 local_space_data.VEM_LocalSpace_Data_Pressure.InternalQuadrature.Weights.asDiagonal() *
                                 pressure_basis_functions_values;
        const VectorXd coeffPolynomial = Hmatrix.llt().solve(rightHandSide);

        result.cell2Ds_numeric_pressure[c] =
            cell2D_internal_quadrature.Weights.transpose() * (pressure_basis_functions_values * pressure_dofs_values);
        result.cell2Ds_exact_pressure[c] = cell2D_internal_quadrature.Weights.transpose() * exact_pressure_values;

        const Eigen::VectorXd local_error_L2_pressure =
            (pressure_basis_functions_values * pressure_dofs_values - exact_pressure_values).array().square();
        const Eigen::VectorXd local_super_error_L2_pressure =
            (pressure_basis_functions_values * pressure_dofs_values - pressure_basis_functions_values * coeffPolynomial)
                .array()
                .square();
        const Eigen::VectorXd local_norm_L2_pressure = (pressure_basis_functions_values * pressure_dofs_values).array().square();

        result.cell2Ds_error_L2_pressure[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_L2_pressure;
        result.cell2Ds_norm_L2_pressure[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_L2_pressure;

        const Eigen::VectorXd local_error_L2_velocity =
            (velocity_basis_functions_values[0] * velocity_dofs_values - exact_velocity_values[0]).array().square() +
            (velocity_basis_functions_values[1] * velocity_dofs_values - exact_velocity_values[1]).array().square();
        const Eigen::VectorXd local_norm_L2_velocity =
            (velocity_basis_functions_values[0] * velocity_dofs_values).array().square() +
            (velocity_basis_functions_values[1] * velocity_dofs_values).array().square();

        result.cell2Ds_super_error_L2_pressure[c] = cell2D_internal_quadrature.Weights.transpose() * local_super_error_L2_pressure;
        result.cell2Ds_error_L2_velocity[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_L2_velocity;
        result.cell2Ds_norm_L2_velocity[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_L2_velocity;

        if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2_pressure = std::sqrt(result.cell2Ds_error_L2_pressure.sum());
    result.super_error_L2_pressure = std::sqrt(result.cell2Ds_super_error_L2_pressure.sum());
    result.norm_L2_pressure = std::sqrt(result.cell2Ds_norm_L2_pressure.sum());
    result.error_L2_velocity = std::sqrt(result.cell2Ds_error_L2_velocity.sum());
    result.norm_L2_velocity = std::sqrt(result.cell2Ds_norm_L2_velocity.sum());

    return result;
}
// ***************************************************************************
} // namespace Elliptic_MCC_2D
} // namespace examples
} // namespace Polydim
