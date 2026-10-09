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

#include "Eigen_LUSolver.hpp"
#include "EllipticEquation.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "Quadrature_Gauss1D.hpp"
#include <cassert>

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{
// ***************************************************************************
void Assembler::ComputeStrongTerm(const unsigned int cell2D_index,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                  const std::vector<size_t> &offsetStrongs,
                                  const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                  const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                                  const test::I_Test &test,
                                  Stokes_DF_PCC_2D_Problem_Data &assembler_data) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int h = 0; h < 2; h++)
    {
        for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(cell2D_index); ++v)
        {
            const unsigned int cell0D_index = mesh.Cell2DVertex(cell2D_index, v);
            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(0).at(cell0D_index);

            if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
                continue;

            const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

            const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, coordinates).at(h);

            const auto local_dofs = dofs_data[h].CellsDOFs.at(0).at(cell0D_index);

            assert(local_dofs.size() == strong_boundary_values.size());

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                    assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index + offsetStrongs[h],
                                                              strong_boundary_values[loc_i]);
                }
                break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    continue;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }

        // Assemble strong boundary condition on Cell1Ds

        for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(cell2D_index); ++ed)
        {
            const unsigned int cell1D_index = mesh.Cell2DEdge(cell2D_index, ed);
            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(cell1D_index);

            const auto local_dofs = dofs_data[h].CellsDOFs.at(1).at(cell1D_index);

            if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
                local_dofs.size() == 0)
                continue;

            const auto edge_dofs_coordinates =
                Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityEdgeDofsCoordinates(reference_element_data, local_space_data, ed);

            const auto strong_boundary_values =
                test.strong_boundary_condition(boundary_info.Marker, edge_dofs_coordinates).at(h);

            assert(local_dofs.size() == strong_boundary_values.size());

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                    assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index + offsetStrongs[h],
                                                              strong_boundary_values[loc_i]);
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
}
// ***************************************************************************
void Assembler::ComputeWeakTerm(const unsigned int cell2D_index,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                                const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test,
                                Stokes_DF_PCC_2D_Problem_Data &assembler_data) const
{

    for (unsigned int h = 0; h < reference_element_data.Dimension; h++)
    {
        const unsigned numVertices = mesh_geometric_data.Cell2DsVertices.at(cell2D_index).cols();

        for (unsigned int ed = 0; ed < numVertices; ed++)
        {
            const unsigned int cell1D_index = mesh.Cell2DEdge(cell2D_index, ed);

            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(cell1D_index);

            if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
                continue;

            // compute vem values
            const auto weakReferenceSegment =
                Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

            const Eigen::VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0);

            // map edge internal quadrature points
            const Eigen::Vector3d &edgeStart = mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index)[ed]
                                                   ? mesh_geometric_data.Cell2DsVertices.at(cell2D_index).col(ed)
                                                   : mesh_geometric_data.Cell2DsVertices.at(cell2D_index).col((ed + 1) % numVertices);

            const Eigen::Vector3d &edgeTangent = mesh_geometric_data.Cell2DsEdgeTangents.at(cell2D_index).col(ed);
            const double direction = mesh_geometric_data.Cell2DsEdgeDirections.at(cell2D_index)[ed] ? 1.0 : -1.0;

            const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegment.Points.cols();
            Eigen::MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
            for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
                weakQuadraturePoints.col(q) = edgeStart + direction * weakReferenceSegment.Points(0, q) * edgeTangent;

            const double absMapDeterminant = std::abs(mesh_geometric_data.Cell2DsEdgeLengths.at(cell2D_index)[ed]);
            const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights * absMapDeterminant;

            const auto neumannValues = test.weak_boundary_condition(boundary_info.Marker, weakQuadraturePoints).at(h);
            const auto weak_basis_function_values =
                Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsValuesOnEdge(ed,
                                                                                            reference_element_data,
                                                                                            local_space_data,
                                                                                            pointsCurvilinearCoordinates);

            // compute values of Neumann condition
            const Eigen::VectorXd neumannContributions =
                -weak_basis_function_values.transpose() * weakQuadratureWeights.asDiagonal() * neumannValues;

            for (unsigned int p = 0; p < 2; ++p)
            {
                const unsigned int cell0D_index = mesh.Cell1DVertex(cell1D_index, p);

                const auto local_dofs = dofs_data[h].CellsDOFs.at(0).at(cell0D_index);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                    const auto &local_dof_i = local_dofs.at(loc_i);

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                        continue;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                        assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index + count_dofs.offsets_DOFs[h],
                                                              neumannContributions[p]);
                    }
                    break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }

            const auto local_dofs = dofs_data[h].CellsDOFs.at(1).at(cell1D_index);
            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                const unsigned int localIndex = loc_i;

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    continue;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index + count_dofs.offsets_DOFs[h],
                                                          neumannContributions[localIndex + 2]);
                }
                break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }
    }
}
// ***************************************************************************
Assembler::Stokes_DF_PCC_2D_Problem_Data Assembler::Assemble(
    const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const
{
    Stokes_DF_PCC_2D_Problem_Data result;
    result.globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    result.dirichletMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    result.rightHandSide.SetSize(count_dofs.num_total_dofs);
    result.solution.SetSize(count_dofs.num_total_dofs);
    result.solutionDirichlet.SetSize(count_dofs.num_total_strong);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                                config.GeometricTolerance2D(),
                                                                                                mesh_geometric_data,
                                                                                                c,
                                                                                                reference_element_data);

        const auto internal_quadrature =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const auto velocity_basis_functions_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsValues(reference_element_data,
                                                                                  local_space_data,
                                                                                  Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

        const auto velocity_basis_functions_divergence_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsDivergenceValues(reference_element_data, local_space_data);

        const auto velocity_basis_functions_derivatives_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsDerivativeValues(
                reference_element_data,
                local_space_data,
                Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

        const MatrixXd pressure_basis_functions_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::PressureBasisFunctionsValues(reference_element_data, local_space_data);

        const auto fluid_viscosity_values = test.fluid_viscosity(internal_quadrature.Points);
        const auto inverse_diffusion_term_values = test.inverse_diffusion_term(internal_quadrature.Points);
        const auto source_term_values = test.source_term(internal_quadrature.Points);
        const auto divergence_term_values = test.divergence_term(internal_quadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(fluid_viscosity_values,
                                                           velocity_basis_functions_derivatives_values,
                                                           internal_quadrature.Weights);

        const double mu_max = fluid_viscosity_values.cwiseAbs().maxCoeff();
        local_A += mu_max * Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityStabilizationMatrix(
                                reference_element_data,
                                local_space_data,
                                Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

        local_A += equation.ComputeCellDiffusionMatrix(inverse_diffusion_term_values,
                                                       velocity_basis_functions_values,
                                                       internal_quadrature.Weights);

        double k_max = 0.0;
        for (const auto &inv_diffusion_term : inverse_diffusion_term_values)
        {
            const double max_k = inv_diffusion_term.cwiseAbs().maxCoeff();
            k_max = k_max < max_k ? max_k : k_max;
        }
        local_A += k_max * Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityStabilizationMatrix(
                               reference_element_data,
                               local_space_data,
                               Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        internal_quadrature.Weights.asDiagonal() * velocity_basis_functions_divergence_values;

        const auto local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, velocity_basis_functions_values, internal_quadrature.Weights);

        const auto local_div =
            equation.ComputeCellForcingTerm(divergence_term_values, pressure_basis_functions_values, internal_quadrature.Weights);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[3].CellsGlobalDOFs[2].at(c).size();

        Eigen::MatrixXd elemental_matrix = MatrixXd::Zero(local_count_dofs.num_total_dofs, local_count_dofs.num_total_dofs);
        Eigen::VectorXd elemental_rhs = VectorXd::Zero(local_count_dofs.num_total_dofs);

        elemental_matrix << local_A, -local_B.transpose(), local_B, MatrixXd::Zero(num_local_dofs_pressure, num_local_dofs_pressure);

        elemental_rhs << local_rhs, local_div;

        assert(Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocitySize(reference_element_data, local_space_data) ==
               local_count_dofs.num_total_dofs - num_local_dofs_pressure);

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data = {
            {std::cref(dofs_data[0]), std::cref(dofs_data[1]), std::cref(dofs_data[2]), std::cref(dofs_data[3])},
            local_count_dofs.offsets_DOFs,
            count_dofs.offsets_DOFs,
            count_dofs.offsets_Strongs};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          elemental_matrix,
                                                                                          elemental_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.dirichletMatrixA,
                                                                                          result.rightHandSide);

        if (count_dofs.num_total_boundary_dofs == 0)
        {
            // Compute mean values
            const VectorXd mean_value_pressure = pressure_basis_functions_values.transpose() * internal_quadrature.Weights;

            // Mean value condition
            const unsigned int h1 = 3;
            const unsigned int num_global_offset_lagrange = count_dofs.num_total_dofs - 1;
            for (unsigned int loc_i = 0; loc_i < dofs_data[h1].CellsGlobalDOFs[2].at(c).size(); loc_i++)
            {
                const auto global_dof_i = dofs_data[h1].CellsGlobalDOFs[2].at(c).at(loc_i);
                const auto local_dof_i =
                    dofs_data[h1].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);
                const unsigned int global_index_i = local_dof_i.Global_Index + count_dofs.offsets_DOFs[h1];

                result.globalMatrixA.Triplet(global_index_i, num_global_offset_lagrange, mean_value_pressure(loc_i));

                result.globalMatrixA.Triplet(num_global_offset_lagrange, global_index_i, mean_value_pressure(loc_i));
            }
        }

        ComputeStrongTerm(c, mesh, mesh_dofs_info, dofs_data, count_dofs.offsets_Strongs, reference_element_data, local_space_data, test, result);

        ComputeWeakTerm(c, mesh, mesh_geometric_data, mesh_dofs_info, dofs_data, count_dofs, reference_element_data, local_space_data, test, result);
    }

    result.rightHandSide.Create();
    result.solutionDirichlet.Create();
    result.globalMatrixA.Create();
    result.dirichletMatrixA.Create();

    if (count_dofs.num_total_strong > 0)
        result.rightHandSide.SubtractionMultiplication(result.dirichletMatrixA, result.solutionDirichlet);

    return result;
}
// ***************************************************************************
Assembler::Performance_Data Assembler::ComputeMethodPerformance(
    const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data) const
{
    Assembler::Performance_Data result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                                config.GeometricTolerance2D(),
                                                                                                mesh_geometric_data,
                                                                                                c,
                                                                                                reference_element_data);

        result.Cell2DsPerformance[c] =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::ComputePerformance(reference_element_data, local_space_data);
    }

    return result;
}
// ***************************************************************************
Assembler::DiscrepancyErrors_Data Assembler::ComputeDiscrepancyErrors(
    const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const vector<PDETools::DOFs::DOFsManager::DOFsData> &reduced_dofs_data,
    const Polydim::PDETools::Assembler_Utilities::count_dofs_data &reduced_count_dofs,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reduced_reference_element_data,
    const Stokes_DF_PCC_2D_Problem_Data &reduced_assembler_data,
    const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const
{

    const auto full_reference_element_data = Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateReferenceElement(
        Polydim::PDETools::LocalSpace_DF_PCC_2D::MethodTypes::VEM_DF_PCC_FULL,
        config.MethodOrder());

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);

    const auto full_mesh_dofs_info =
        Polydim::PDETools::LocalSpace_DF_PCC_2D::SetMeshDOFsInfo(full_reference_element_data, mesh, test.boundary_info());
    const unsigned int full_num_mesh_dofs_info = full_mesh_dofs_info.size();
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> full_dofs_data(full_num_mesh_dofs_info);

    for (unsigned int i = 0; i < full_num_mesh_dofs_info; i++)
        full_dofs_data[i] = dofManager.CreateDOFs_2D(full_mesh_dofs_info[i], mesh_connectivity_data);

    auto full_count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(full_dofs_data);
    if (full_count_dofs.num_total_boundary_dofs == 0)
        full_count_dofs.num_total_dofs += 1; // lagrange

    auto full_assembler_data =
        Assemble(config, mesh, mesh_geometric_data, full_mesh_dofs_info, full_dofs_data, full_count_dofs, full_reference_element_data, test);

    if (full_count_dofs.num_total_dofs > 0)
    {
        Gedim::Eigen_LUSolver solver;
        solver.Initialize(full_assembler_data.globalMatrixA);

        solver.Solve(full_assembler_data.rightHandSide, full_assembler_data.solution);
    }

    DiscrepancyErrors_Data result;

    result.residual_norm = 0.0;
    if (full_count_dofs.num_total_dofs > 0)
    {
        Gedim::Eigen_Array<> residual;
        residual.SetSize(full_count_dofs.num_total_dofs);
        residual.SumMultiplication(full_assembler_data.globalMatrixA, full_assembler_data.solution);
        residual -= full_assembler_data.rightHandSide;

        result.residual_norm = residual.Norm();
    }

    result.reduced_residual_norm = 0.0;
    if (reduced_count_dofs.num_total_dofs > 0)
    {
        Gedim::Eigen_Array<> residual;
        residual.SetSize(reduced_count_dofs.num_total_dofs);
        residual.SumMultiplication(reduced_assembler_data.globalMatrixA, reduced_assembler_data.solution);
        residual -= reduced_assembler_data.rightHandSide;

        result.reduced_residual_norm = residual.Norm();
    }

    result.cell2Ds_discrepancy_error_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_discrepancy_error_H1_velocity.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_full_norm_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_full_norm_H1_velocity.setZero(mesh.Cell2DTotalNumber());

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto reduced_local_space_data =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                      config.GeometricTolerance2D(),
                                                                      mesh_geometric_data,
                                                                      c,
                                                                      reduced_reference_element_data);

        const auto reduced_internal_quadrature =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::InternalQuadrature(reduced_reference_element_data, reduced_local_space_data);

        const auto reduced_velocity_basis_functions_derivatives_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsDerivativeValues(
                reduced_reference_element_data,
                reduced_local_space_data,
                Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);

        const MatrixXd reduced_pressure_basis_functions_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::PressureBasisFunctionsValues(reduced_reference_element_data,
                                                                                  reduced_local_space_data);

        const auto full_local_space_data =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                      config.GeometricTolerance2D(),
                                                                      mesh_geometric_data,
                                                                      c,
                                                                      full_reference_element_data);

        const auto full_internal_quadrature =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::InternalQuadrature(full_reference_element_data, full_local_space_data);

        const auto full_velocity_basis_functions_derivatives_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsDerivativeValues(
                full_reference_element_data,
                full_local_space_data,
                Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);

        const MatrixXd full_pressure_basis_functions_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::PressureBasisFunctionsValues(full_reference_element_data, full_local_space_data);

        const auto full_local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, full_dofs_data);
        const unsigned int full_num_local_dofs_pressure = full_dofs_data[3].CellsGlobalDOFs[2].at(c).size();

        const Eigen::VectorXd full_dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                full_dofs_data,
                                                                                full_local_count_dofs.num_total_dofs,
                                                                                full_local_count_dofs.offsets_DOFs,
                                                                                full_count_dofs.offsets_DOFs,
                                                                                full_count_dofs.offsets_Strongs,
                                                                                full_assembler_data.solution,
                                                                                full_assembler_data.solutionDirichlet);

        const Eigen::VectorXd full_velocity_dofs_values =
            full_dofs_values.segment(0, full_local_count_dofs.num_total_dofs - full_num_local_dofs_pressure);
        const Eigen::VectorXd full_pressure_dofs_values =
            full_dofs_values.segment(full_local_count_dofs.num_total_dofs - full_num_local_dofs_pressure, full_num_local_dofs_pressure);

        const auto reduced_local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, reduced_dofs_data);
        const unsigned int reduced_num_local_dofs_pressure = reduced_dofs_data[3].CellsGlobalDOFs[2].at(c).size();

        const Eigen::VectorXd reduced_dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                reduced_dofs_data,
                                                                                reduced_local_count_dofs.num_total_dofs,
                                                                                reduced_local_count_dofs.offsets_DOFs,
                                                                                reduced_count_dofs.offsets_DOFs,
                                                                                reduced_count_dofs.offsets_Strongs,
                                                                                reduced_assembler_data.solution,
                                                                                reduced_assembler_data.solutionDirichlet);

        const Eigen::VectorXd reduced_velocity_dofs_values =
            reduced_dofs_values.segment(0, reduced_local_count_dofs.num_total_dofs - reduced_num_local_dofs_pressure);
        const Eigen::VectorXd reduced_pressure_dofs_values =
            reduced_dofs_values.segment(reduced_local_count_dofs.num_total_dofs - reduced_num_local_dofs_pressure,
                                        reduced_num_local_dofs_pressure);

        const VectorXd full_numeric_pressure_values = full_pressure_basis_functions_values * full_pressure_dofs_values;
        const VectorXd reduced_numeric_pressure_values = reduced_pressure_basis_functions_values * reduced_pressure_dofs_values;

        const VectorXd full_numeric_projected_pressure_values =
            ((1.0 / mesh_geometric_data.Cell2DsAreas[c]) * full_numeric_pressure_values.transpose() *
             full_internal_quadrature.Weights) *
            reduced_pressure_basis_functions_values;

        const Eigen::VectorXd local_error_L2_pressure =
            (full_numeric_projected_pressure_values - reduced_numeric_pressure_values).array().square();
        const Eigen::VectorXd local_norm_L2_pressure = (full_numeric_projected_pressure_values).array().square();

        result.cell2Ds_discrepancy_error_L2_pressure[c] = full_internal_quadrature.Weights.transpose() * local_error_L2_pressure;
        result.cell2Ds_full_norm_L2_pressure[c] = full_internal_quadrature.Weights.transpose() * local_norm_L2_pressure;

        const unsigned int numQuadraturePoints = full_internal_quadrature.Points.cols();
        Eigen::VectorXd local_error_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        Eigen::VectorXd local_norm_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        for (unsigned int d1 = 0; d1 < full_reference_element_data.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < full_reference_element_data.Dimension; d2++)
            {
                local_error_H1_velocity.array() +=
                    (full_velocity_basis_functions_derivatives_values[full_reference_element_data.Dimension * d1 + d2] * full_velocity_dofs_values -
                     reduced_velocity_basis_functions_derivatives_values[full_reference_element_data.Dimension * d1 + d2] * reduced_velocity_dofs_values)
                        .array()
                        .square();

                local_norm_H1_velocity.array() +=
                    (full_velocity_basis_functions_derivatives_values[full_reference_element_data.Dimension * d1 + d2] * full_velocity_dofs_values)
                        .array()
                        .square();
            }
        }
        result.cell2Ds_discrepancy_error_H1_velocity[c] = full_internal_quadrature.Weights.transpose() * local_error_H1_velocity;
        result.cell2Ds_full_norm_H1_velocity[c] = full_internal_quadrature.Weights.transpose() * local_norm_H1_velocity;
    }

    result.discrepancy_error_L2_pressure = std::sqrt(result.cell2Ds_discrepancy_error_L2_pressure.sum());
    result.discrepancy_error_H1_velocity = std::sqrt(result.cell2Ds_discrepancy_error_H1_velocity.sum());
    result.full_norm_H1_velocity = std::sqrt(result.cell2Ds_full_norm_H1_velocity.sum());
    result.full_norm_L2_pressure = std::sqrt(result.cell2Ds_full_norm_L2_pressure.sum());

    result.pressure_dofs_ratio = ((double)reduced_dofs_data[3].NumberDOFs) / full_dofs_data[3].NumberDOFs;
    result.velocity_dofs_ratio = ((double)reduced_count_dofs.offsets_DOFs[3]) / full_count_dofs.offsets_DOFs[3];

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(
    const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Stokes_DF_PCC_2D_Problem_Data &assembler_data,
    const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test) const
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

    for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
    {
        result.cell0Ds_numeric_velocity[d].resize(mesh.Cell0DTotalNumber());
        result.cell0Ds_exact_velocity[d].resize(mesh.Cell0DTotalNumber());
    }

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        try
        {
            const auto exact_velocity = test.exact_velocity(mesh.Cell0DCoordinates(p));

            for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
                result.cell0Ds_exact_velocity[d](p) = exact_velocity[d](0);
        }
        catch (...)
        {
            for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
                result.cell0Ds_exact_velocity[d](p) = std::nan("");
        }

        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
        {
            const auto local_dofs = dofs_data[d].CellsDOFs.at(0).at(p);

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    result.cell0Ds_numeric_velocity[d][p] =
                        assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index + count_dofs.offsets_Strongs[d]);
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    result.cell0Ds_numeric_velocity[d][p] =
                        assembler_data.solution.GetValue(local_dof_i.Global_Index + count_dofs.offsets_DOFs[d]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }
    }

    result.repeated_connectivity = mesh.Cell2DsVertices();
    unsigned int num_repetead_points = 0;
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(c); v++)
        {
            result.repeated_connectivity[c][v] = num_repetead_points;
            num_repetead_points++;
        }
    }
    result.repeated_vertices_coordinates = Eigen::MatrixXd::Zero(3, num_repetead_points);
    result.cell0Ds_exact_pressure = Eigen::VectorXd::Constant(num_repetead_points, std::nan(""));
    result.cell0Ds_numeric_pressure = Eigen::VectorXd::Zero(num_repetead_points);

    result.cell2Ds_error_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_H1_velocity.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_H1_velocity.setZero(mesh.Cell2DTotalNumber());
    result.error_L2_pressure = 0.0;
    result.norm_L2_pressure = 0.0;
    result.error_H1_velocity = 0.0;
    result.norm_H1_velocity = 0.0;
    result.mesh_size = 0.0;

    result.inverse_diffusion_coeff_values.setZero(mesh.Cell2DTotalNumber());
    result.viscosity_values.setZero(mesh.Cell2DTotalNumber());

    result.flux = ComputeFlux(config, mesh, mesh_geometric_data, dofs_data, count_dofs, reference_element_data, test, assembler_data);

    num_repetead_points = 0;
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                                config.GeometricTolerance2D(),
                                                                                                mesh_geometric_data,
                                                                                                c,
                                                                                                reference_element_data);

        const auto internal_quadrature =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const auto velocity_basis_functions_derivatives_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsDerivativeValues(
                reference_element_data,
                local_space_data,
                Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);

        const MatrixXd pressure_basis_functions_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::PressureBasisFunctionsValues(reference_element_data, local_space_data);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[3].CellsGlobalDOFs[2].at(c).size();

        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                assembler_data.solution,
                                                                                assembler_data.solutionDirichlet);

        result.viscosity_values[c] = test.fluid_viscosity(mesh_geometric_data.Cell2DsCentroids.at(c))[0];
        const auto inverse_diffusion_term_values = test.inverse_diffusion_term(mesh_geometric_data.Cell2DsCentroids.at(c));
        double k_max = 0.0;
        for (const auto &inv_diffusion_term : inverse_diffusion_term_values)
        {
            const double max_k = inv_diffusion_term.cwiseAbs().maxCoeff();
            k_max = k_max < max_k ? max_k : k_max;
        }
        result.inverse_diffusion_coeff_values[c] = k_max;

        const Eigen::VectorXd velocity_dofs_values =
            dofs_values.segment(0, local_count_dofs.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd pressure_dofs_values =
            dofs_values.segment(local_count_dofs.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        const unsigned int num_cell_vertices = mesh.Cell2DNumberVertices(c);
        const Eigen::MatrixXd pressure_basis_functions_vertices_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::PressureBasisFunctionsValues(reference_element_data,
                                                                                  local_space_data,
                                                                                  mesh_geometric_data.Cell2DsVertices.at(c));
        result.repeated_vertices_coordinates.middleCols(num_repetead_points, num_cell_vertices) =
            mesh_geometric_data.Cell2DsVertices.at(c);
        result.cell0Ds_numeric_pressure.segment(num_repetead_points, num_cell_vertices) =
            pressure_basis_functions_vertices_values * pressure_dofs_values;

        try
        {
            result.cell0Ds_exact_pressure.segment(num_repetead_points, num_cell_vertices) =
                test.exact_pressure(mesh_geometric_data.Cell2DsVertices.at(c));
        }
        catch (...)
        {
        }

        num_repetead_points += num_cell_vertices;

        const VectorXd numeric_pressure_values = pressure_basis_functions_values * pressure_dofs_values;

        try
        {
            const VectorXd exact_pressure_values = test.exact_pressure(internal_quadrature.Points);
            const auto exact_velocity_derivatives_values = test.exact_derivatives_velocity(internal_quadrature.Points);

            const Eigen::VectorXd local_error_L2_pressure = (numeric_pressure_values - exact_pressure_values).array().square();
            const Eigen::VectorXd local_norm_L2_pressure = (numeric_pressure_values).array().square();

            result.cell2Ds_error_L2_pressure[c] = internal_quadrature.Weights.transpose() * local_error_L2_pressure;
            result.cell2Ds_norm_L2_pressure[c] = internal_quadrature.Weights.transpose() * local_norm_L2_pressure;

            const unsigned int numQuadraturePoints = internal_quadrature.Points.cols();
            Eigen::VectorXd local_error_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
            Eigen::VectorXd local_norm_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
            for (unsigned int d1 = 0; d1 < reference_element_data.Dimension; d1++)
            {
                for (unsigned int d2 = 0; d2 < reference_element_data.Dimension; d2++)
                {
                    local_error_H1_velocity.array() +=
                        (velocity_basis_functions_derivatives_values[reference_element_data.Dimension * d1 + d2] * velocity_dofs_values -
                         exact_velocity_derivatives_values[3 * d1 + d2])
                            .array()
                            .square();

                    local_norm_H1_velocity.array() +=
                        (velocity_basis_functions_derivatives_values[reference_element_data.Dimension * d1 + d2] * velocity_dofs_values)
                            .array()
                            .square();
                }
            }
            result.cell2Ds_error_H1_velocity[c] = internal_quadrature.Weights.transpose() * local_error_H1_velocity;
            result.cell2Ds_norm_H1_velocity[c] = internal_quadrature.Weights.transpose() * local_norm_H1_velocity;
        }
        catch (...)
        {
            const Eigen::VectorXd local_norm_L2_pressure = (numeric_pressure_values).array().square();
            result.cell2Ds_norm_L2_pressure[c] = internal_quadrature.Weights.transpose() * local_norm_L2_pressure;

            const unsigned int numQuadraturePoints = internal_quadrature.Points.cols();
            Eigen::VectorXd local_norm_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
            for (unsigned int d1 = 0; d1 < reference_element_data.Dimension; d1++)
            {
                for (unsigned int d2 = 0; d2 < reference_element_data.Dimension; d2++)
                {
                    local_norm_H1_velocity.array() +=
                        (velocity_basis_functions_derivatives_values[reference_element_data.Dimension * d1 + d2] * velocity_dofs_values)
                            .array()
                            .square();
                }
            }
            result.cell2Ds_norm_H1_velocity[c] = internal_quadrature.Weights.transpose() * local_norm_H1_velocity;
        }

        if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2_pressure = std::sqrt(result.cell2Ds_error_L2_pressure.sum());
    result.norm_L2_pressure = std::sqrt(result.cell2Ds_norm_L2_pressure.sum());
    result.error_H1_velocity = std::sqrt(result.cell2Ds_error_H1_velocity.sum());
    result.norm_H1_velocity = std::sqrt(result.cell2Ds_norm_H1_velocity.sum());

    return result;
}
// ***************************************************************************
std::map<unsigned int, double> Assembler::ComputeFlux(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                                                      const Gedim::MeshMatricesDAO &mesh,
                                                      const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                      const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                                      const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                                      const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                      const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test,
                                                      const Stokes_DF_PCC_2D_Problem_Data &assembler_data) const
{

    std::map<unsigned int, double> flux;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                                config.GeometricTolerance2D(),
                                                                                                mesh_geometric_data,
                                                                                                c,
                                                                                                reference_element_data);

        const unsigned numVertices = mesh_geometric_data.Cell2DsVertices.at(c).cols();

        for (unsigned int ed = 0; ed < numVertices; ed++)
        {
            const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

            const unsigned int edge_marker = mesh.Cell1DMarker(cell1D_index);

            if (edge_marker == 0)
                continue;

            // compute vem values
            const auto weakReferenceSegment =
                Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

            const Eigen::VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0);

            const auto weak_basis_function_values =
                Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsValuesOnEdge(ed,
                                                                                            reference_element_data,
                                                                                            local_space_data,
                                                                                            pointsCurvilinearCoordinates);

            const unsigned int num_edge_dofs = weak_basis_function_values.cols();

            Eigen::VectorXd local_solution_dofs = Eigen::VectorXd::Zero(num_edge_dofs * 2);
            unsigned int offset_dofs = 0;
            for (unsigned int h = 0; h < 2; h++)
            {
                const auto &global_dof = dofs_data[h].CellsGlobalDOFs[1].at(cell1D_index);
                for (unsigned int loc_i = 0; loc_i < global_dof.size(); ++loc_i)
                {
                    const auto &global_dof_i = global_dof.at(loc_i);
                    const auto &local_dof_i =
                        dofs_data[h].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                        local_solution_dofs[loc_i + offset_dofs] =
                            assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index + count_dofs.offsets_Strongs[h]);
                        break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                        local_solution_dofs[loc_i + offset_dofs] =
                            assembler_data.solution.GetValue(local_dof_i.Global_Index + count_dofs.offsets_DOFs[h]);
                        break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }

                offset_dofs += num_edge_dofs;
            }

            const double absMapDeterminant = std::abs(mesh_geometric_data.Cell2DsEdgeLengths.at(c)[ed]);
            const Eigen::VectorXd weakQuadratureWeights = weakReferenceSegment.Weights * absMapDeterminant;

            const Eigen::VectorXd edge_normal = mesh_geometric_data.Cell2DsEdgeNormals.at(c).col(ed);
            const Eigen::VectorXd local_edge_flux =
                edge_normal[0] * weak_basis_function_values * local_solution_dofs.segment(0, num_edge_dofs) +
                edge_normal[1] * weak_basis_function_values * local_solution_dofs.segment(num_edge_dofs, num_edge_dofs);
            const double edge_flux = local_edge_flux.transpose() * weakQuadratureWeights;

            auto it = flux.find(edge_marker);
            if (it != flux.end())
                it->second += edge_flux;
            else
                flux.insert({edge_marker, edge_flux});
        }
    }

    return flux;
}
// ***************************************************************************
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim
