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

#include "Assembler_Utilities.hpp"
#include "EllipticEquation.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "VEM_PCC_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elastic_PCC_2D
{
//***************************************************************************
void Assembler::ComputeStrongTerm(const unsigned int cell2D_index,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                  const local_space::ReferenceElement_Data &reference_element_data,
                                  const local_space::LocalSpace_Data &local_space_data,
                                  const test::I_Test &test,
                                  Elastic_PCC_2D_Problem_Data &assembler_data) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(cell2D_index); ++v)
    {
        const unsigned int cell0D_index = mesh.Cell2DVertex(cell2D_index, v);
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
            continue;

        const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, coordinates);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

        for (unsigned int h = 0; h < reference_element_data.Dimension; h++)
        {
            assert(local_dofs.size() == strong_boundary_values[h].size());

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                    assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index + count_dofs.offsets_Strongs[h],
                                                              strong_boundary_values[h][loc_i]);
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

    // Assemble strong boundary condition on Cell1Ds
    for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(cell2D_index); ++ed)
    {
        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2D_index, ed);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);
        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
            local_dofs.size() == 0)
            continue;

        const auto edge_dofs_coordinates = local_space::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, edge_dofs_coordinates);

        for (unsigned int h = 0; h < reference_element_data.Dimension; h++)
        {
            assert(local_dofs.size() == strong_boundary_values[h].size());

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                    assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index + count_dofs.offsets_Strongs[h],
                                                              strong_boundary_values[h][loc_i]);
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
void Assembler::ComputeWeakTerm(const unsigned int cell2DIndex,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                const local_space::ReferenceElement_Data &reference_element_data,
                                const local_space::LocalSpace_Data &local_space_data,
                                const Polydim::examples::Elastic_PCC_2D::test::I_Test &test,
                                Elastic_PCC_2D_Problem_Data &assembler_data) const
{
    const unsigned numVertices = mesh_geometric_data.Cell2DsVertices.at(cell2DIndex).cols();

    for (unsigned int ed = 0; ed < numVertices; ed++)
    {
        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex, ed);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
            continue;

        // compute vem values
        const auto weakReferenceSegment =
            Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

        const Eigen::VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0);

        // map edge internal quadrature points
        const Eigen::Vector3d &edgeStart = mesh_geometric_data.Cell2DsEdgeDirections.at(cell2DIndex)[ed]
                                               ? mesh_geometric_data.Cell2DsVertices.at(cell2DIndex).col(ed)
                                               : mesh_geometric_data.Cell2DsVertices.at(cell2DIndex).col((ed + 1) % numVertices);

        const Eigen::Vector3d &edgeTangent = mesh_geometric_data.Cell2DsEdgeTangents.at(cell2DIndex).col(ed);
        const double direction = mesh_geometric_data.Cell2DsEdgeDirections.at(cell2DIndex)[ed] ? 1.0 : -1.0;

        const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegment.Points.cols();
        Eigen::MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
        for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
            weakQuadraturePoints.col(q) = edgeStart + direction * weakReferenceSegment.Points(0, q) * edgeTangent;

        const double absMapDeterminant = std::abs(mesh_geometric_data.Cell2DsEdgeLengths.at(cell2DIndex)[ed]);
        const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights * absMapDeterminant;

        const auto neumannValues = test.weak_boundary_condition(boundary_info.Marker, weakQuadraturePoints);
        const auto weak_basis_function_values =
            local_space::BasisFunctionsValuesOnEdges(ed, reference_element_data, local_space_data, pointsCurvilinearCoordinates);

        for (unsigned int h = 0; h < 2; h++)
        {
            // compute values of Neumann condition
            const Eigen::VectorXd neumannContributions =
                weak_basis_function_values.transpose() * weakQuadratureWeights.asDiagonal() * neumannValues[h];

            for (unsigned int p = 0; p < 2; ++p)
            {
                const unsigned int cell0D_index = mesh.Cell1DVertex(cell1D_index, p);

                const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

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

            const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);
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
Assembler::Elastic_PCC_2D_Problem_Data Assembler::Assemble(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
                                                           const Gedim::MeshMatricesDAO &mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                           const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                                           const local_space::ReferenceElement_Data &reference_element_data,
                                                           const Polydim::examples::Elastic_PCC_2D::test::I_Test &test) const
{
    Elastic_PCC_2D_Problem_Data result;

    result.globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
    result.dirichletMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    result.rightHandSide.SetSize(count_dofs.num_total_dofs);
    result.solution.SetSize(count_dofs.num_total_dofs);
    result.solutionDirichlet.SetSize(count_dofs.num_total_strong);

    Polydim::PDETools::Equations::EllipticEquation equation;

    std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> dofs_data_ref;
    dofs_data_ref.reserve(2);
    dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data)));
    dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data)));

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const auto local_space_data = local_space::CreateLocalSpace(config, mesh_geometric_data, c, reference_element_data);

        const auto basis_functions_values = local_space::BasisFunctionsValues(reference_element_data, local_space_data);

        const auto basis_functions_derivative_values =
            local_space::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        std::vector<Eigen::MatrixXd> symmetric_gardient(4);
        for (unsigned int d1 = 0; d1 < 2; d1++)
            for (unsigned int d2 = 0; d2 < 2; d2++)
                symmetric_gardient[2 * d1 + d2] =
                    0.5 * (basis_functions_derivative_values[2 * d1 + d2] + basis_functions_derivative_values[2 * d2 + d1]);

        const auto cell2D_internal_quadrature = local_space::InternalQuadrature(reference_element_data, local_space_data);

        const auto lame_coefficients_values = test.lame_coefficients(cell2D_internal_quadrature.Points);
        const auto source_term_values = test.source_term(cell2D_internal_quadrature.Points);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data_ref);

        const Eigen::MatrixXd local_A =
            2.0 * equation.ComputeCellDiffusionMatrix(lame_coefficients_values[0],
                                                      symmetric_gardient,
                                                      cell2D_internal_quadrature.Weights) +
            equation.ComputeCellReactionMatrix(lame_coefficients_values[1],
                                               basis_functions_derivative_values[0] + basis_functions_derivative_values[3],
                                               cell2D_internal_quadrature.Weights);

        const Eigen::VectorXd local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        const Eigen::MatrixXd local_A_stab = lame_coefficients_values[0].cwiseAbs().maxCoeff() *
                                             local_space::StabilizationMatrix(reference_element_data, local_space_data);

        {
            const auto &global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);
            assert(local_space::Size(reference_element_data, local_space_data) == global_dofs.size());
        }

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data = {
            {std::cref(dofs_data), std::cref(dofs_data)},
            local_count_dofs.offsets_DOFs,
            count_dofs.offsets_DOFs,
            count_dofs.offsets_Strongs};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A + local_A_stab,
                                                                                          local_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.dirichletMatrixA,
                                                                                          result.rightHandSide);

        ComputeStrongTerm(c, mesh, mesh_dofs_info, dofs_data, count_dofs, reference_element_data, local_space_data, test, result);

        ComputeWeakTerm(c, mesh, mesh_geometric_data, mesh_dofs_info, dofs_data, count_dofs, reference_element_data, local_space_data, test, result);
    }

    result.rightHandSide.Create();
    result.solutionDirichlet.Create();
    result.globalMatrixA.Create();
    result.dirichletMatrixA.Create();

    if (dofs_data.NumberStrongs > 0)
        result.rightHandSide.SubtractionMultiplication(result.dirichletMatrixA, result.solutionDirichlet);

    return result;
}
// ***************************************************************************
Assembler::Performance_Data Assembler::ComputePerformance(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
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
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Polydim::examples::Elastic_PCC_2D::Program_configuration &config,
                                                           const Gedim::MeshMatricesDAO &mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                           const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                                           const local_space::ReferenceElement_Data &reference_element_data,
                                                           const Elastic_PCC_2D_Problem_Data &assembler_data,
                                                           const Polydim::examples::Elastic_PCC_2D::test::I_Test &test) const
{
    std::vector<std::reference_wrapper<const Polydim::PDETools::DOFs::DOFsManager::DOFsData>> dofs_data_ref;
    dofs_data_ref.reserve(2);
    dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data)));
    dofs_data_ref.push_back(std::cref(static_cast<const Polydim::PDETools::DOFs::DOFsManager::DOFsData &>(dofs_data)));

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
        result.cell0Ds_numeric_displacement[d].resize(mesh.Cell0DTotalNumber());
        result.cell0Ds_exact_displacement[d].resize(mesh.Cell0DTotalNumber());
    }

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {

        try
        {
            const auto exact_displacement = test.exact_displacement(mesh.Cell0DCoordinates(p));

            for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
                result.cell0Ds_exact_displacement[d](p) = exact_displacement[d](0);
        }
        catch (...)
        {
            for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
                result.cell0Ds_exact_displacement[d](p) = std::nan("");
        }

        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
        {
            const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    result.cell0Ds_numeric_displacement[d][p] =
                        assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index + count_dofs.offsets_Strongs[d]);
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    result.cell0Ds_numeric_displacement[d][p] =
                        assembler_data.solution.GetValue(local_dof_i.Global_Index + count_dofs.offsets_DOFs[d]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }
    }

    result.cell2Ds_error_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_H1.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_H1.setZero(mesh.Cell2DTotalNumber());
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = local_space::CreateLocalSpace(config, mesh_geometric_data, c, reference_element_data);

        const auto basis_functions_values = local_space::BasisFunctionsValues(reference_element_data, local_space_data);

        const auto basis_functions_derivative_values =
            local_space::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature = local_space::InternalQuadrature(reference_element_data, local_space_data);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data_ref);

        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data_ref,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                assembler_data.solution,
                                                                                assembler_data.solutionDirichlet);

        const unsigned int numQuadraturePoints = cell2D_internal_quadrature.Points.cols();
        Eigen::VectorXd local_error_L2_displacement = Eigen::VectorXd::Zero(numQuadraturePoints);
        Eigen::VectorXd local_norm_L2_displacement = Eigen::VectorXd::Zero(numQuadraturePoints);
        Eigen::VectorXd local_error_H1_displacement = Eigen::VectorXd::Zero(numQuadraturePoints);
        Eigen::VectorXd local_norm_H1_displacement = Eigen::VectorXd::Zero(numQuadraturePoints);

        try
        {
            const auto exact_solution_values = test.exact_displacement(cell2D_internal_quadrature.Points);
            const auto exact_derivative_displacement_values =
                test.exact_derivatives_displacement(cell2D_internal_quadrature.Points);

            for (unsigned int d1 = 0; d1 < reference_element_data.Dimension; d1++)
            {
                local_error_L2_displacement.array() +=
                    (basis_functions_values[d1] * dofs_values - exact_solution_values[d1]).array().square();
                local_norm_L2_displacement.array() += (basis_functions_values[d1] * dofs_values).array().square();

                for (unsigned int d2 = 0; d2 < reference_element_data.Dimension; d2++)
                {
                    local_error_H1_displacement.array() +=
                        (basis_functions_derivative_values[reference_element_data.Dimension * d1 + d2] * dofs_values -
                         exact_derivative_displacement_values[3 * d1 + d2])
                            .array()
                            .square();

                    local_norm_H1_displacement.array() +=
                        (basis_functions_derivative_values[reference_element_data.Dimension * d1 + d2] * dofs_values)
                            .array()
                            .square();
                }
            }

            result.cell2Ds_error_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_L2_displacement;
            result.cell2Ds_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_L2_displacement;

            result.cell2Ds_error_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_H1_displacement;
            result.cell2Ds_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_H1_displacement;
        }
        catch (...)
        {

            for (unsigned int d1 = 0; d1 < reference_element_data.Dimension; d1++)
            {

                local_norm_L2_displacement.array() += (basis_functions_values[d1] * dofs_values).array().square();

                for (unsigned int d2 = 0; d2 < reference_element_data.Dimension; d2++)
                    local_norm_H1_displacement.array() +=
                        (basis_functions_derivative_values[reference_element_data.Dimension * d1 + d2] * dofs_values)
                            .array()
                            .square();
            }

            result.cell2Ds_error_L2[c] = std::nan("");
            result.cell2Ds_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_L2_displacement;

            result.cell2Ds_error_H1[c] = std::nan("");
            result.cell2Ds_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_H1_displacement;
        }

        if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2 = std::sqrt(result.cell2Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell2Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell2Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell2Ds_norm_H1.sum());

    return result;
}
} // namespace Elastic_PCC_2D
} // namespace examples
} // namespace Polydim
