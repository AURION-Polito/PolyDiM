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
#include "Quadrature_Gauss1D.hpp"
#include <cassert>

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace examples
{
namespace NavierStokes_DF_PCC_2D
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
                                  NavierStokes_DF_PCC_2D_Problem_Data &assembler_data) const
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
                                const Polydim::examples::NavierStokes_DF_PCC_2D::test::I_Test &test,
                                NavierStokes_DF_PCC_2D_Problem_Data &assembler_data) const
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
Assembler::NavierStokes_DF_PCC_2D_Problem_Data Assembler::AssembleStokes(
    const Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::examples::NavierStokes_DF_PCC_2D::test::I_Test &test) const
{
    NavierStokes_DF_PCC_2D_Problem_Data result;
    result.globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    result.dirichletMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    result.rightHandSide.SetSize(count_dofs.num_total_dofs);
    result.solution.SetSize(count_dofs.num_total_dofs);
    result.previousIteration.SetSize(count_dofs.num_total_dofs);
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
        const auto source_term_values = test.source_term(internal_quadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(fluid_viscosity_values,
                                                           velocity_basis_functions_derivatives_values,
                                                           internal_quadrature.Weights);

        const double mu_max = fluid_viscosity_values.cwiseAbs().maxCoeff();
        local_A += mu_max * Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityStabilizationMatrix(
                                reference_element_data,
                                local_space_data,
                                Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        internal_quadrature.Weights.asDiagonal() * velocity_basis_functions_divergence_values;

        const auto local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, velocity_basis_functions_values, internal_quadrature.Weights);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[3].CellsGlobalDOFs[2].at(c).size();

        Eigen::MatrixXd elemental_matrix = MatrixXd::Zero(local_count_dofs.num_total_dofs, local_count_dofs.num_total_dofs);
        Eigen::VectorXd elemental_rhs = VectorXd::Zero(local_count_dofs.num_total_dofs);

        elemental_matrix << local_A, local_B.transpose(), local_B, MatrixXd::Zero(num_local_dofs_pressure, num_local_dofs_pressure);

        elemental_rhs << local_rhs, VectorXd::Zero(num_local_dofs_pressure);

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
Eigen::MatrixXd Assembler::ComputeConvectiveMatrix(const std::vector<Eigen::VectorXd> &previous_iteration_values,
                                                   const std::vector<Eigen::VectorXd> &previous_iteration_derivatives_values,
                                                   const std::vector<Eigen::MatrixXd> &basis_functions_values,
                                                   const std::vector<Eigen::MatrixXd> &basis_functions_derivatives_values,
                                                   const Eigen::VectorXd &quadrature_weights) const
{
    MatrixXd cellCMatrix = MatrixXd::Zero(basis_functions_values[0].cols(), basis_functions_values[0].cols());

    for (unsigned int d1 = 0; d1 < 2; d1++)
    {
        for (unsigned int d2 = 0; d2 < 2; d2++)
            cellCMatrix += basis_functions_values[d1].transpose() *
                           quadrature_weights.cwiseProduct(previous_iteration_values[d2]).asDiagonal() *
                           basis_functions_derivatives_values[2 * d1 + d2];
    }

    for (unsigned int d1 = 0; d1 < 2; d1++)
    {
        for (unsigned int d2 = 0; d2 < 2; d2++)
            cellCMatrix += basis_functions_values[d1].transpose() *
                           quadrature_weights.cwiseProduct(previous_iteration_derivatives_values[2 * d1 + d2]).asDiagonal() *
                           basis_functions_values[d2];
    }

    return cellCMatrix;
}
// ***************************************************************************
Eigen::VectorXd Assembler::ComputeConvectiveRightHandSideTerm(const std::vector<Eigen::VectorXd> &previous_iteration_values,
                                                              const std::vector<Eigen::VectorXd> &previous_iteration_derivatives_values,
                                                              const std::vector<Eigen::MatrixXd> &basis_functions_values,
                                                              const Eigen::VectorXd &quadrature_weights) const
{
    VectorXd cellRightHandSide = VectorXd::Zero(basis_functions_values[0].cols());

    for (unsigned int d1 = 0; d1 < 2; d1++)
    {
        for (unsigned int d2 = 0; d2 < 2; d2++)
            cellRightHandSide += basis_functions_values[d1].transpose() * quadrature_weights.asDiagonal() *
                                 previous_iteration_values[d2].cwiseProduct(previous_iteration_derivatives_values[2 * d1 + d2]);
    }

    return cellRightHandSide;
}
// ***************************************************************************
Eigen::MatrixXd Assembler::ComputeSkewMatrix(const std::vector<Eigen::VectorXd> &previous_iteration_values,
                                             const std::vector<Eigen::VectorXd> &previous_iteration_derivatives_values,
                                             const std::vector<Eigen::MatrixXd> &basis_functions_values,
                                             const std::vector<Eigen::MatrixXd> &basis_functions_derivatives_values,
                                             const Eigen::VectorXd &quadrature_weights) const
{
    MatrixXd cellCMatrix = MatrixXd::Zero(basis_functions_values[0].cols(), basis_functions_values[0].cols());

    for (unsigned int d1 = 0; d1 < 2; d1++)
    {
        for (unsigned int d2 = 0; d2 < 2; d2++)
            cellCMatrix += basis_functions_derivatives_values[2 * d1 + d2].transpose() *
                           quadrature_weights.cwiseProduct(previous_iteration_values[d2]).asDiagonal() *
                           basis_functions_values[d1];
    }

    for (unsigned int d1 = 0; d1 < 2; d1++)
    {
        for (unsigned int d2 = 0; d2 < 2; d2++)
            cellCMatrix += basis_functions_derivatives_values[2 * d1 + d2].transpose() *
                           quadrature_weights.cwiseProduct(previous_iteration_values[d1]).asDiagonal() *
                           basis_functions_values[d2];
    }

    return cellCMatrix;
}
// ***************************************************************************
Eigen::VectorXd Assembler::ComputeSkewRightHandSideTerm(const std::vector<Eigen::VectorXd> &previous_iteration_values,
                                                        const std::vector<Eigen::VectorXd> &previous_iteration_derivatives_values,
                                                        const std::vector<Eigen::MatrixXd> &basis_functions_derivatives_values,
                                                        const Eigen::VectorXd &quadrature_weights) const
{
    VectorXd cellRightHandSide = VectorXd::Zero(basis_functions_derivatives_values[0].cols());

    for (unsigned int d1 = 0; d1 < 2; d1++)
    {
        for (unsigned int d2 = 0; d2 < 2; d2++)
            cellRightHandSide += basis_functions_derivatives_values[2 * d1 + d2].transpose() * quadrature_weights.asDiagonal() *
                                 previous_iteration_values[d2].cwiseProduct(previous_iteration_values[d1]);
    }

    return cellRightHandSide;
}
// ***************************************************************************
void Assembler::AssembleNavierStokes(const Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration &config,
                                     const Gedim::MeshMatricesDAO &mesh,
                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                     const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                     const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                     NavierStokes_DF_PCC_2D_Problem_Data &result)
{

    result.globalMatrixC.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    result.dirichletMatrixC.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    result.rightHandSideC.SetSize(count_dofs.num_total_dofs);

    // Assemble equation elements
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

        const auto velocity_basis_functions_derivatives_values =
            Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocityBasisFunctionsDerivativeValues(
                reference_element_data,
                local_space_data,
                Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[3].CellsGlobalDOFs[2].at(c).size();
        const unsigned int num_local_dofs_velocity = local_count_dofs.num_total_dofs - num_local_dofs_pressure;
        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                result.previousIteration,
                                                                                result.solutionDirichlet);

        const Eigen::VectorXd velocity_dofs_values = dofs_values.segment(0, num_local_dofs_velocity);

        vector<VectorXd> previous_iteration_values(2);
        for (unsigned int d = 0; d < 2; d++)
            previous_iteration_values[d] = velocity_basis_functions_values[d] * velocity_dofs_values;

        vector<VectorXd> previous_iteration_derivatives_values(4);
        for (unsigned int d = 0; d < 4; d++)
            previous_iteration_derivatives_values[d] = velocity_basis_functions_derivatives_values[d] * velocity_dofs_values;

        // Compute full cell matrix
        MatrixXd cellMatrixC = MatrixXd::Zero(num_local_dofs_velocity, num_local_dofs_velocity);
        VectorXd cellRightHandSideC = VectorXd::Zero(num_local_dofs_velocity);

        switch (config.ConvectiveForm())
        {
        case Program_configuration::ConvectiveFormType::None:
            break;
        case Program_configuration::ConvectiveFormType::Conv: {

            cellMatrixC = ComputeConvectiveMatrix(previous_iteration_values,
                                                  previous_iteration_derivatives_values,
                                                  velocity_basis_functions_values,
                                                  velocity_basis_functions_derivatives_values,
                                                  internal_quadrature.Weights);

            cellRightHandSideC = ComputeConvectiveRightHandSideTerm(previous_iteration_values,
                                                                    previous_iteration_derivatives_values,
                                                                    velocity_basis_functions_values,
                                                                    internal_quadrature.Weights);
        }
        break;
        case Program_configuration::ConvectiveFormType::Skew: {

            const MatrixXd cellCMatrix1 = ComputeConvectiveMatrix(previous_iteration_values,
                                                                  previous_iteration_derivatives_values,
                                                                  velocity_basis_functions_values,
                                                                  velocity_basis_functions_derivatives_values,
                                                                  internal_quadrature.Weights);

            const MatrixXd cellCMatrix2 = ComputeSkewMatrix(previous_iteration_values,
                                                            previous_iteration_derivatives_values,
                                                            velocity_basis_functions_values,
                                                            velocity_basis_functions_derivatives_values,
                                                            internal_quadrature.Weights);

            cellMatrixC = 0.5 * (cellCMatrix1 - cellCMatrix2);

            const VectorXd cellRightHandSide1 = ComputeConvectiveRightHandSideTerm(previous_iteration_values,
                                                                                   previous_iteration_derivatives_values,
                                                                                   velocity_basis_functions_values,
                                                                                   internal_quadrature.Weights);
            const VectorXd cellRightHandSide2 = ComputeSkewRightHandSideTerm(previous_iteration_values,
                                                                             previous_iteration_derivatives_values,
                                                                             velocity_basis_functions_derivatives_values,
                                                                             internal_quadrature.Weights);

            cellRightHandSideC = 0.5 * (cellRightHandSide1 - cellRightHandSide2);
        }
        break;
        default:
            throw runtime_error("Not valid convective form type");
        }

        assert(Polydim::PDETools::LocalSpace_DF_PCC_2D::VelocitySize(reference_element_data, local_space_data) == num_local_dofs_velocity);

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data = {
            {std::cref(dofs_data[0]), std::cref(dofs_data[1]), std::cref(dofs_data[2])},
            local_count_dofs.offsets_DOFs,
            count_dofs.offsets_DOFs,
            count_dofs.offsets_Strongs};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          cellMatrixC,
                                                                                          cellRightHandSideC,
                                                                                          result.globalMatrixC,
                                                                                          result.dirichletMatrixC,
                                                                                          result.rightHandSideC);
    }

    result.rightHandSideC.Create();
    result.globalMatrixC.Create();
    result.dirichletMatrixC.Create();

    if (count_dofs.num_total_strong > 0)
        result.rightHandSideC.SubtractionMultiplication(result.dirichletMatrixC, result.solutionDirichlet);

    result.rightHandSideC += result.rightHandSide;
    result.globalMatrixC += result.globalMatrixA;
}
// ***************************************************************************
Assembler::Performance_Data Assembler::ComputeMethodPerformance(
    const Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration &config,
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
Assembler::PostProcess_Data Assembler::PostProcessSolution(
    const Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const NavierStokes_DF_PCC_2D_Problem_Data &assembler_data,
    const double &residual_norm,
    const Polydim::examples::NavierStokes_DF_PCC_2D::test::I_Test &test) const
{
    PostProcess_Data result;

    result.residual_norm = residual_norm;

    for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
    {
        result.cell0Ds_numeric_velocity[d].resize(mesh.Cell0DTotalNumber());
        result.cell0Ds_exact_velocity[d].resize(mesh.Cell0DTotalNumber());
    }

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        const auto exact_velocity = test.exact_velocity(mesh.Cell0DCoordinates(p));

        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
            result.cell0Ds_exact_velocity[d](p) = exact_velocity[d](0);

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

    result.cell2Ds_error_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_H1_velocity.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_H1_velocity.setZero(mesh.Cell2DTotalNumber());
    result.error_L2_pressure = 0.0;
    result.norm_L2_pressure = 0.0;
    result.error_H1_velocity = 0.0;
    result.norm_H1_velocity = 0.0;
    result.mesh_size = 0.0;

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

        const Eigen::VectorXd velocity_dofs_values =
            dofs_values.segment(0, local_count_dofs.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd pressure_dofs_values =
            dofs_values.segment(local_count_dofs.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        const VectorXd numeric_pressure_values = pressure_basis_functions_values * pressure_dofs_values;
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
} // namespace NavierStokes_DF_PCC_2D
} // namespace examples
} // namespace Polydim
