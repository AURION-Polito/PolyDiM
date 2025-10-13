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
#include "Eigen_LUSolver.hpp"
#include "EllipticEquation.hpp"
#include "FEM_PCC_1D_Creator.hpp"
#include "VEM_PCC_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_BulkFace_2D
{
// ***************************************************************************
Assembler::Elliptic_PCC_BF_2D_Problem_Data Assembler::Solve(
    const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh_2D,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
    const Gedim::MeshMatricesDAO &mesh_1D,
    const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
    const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
    const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test) const
{
    Assembler::Elliptic_PCC_BF_2D_Problem_Data assembler_data;

    assembler_data.globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    assembler_data.rightHandSide.SetSize(count_dofs.num_total_dofs);

    Gedim::Output::PrintGenericMessage("AssembleSystem Discrete Type " +
                                           std::to_string(static_cast<unsigned int>(config.MethodType())) + "...",
                                       true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Assemble_2D(config, mesh_2D, mesh_geometric_data_2D, mesh_dofs_info[0], dofs_data[0], count_dofs, reference_element_data_2D, test, assembler_data);

    Assemble_1D(config, mesh_1D, mesh_geometric_data_1D, mesh_dofs_info[1], dofs_data[1], count_dofs, reference_element_data_1D, test, assembler_data);

    assembler_data.rightHandSide.Create();
    assembler_data.globalMatrixA.Create();

    Gedim::Profiler::StopTime("AssembleSystem");
    Gedim::Output::PrintStatusProgram("AssembleSystem");

    assembler_data.solution.SetSize(count_dofs.num_total_dofs);

    Gedim::Output::PrintGenericMessage("Factorize...", true);
    Gedim::Profiler::StartTime("Factorize");

    Gedim::Eigen_LUSolver solver;
    solver.Initialize(assembler_data.globalMatrixA);

    Gedim::Profiler::StopTime("Factorize");
    Gedim::Output::PrintStatusProgram("Factorize");

    Gedim::Output::PrintGenericMessage("Solve...", true);
    Gedim::Profiler::StartTime("Solve");

    solver.Solve(assembler_data.rightHandSide, assembler_data.solution);

    Gedim::Profiler::StopTime("Solve");
    Gedim::Output::PrintStatusProgram("Solve");

    return assembler_data;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(
    const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh_2D,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
    const Gedim::MeshMatricesDAO &mesh_1D,
    const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
    const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
    const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
    const Elliptic_PCC_BF_2D_Problem_Data &assembler_data) const
{
    PostProcess_Data post_process_data;

    post_process_data.residual_norm = 0.0;

    Gedim::Eigen_Array<> residual;
    residual.SetSize(count_dofs.num_total_dofs);
    residual.SumMultiplication(assembler_data.globalMatrixA, assembler_data.solution);
    residual -= assembler_data.rightHandSide;

    post_process_data.residual_norm = residual.Norm();

    PostProcessSolution_2D(config,
                           mesh_2D,
                           mesh_geometric_data_2D,
                           dofs_data[0],
                           count_dofs,
                           reference_element_data_2D,
                           assembler_data,
                           test,
                           post_process_data.post_process_data_2D);

    PostProcessSolution_1D(config,
                           mesh_1D,
                           mesh_geometric_data_1D,
                           dofs_data[1],
                           count_dofs,
                           reference_element_data_1D,
                           assembler_data,
                           test,
                           post_process_data.post_process_data_1D);

    return post_process_data;
}
// ***************************************************************************
void Assembler::Assemble_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                            const Gedim::MeshMatricesDAO &mesh,
                            const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                            const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                            const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                            const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                            const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                            const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                            Assembler::Elliptic_PCC_BF_2D_Problem_Data &assembler_data) const
{

    const unsigned int dimension = 2;
    const unsigned int id_domain = 0;

    Gedim::Eigen_SparseArray<> dirichletMatrixA;

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const auto basis_functions_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data,
                                                                       local_space_data,
                                                                       Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto basis_functions_derivative_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const auto diffusion_term_values = test.diffusion_term(dimension, id_domain, cell2D_internal_quadrature.Points);
        const auto reaction_term_values = test.reaction_term(dimension, id_domain, cell2D_internal_quadrature.Points);
        const auto source_term_values = test.source_term(dimension, id_domain, cell2D_internal_quadrature.Points);

        const Eigen::MatrixXd local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                            basis_functions_derivative_values,
                                                                            cell2D_internal_quadrature.Weights);

        Eigen::MatrixXd local_C =
            equation.ComputeCellReactionMatrix(reaction_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        Eigen::VectorXd local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        const double k_max = diffusion_term_values.cwiseAbs().maxCoeff();
        const double g_max = reaction_term_values.cwiseAbs().maxCoeff();

        const Eigen::MatrixXd local_A_stab =
            k_max * Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data, local_space_data);

        const Eigen::MatrixXd local_C_stab =
            g_max * Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data,
                                                                              local_space_data,
                                                                              Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

        assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {0}, {0}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(
            c,
            local_matrix_to_global_matrix_dofs_data,
            local_matrix_to_global_matrix_dofs_data,
            local_A + local_A_stab + local_C + local_C_stab,
            local_rhs,
            assembler_data.globalMatrixA,
            dirichletMatrixA,
            assembler_data.rightHandSide);
    }
}
// ***************************************************************************
void Assembler::Assemble_1D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                            const Gedim::MeshMatricesDAO &mesh,
                            const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                            const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                            const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                            const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                            const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                            const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                            Assembler::Elliptic_PCC_BF_2D_Problem_Data &assembler_data) const
{
    const unsigned int dimension = 1;

    Gedim::Eigen_SparseArray<> dirichletMatrixA;
    const auto local_space =
        Polydim::FEM::PCC::create_FEM_PCC_1D_local_space(FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(),
                                                                        mesh_geometric_data.Cell1DsVertices.at(c).col(0),
                                                                        mesh_geometric_data.Cell1DsTangents.at(c),
                                                                        mesh_geometric_data.Cell1DsLengths.at(c)};

        const auto local_space_data = local_space->CreateLocalSpace(reference_element_data, segment);

        const auto basis_functions_values = local_space->ComputeBasisFunctionsValues(reference_element_data, local_space_data);

        const auto basis_functions_derivative_values =
            local_space->ComputeBasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto diffusion_term_values =
            test.diffusion_term(dimension, mesh.Cell1DMarker(c), local_space_data.InternalQuadrature.Points);
        const auto reaction_term_values =
            test.reaction_term(dimension, mesh.Cell1DMarker(c), local_space_data.InternalQuadrature.Points);
        const auto source_term_values =
            test.source_term(dimension, mesh.Cell1DMarker(c), local_space_data.InternalQuadrature.Points);

        const auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                 basis_functions_derivative_values,
                                                                 local_space_data.InternalQuadrature.Weights);

        const auto local_C = equation.ComputeCellReactionMatrix(reaction_term_values,
                                                                basis_functions_values,
                                                                local_space_data.InternalQuadrature.Weights);

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               basis_functions_values,
                                                               local_space_data.InternalQuadrature.Weights);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[1].at(c);

        assert(local_space_data.NumberOfBasisFunctions == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {count_dofs.offsets_DOFs[1]}, {count_dofs.offsets_Strongs[1]}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<1>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A + local_C,
                                                                                          local_rhs,
                                                                                          assembler_data.globalMatrixA,
                                                                                          dirichletMatrixA,
                                                                                          assembler_data.rightHandSide);
    }
}
// ***************************************************************************
Assembler::Performance_Data_2D Assembler::ComputePerformance_2D(
    const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data) const
{
    Assembler::Performance_Data_2D result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        result.Cell2DsPerformance[c] =
            Polydim::PDETools::LocalSpace_PCC_2D::ComputePerformance(reference_element_data, local_space_data);
    }

    return result;
}
// ***************************************************************************
void Assembler::PostProcessSolution_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                       const Gedim::MeshMatricesDAO &mesh,
                                       const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                       const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                       const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                       const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                       const Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                                       const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                                       Assembler::PostProcess_Data_2D &result) const
{

    const unsigned int dimension = 2;
    const unsigned int id_domain = 0;

    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        result.cell0Ds_exact[p] = test.exact_solution(dimension, id_domain, mesh.Cell0DCoordinates(p))[0];

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                result.cell0Ds_numeric[p] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    result.cell2Ds_error_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_H1.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_H1.setZero(mesh.Cell2DTotalNumber());
    result.error_L2 = 0.0;
    result.norm_L2 = 0.0;
    result.error_H1 = 0.0;
    result.norm_H1 = 0.0;
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const auto basis_functions_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data,
                                                                       local_space_data,
                                                                       Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto basis_functions_derivative_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const auto exact_solution_values = test.exact_solution(dimension, id_domain, cell2D_internal_quadrature.Points);
        const auto exact_derivative_solution_values =
            test.exact_derivative_solution(dimension, id_domain, cell2D_internal_quadrature.Points);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                {0},
                                                                                {0},
                                                                                assembler_data.solution,
                                                                                {});

        const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values - exact_solution_values).array().square();
        const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();

        result.cell2Ds_error_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_L2;
        result.cell2Ds_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_L2;

        const Eigen::VectorXd local_error_H1 =
            (basis_functions_derivative_values[0] * dofs_values - exact_derivative_solution_values[0]).array().square() +
            (basis_functions_derivative_values[1] * dofs_values - exact_derivative_solution_values[1]).array().square();

        const Eigen::VectorXd local_norm_H1 = (basis_functions_derivative_values[0] * dofs_values).array().square() +
                                              (basis_functions_derivative_values[1] * dofs_values).array().square();

        result.cell2Ds_error_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_H1;
        result.cell2Ds_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_H1;

        if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2 = std::sqrt(result.cell2Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell2Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell2Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell2Ds_norm_H1.sum());
}
// ***************************************************************************
void Assembler::PostProcessSolution_1D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                       const Gedim::MeshMatricesDAO &mesh,
                                       const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                       const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                       const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                       const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                       const Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                                       const Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test &test,
                                       Assembler::PostProcess_Data_1D &result) const
{
    const unsigned int dimension = 1;

    const auto local_space =
        Polydim::FEM::PCC::create_FEM_PCC_1D_local_space(FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        result.cell0Ds_exact[p] = test.exact_solution(dimension, mesh.Cell0DMarker(p), mesh.Cell0DCoordinates(p))[0];

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                result.cell0Ds_numeric[p] =
                    assembler_data.solution.GetValue(local_dof_i.Global_Index + count_dofs.offsets_DOFs[1]);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    result.cell1Ds_error_L2.setZero(mesh.Cell1DTotalNumber());
    result.cell1Ds_norm_L2.setZero(mesh.Cell1DTotalNumber());
    result.cell1Ds_error_H1.setZero(mesh.Cell1DTotalNumber());
    result.cell1Ds_norm_H1.setZero(mesh.Cell1DTotalNumber());
    result.error_L2 = 0.0;
    result.norm_L2 = 0.0;
    result.error_H1 = 0.0;
    result.norm_H1 = 0.0;
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); c++)
    {
        const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(),
                                                                        mesh_geometric_data.Cell1DsVertices.at(c).col(0),
                                                                        mesh_geometric_data.Cell1DsTangents.at(c),
                                                                        mesh_geometric_data.Cell1DsLengths.at(c)};

        const auto local_space_data = local_space->CreateLocalSpace(reference_element_data, segment);

        const auto basis_functions_values = local_space->ComputeBasisFunctionsValues(reference_element_data, local_space_data);

        const auto basis_functions_derivative_values =
            local_space->ComputeBasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto exact_solution_values =
            test.exact_solution(dimension, mesh.Cell1DMarker(c), local_space_data.InternalQuadrature.Points);
        const auto exact_derivative_solution_values =
            test.exact_derivative_solution(dimension, mesh.Cell1DMarker(c), local_space_data.InternalQuadrature.Points);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<1>(c, dofs_data);
        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<1>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                {count_dofs.offsets_DOFs[1]},
                                                                                {count_dofs.offsets_Strongs[1]},
                                                                                assembler_data.solution,
                                                                                {});

        const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values - exact_solution_values).array().square();
        const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();

        result.cell1Ds_error_L2[c] = local_space_data.InternalQuadrature.Weights.transpose() * local_error_L2;
        result.cell1Ds_norm_L2[c] = local_space_data.InternalQuadrature.Weights.transpose() * local_norm_L2;

        const Eigen::VectorXd local_error_H1 =
            (basis_functions_derivative_values[0] * dofs_values - exact_derivative_solution_values[0]).array().square();
        const Eigen::VectorXd local_norm_H1 = (basis_functions_derivative_values[0] * dofs_values).array().square();

        result.cell1Ds_error_H1[c] = local_space_data.InternalQuadrature.Weights.transpose() * local_error_H1;
        result.cell1Ds_norm_H1[c] = local_space_data.InternalQuadrature.Weights.transpose() * local_norm_H1;

        if (mesh_geometric_data.Cell1DsLengths.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell1DsLengths.at(c);
    }

    result.error_L2 = std::sqrt(result.cell1Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell1Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell1Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell1Ds_norm_H1.sum());
}
} // namespace Elliptic_PCC_BulkFace_2D
} // namespace examples
} // namespace Polydim
