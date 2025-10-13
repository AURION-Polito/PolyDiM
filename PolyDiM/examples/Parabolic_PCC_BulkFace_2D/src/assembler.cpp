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
#include "FEM_PCC_1D_Creator.hpp"
#include "VEM_PCC_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_BulkFace_2D
{
// ***************************************************************************
void Assembler::ComputeInitialCondition(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                        const Gedim::MeshMatricesDAO &mesh_2D,
                                        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
                                        const Gedim::MeshMatricesDAO &mesh_1D,
                                        const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
                                        const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                        const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                        const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
                                        const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
                                        const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                                        Gedim::Eigen_Array<> &initial_condition) const
{
    initial_condition.SetSize(count_dofs.num_total_dofs);

    ComputeInitialCondition_2D(config, mesh_2D, mesh_geometric_data_2D, dofs_data[0], count_dofs, reference_element_data_2D, test, initial_condition);

    ComputeInitialCondition_1D(config, mesh_1D, mesh_geometric_data_1D, dofs_data[1], count_dofs, reference_element_data_1D, test, initial_condition);

    initial_condition.Create();
}
// ***************************************************************************
void Assembler::AssembleMatrix(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                               const Gedim::MeshMatricesDAO &mesh_2D,
                               const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
                               const Gedim::MeshMatricesDAO &mesh_1D,
                               const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
                               const Gedim::MeshUtilities::ExtractMeshData &extract_data,
                               const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                               const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                               const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                               const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
                               const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
                               const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                               Gedim::Eigen_SparseArray<> &globalMatrixA,
                               Gedim::Eigen_SparseArray<> &globalMatrixM) const
{
    globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    globalMatrixM.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);

    AssembleMatrix_2D(config, mesh_2D, mesh_geometric_data_2D, mesh_dofs_info[0], dofs_data[0], count_dofs, reference_element_data_2D, test, globalMatrixA, globalMatrixM);

    AssembleMatrix_1D(config, mesh_1D, mesh_geometric_data_1D, mesh_dofs_info[1], dofs_data[1], count_dofs, reference_element_data_1D, test, globalMatrixA, globalMatrixM);

    ComputeTransitionMatrices(config, mesh_1D, mesh_geometric_data_1D, extract_data, mesh_dofs_info, dofs_data, count_dofs, reference_element_data_1D, test, globalMatrixA);

    globalMatrixA.Create();
    globalMatrixM.Create();
}
// ***************************************************************************
void Assembler::AssembleRhs(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                            const double &value_time,
                            const Gedim::MeshMatricesDAO &mesh_2D,
                            const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
                            const Gedim::MeshMatricesDAO &mesh_1D,
                            const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
                            const Gedim::MeshUtilities::ExtractMeshData &extract_data,
                            const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                            const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                            const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                            const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
                            const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
                            const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                            Gedim::Eigen_Array<> &rightHandSide) const
{
    rightHandSide.SetSize(count_dofs.num_total_dofs);

    AssembleRhs_2D(config, value_time, mesh_2D, mesh_geometric_data_2D, mesh_dofs_info[0], dofs_data[0], count_dofs, reference_element_data_2D, test, rightHandSide);

    AssembleRhs_1D(config, value_time, mesh_1D, mesh_geometric_data_1D, mesh_dofs_info[1], dofs_data[1], count_dofs, reference_element_data_1D, test, rightHandSide);

    rightHandSide.Create();
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(
    const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
    const double &value_time,
    const Gedim::MeshMatricesDAO &mesh_2D,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data_2D,
    const Gedim::MeshMatricesDAO &mesh_1D,
    const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data_2D,
    const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
    const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
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
                           value_time,
                           mesh_2D,
                           mesh_geometric_data_2D,
                           dofs_data[0],
                           count_dofs,
                           reference_element_data_2D,
                           assembler_data,
                           test,
                           post_process_data.post_process_data_2D);

    PostProcessSolution_1D(config,
                           value_time,
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
void Assembler::AssembleMatrix_2D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                  const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                  const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                                  Gedim::Eigen_SparseArray<> &globalMatrixA,
                                  Gedim::Eigen_SparseArray<> &globalMatrixM) const
{

    const unsigned int dimension = 2;
    const unsigned int id_domain = 0;

    Gedim::Eigen_SparseArray<> dirichletMatrixA;

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

        const Eigen::MatrixXd local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                            basis_functions_derivative_values,
                                                                            cell2D_internal_quadrature.Weights);

        Eigen::MatrixXd local_C =
            equation.ComputeCellReactionMatrix(reaction_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        Eigen::MatrixXd local_M =
            equation.ComputeCellReactionMatrix(Eigen::VectorXd::Ones(cell2D_internal_quadrature.Weights.size()),
                                               basis_functions_values,
                                               cell2D_internal_quadrature.Weights);

        const double k_max = diffusion_term_values.cwiseAbs().maxCoeff();
        const double g_max = reaction_term_values.cwiseAbs().maxCoeff();

        const Eigen::MatrixXd local_A_stab =
            k_max * Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data, local_space_data);

        const Eigen::MatrixXd local_C_stab =
            Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data,
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
            local_A + local_A_stab + local_C + g_max * local_C_stab,
            globalMatrixA,
            dirichletMatrixA);

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_M + local_C_stab,
                                                                                          globalMatrixM,
                                                                                          dirichletMatrixA);
    }
}
// ***************************************************************************
void Assembler::AssembleRhs_2D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                               const double &value_time,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                               const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                               const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                               const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                               const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                               const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                               Gedim::Eigen_Array<> &rightHandSide) const
{

    const unsigned int dimension = 2;
    const unsigned int id_domain = 0;

    Gedim::Eigen_SparseArray<> dirichletMatrixA;

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
        const auto cell2D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const auto source_term_values = test.source_term(dimension, id_domain, cell2D_internal_quadrature.Points, value_time);

        Eigen::VectorXd local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

        assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {0}, {0}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_rhs,
                                                                                          rightHandSide);
    }
}
// ***************************************************************************
void Assembler::AssembleMatrix_1D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                  const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                  const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                                  Gedim::Eigen_SparseArray<> &globalMatrixA,
                                  Gedim::Eigen_SparseArray<> &globalMatrixM) const
{
    const unsigned int dimension = 1;

    Gedim::Eigen_SparseArray<> dirichletMatrixA;
    const auto local_space =
        Polydim::FEM::PCC::create_FEM_PCC_1D_local_space(FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

    Eigen::Vector3d reference_origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d reference_tangent = Eigen::Vector3d::Zero();
    reference_tangent(0) = 1.0;

    const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(), reference_origin, reference_tangent, 1.0};

    // compute vem values
    const Gedim::Quadrature::QuadratureData reference_quadrature =
        Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);
    const Eigen::VectorXd points_curvilinear_coordinates = reference_quadrature.Points.row(0);
    const unsigned int num_ref_quadrature_points = reference_quadrature.Points.cols();

    const auto local_space_data = local_space->CreateLocalSpace(reference_element_data, segment);

    const auto basis_functions_values =
        local_space->ComputeBasisFunctionsValues(reference_element_data, local_space_data, reference_quadrature.Points);

    const auto basis_functions_derivative_values =
        local_space->ComputeBasisFunctionsDerivativeValues(reference_element_data,
                                                           local_space_data,
                                                           reference_quadrature.Points);

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {

        // map edge internal quadrature points
        const Eigen::Vector3d &edgeStart = mesh_geometric_data.Cell1DsVertices.at(c).col(0);
        const Eigen::Vector3d &edgeTangent = mesh_geometric_data.Cell1DsTangents.at(c);

        Eigen::MatrixXd quadrature_points_2D(3, num_ref_quadrature_points);
        for (unsigned int q = 0; q < num_ref_quadrature_points; q++)
            quadrature_points_2D.col(q) = edgeStart + reference_quadrature.Points(0, q) * edgeTangent;

        const Eigen::VectorXd quadrature_weights_2D = reference_quadrature.Weights * mesh_geometric_data.Cell1DsLengths.at(c);

        const auto diffusion_term_values = test.diffusion_term(dimension, mesh.Cell1DMarker(c), quadrature_points_2D);
        const auto reaction_term_values = test.reaction_term(dimension, mesh.Cell1DMarker(c), quadrature_points_2D);
        const auto beta_term_values = test.beta(dimension, mesh.Cell1DMarker(c), quadrature_points_2D);

        const auto local_A =
            equation.ComputeCellDiffusionMatrix(diffusion_term_values, basis_functions_derivative_values, quadrature_weights_2D);

        const auto local_C =
            equation.ComputeCellReactionMatrix(reaction_term_values + beta_term_values, basis_functions_values, quadrature_weights_2D);

        const auto local_M = equation.ComputeCellReactionMatrix(Eigen::VectorXd::Ones(num_ref_quadrature_points),
                                                                basis_functions_values,
                                                                quadrature_weights_2D);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[1].at(c);

        assert(local_space_data.NumberOfBasisFunctions == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {count_dofs.offsets_DOFs[1]}, {count_dofs.offsets_Strongs[1]}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<1>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A + local_C,
                                                                                          globalMatrixA,
                                                                                          dirichletMatrixA);

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<1>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_M,
                                                                                          globalMatrixM,
                                                                                          dirichletMatrixA);
    }
}
// ***************************************************************************
void Assembler::AssembleRhs_1D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                               const double &value_time,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                               const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                               const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                               const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                               const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                               const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                               Gedim::Eigen_Array<> &rightHandSide) const
{
    const unsigned int dimension = 1;

    Gedim::Eigen_SparseArray<> dirichletMatrixA;
    const auto local_space =
        Polydim::FEM::PCC::create_FEM_PCC_1D_local_space(FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

    Eigen::Vector3d reference_origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d reference_tangent = Eigen::Vector3d::Zero();
    reference_tangent(0) = 1.0;

    const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(), reference_origin, reference_tangent, 1.0};

    // compute vem values
    const Gedim::Quadrature::QuadratureData reference_quadrature =
        Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);
    const Eigen::VectorXd points_curvilinear_coordinates = reference_quadrature.Points.row(0);
    const unsigned int num_ref_quadrature_points = reference_quadrature.Points.cols();

    const auto local_space_data = local_space->CreateLocalSpace(reference_element_data, segment);

    const auto basis_functions_values =
        local_space->ComputeBasisFunctionsValues(reference_element_data, local_space_data, reference_quadrature.Points);

    const auto basis_functions_derivative_values =
        local_space->ComputeBasisFunctionsDerivativeValues(reference_element_data,
                                                           local_space_data,
                                                           reference_quadrature.Points);

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {

        // map edge internal quadrature points
        const Eigen::Vector3d &edgeStart = mesh_geometric_data.Cell1DsVertices.at(c).col(0);
        const Eigen::Vector3d &edgeTangent = mesh_geometric_data.Cell1DsTangents.at(c);

        Eigen::MatrixXd quadrature_points_2D(3, num_ref_quadrature_points);
        for (unsigned int q = 0; q < num_ref_quadrature_points; q++)
            quadrature_points_2D.col(q) = edgeStart + reference_quadrature.Points(0, q) * edgeTangent;

        const Eigen::VectorXd quadrature_weights_2D = reference_quadrature.Weights * mesh_geometric_data.Cell1DsLengths.at(c);

        const auto source_term_values = test.source_term(dimension, mesh.Cell1DMarker(c), quadrature_points_2D, value_time);

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, quadrature_weights_2D);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[1].at(c);

        assert(local_space_data.NumberOfBasisFunctions == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {count_dofs.offsets_DOFs[1]}, {count_dofs.offsets_Strongs[1]}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<1>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_rhs,
                                                                                          rightHandSide);
    }
}
// ***************************************************************************
void Assembler::ComputeTransitionMatrices(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                          const Gedim::MeshMatricesDAO &mesh_1D,
                                          const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data_1D,
                                          const Gedim::MeshUtilities::ExtractMeshData &extract_data,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                          const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                          const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data_1D,
                                          const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                                          Gedim::Eigen_SparseArray<> &globalMatrixA) const
{
    const unsigned int dimension = 1;

    Gedim::Eigen_SparseArray<> dirichletMatrixA;
    const auto local_space =
        Polydim::FEM::PCC::create_FEM_PCC_1D_local_space(FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

    Eigen::Vector3d reference_origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d reference_tangent = Eigen::Vector3d::Zero();
    reference_tangent(0) = 1.0;

    const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(), reference_origin, reference_tangent, 1.0};

    // compute vem values
    const Gedim::Quadrature::QuadratureData reference_quadrature =
        Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data_1D.Order);
    const Eigen::VectorXd points_curvilinear_coordinates = reference_quadrature.Points.row(0);
    const unsigned int num_ref_quadrature_points = reference_quadrature.Points.cols();

    const auto local_space_data = local_space->CreateLocalSpace(reference_element_data_1D, segment);

    const auto basis_functions_values =
        local_space->ComputeBasisFunctionsValues(reference_element_data_1D, local_space_data, reference_quadrature.Points);

    for (unsigned int c = 0; c < mesh_1D.Cell1DTotalNumber(); ++c)
    {

        // map edge internal quadrature points
        const Eigen::Vector3d &edgeStart = mesh_geometric_data_1D.Cell1DsVertices.at(c).col(0);
        const Eigen::Vector3d &edgeTangent = mesh_geometric_data_1D.Cell1DsTangents.at(c);

        Eigen::MatrixXd quadrature_points_2D(3, num_ref_quadrature_points);
        for (unsigned int q = 0; q < num_ref_quadrature_points; q++)
            quadrature_points_2D.col(q) = edgeStart + reference_quadrature.Points(0, q) * edgeTangent;

        const Eigen::VectorXd quadrature_weights_2D =
            reference_quadrature.Weights * mesh_geometric_data_1D.Cell1DsLengths.at(c);

        const auto alpha_term_values = test.alpha(dimension, mesh_1D.Cell1DMarker(c), quadrature_points_2D);
        const auto beta_term_values = test.beta(dimension, mesh_1D.Cell1DMarker(c), quadrature_points_2D);

        const auto &global_dofs = dofs_data[1].CellsGlobalDOFs[1].at(c);

        assert(local_space_data.NumberOfBasisFunctions == global_dofs.size());

        const auto local_alpha = equation.ComputeCellReactionMatrix(alpha_term_values, basis_functions_values, quadrature_weights_2D);
        const auto local_beta = equation.ComputeCellReactionMatrix(beta_term_values, basis_functions_values, quadrature_weights_2D);

        const unsigned int cell1D_index = extract_data.NewCell1DToOldCell1D[c];

        // Assemble system
        for (unsigned int i = 0; i < global_dofs.size(); i++)
        {
            const auto global_dof_i = dofs_data[1].CellsGlobalDOFs[1].at(c).at(i);
            const auto local_dof_i =
                dofs_data[1].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            unsigned int global_index_i = local_dof_i.Global_Index + count_dofs.offsets_DOFs[1];

            for (unsigned int j = 0; j < global_dofs.size(); j++)
            {
                const auto &global_dof_j = dofs_data[0].CellsGlobalDOFs[1].at(cell1D_index).at(j);
                const auto &local_dof_j =
                    dofs_data[0].CellsDOFs.at(global_dof_j.Dimension).at(global_dof_j.CellIndex).at(global_dof_j.DOFIndex);

                const unsigned int global_index_j = local_dof_j.Global_Index;

                globalMatrixA.Triplet(global_index_j, global_index_i, -local_beta(i, j));
                globalMatrixA.Triplet(global_index_i, global_index_j, -local_alpha(i, j));
            }
        }

        // Assemble system
        for (unsigned int i = 0; i < global_dofs.size(); i++)
        {
            const auto global_dof_i = dofs_data[0].CellsGlobalDOFs[1].at(cell1D_index).at(i);
            const auto local_dof_i =
                dofs_data[0].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            unsigned int global_index_i = local_dof_i.Global_Index;

            for (unsigned int j = 0; j < global_dofs.size(); j++)
            {
                const auto &global_dof_j = dofs_data[0].CellsGlobalDOFs[1].at(cell1D_index).at(j);
                const auto &local_dof_j =
                    dofs_data[0].CellsDOFs.at(global_dof_j.Dimension).at(global_dof_j.CellIndex).at(global_dof_j.DOFIndex);

                const unsigned int global_index_j = local_dof_j.Global_Index;

                globalMatrixA.Triplet(global_index_i, global_index_j, local_alpha(i, j));
            }
        }
    }
}
// ***************************************************************************
void Assembler::ComputeInitialCondition_2D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                           const Gedim::IMeshDAO &mesh,
                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                           const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                           const test::I_Test &test,
                                           Gedim::Eigen_Array<> &initial_condition) const
{

    const unsigned int dimension = 2;
    const unsigned int id_domain = 0;

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        // DOFs: vertices
        const Eigen::MatrixXd coordinates = mesh.Cell2DVerticesCoordinates(c);
        const Eigen::VectorXd dofs_vertices = test.initial_solution(dimension, id_domain, coordinates);

        // Assemble local numerical solution
        unsigned int count = 0;
        for (unsigned int p = 0; p < mesh.Cell2DNumberVertices(c); p++)
        {
            const unsigned int cell0D_index = mesh.Cell2DVertex(c, p);

            const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);
            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); loc_i++)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);
                const int global_i = local_dof_i.Global_Index;

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    initial_condition.SetValue(global_i, dofs_vertices(count++));
                }
                break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }

        // Assemble strong boundary condition on Cell1Ds
        if (reference_element_data.Order > 1)
        {

            const auto local_space_data =
                Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                       config.GeometricTolerance2D(),
                                                                       mesh_geometric_data,
                                                                       c,
                                                                       reference_element_data);

            // Assemble strong boundary condition on Cell1Ds
            for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(c); ++ed)
            {
                const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

                const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

                const auto edge_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);

                const Eigen::VectorXd dofs_edge = test.initial_solution(dimension, id_domain, edge_dofs_coordinates);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                    const auto &local_dof_i = local_dofs.at(loc_i);
                    const int global_i = local_dof_i.Global_Index;

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                        initial_condition.SetValue(global_i, dofs_edge(loc_i));
                    }
                    break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }

            const auto local_dofs = dofs_data.CellsDOFs.at(2).at(c);

            if (local_dofs.size())
            {
                const auto internal_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofsCoordinates(reference_element_data, local_space_data);

                const Eigen::VectorXd initial_values_at_dofs =
                    test.initial_solution(dimension, id_domain, internal_dofs_coordinates.Points);

                const Eigen::VectorXd dofs_internal =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofs(reference_element_data,
                                                                       local_space_data,
                                                                       initial_values_at_dofs,
                                                                       internal_dofs_coordinates);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                    const auto &local_dof_i = local_dofs.at(loc_i);
                    const int global_i = local_dof_i.Global_Index;

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                        initial_condition.SetValue(global_i, dofs_internal(loc_i));
                    }
                    break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }
        }
    }
}
// ***************************************************************************
void Assembler::ComputeInitialCondition_1D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                           const Gedim::IMeshDAO &mesh,
                                           const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                           const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                           const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                           const test::I_Test &test,
                                           Gedim::Eigen_Array<> &initial_condition) const
{

    const unsigned int dimension = 1;

    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
    {
        const auto coordinates = mesh.Cell0DCoordinates(p);

        const Eigen::VectorXd dofs_vertices = test.initial_solution(dimension, mesh.Cell0DMarker(p), coordinates);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                initial_condition.SetValue(local_dof_i.Global_Index + count_dofs.offsets_DOFs[1], dofs_vertices(0));
            }
            break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    // Assemble strong boundary condition on Cell1Ds
    if (reference_element_data.Order > 1)
    {

        const auto local_space =
            Polydim::FEM::PCC::create_FEM_PCC_1D_local_space(FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

        // Assemble equation elements
        for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); c++)
        {
            const auto local_dofs = dofs_data.CellsDOFs.at(1).at(c);

            const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(),
                                                                            mesh_geometric_data.Cell1DsVertices.at(c).col(0),
                                                                            mesh_geometric_data.Cell1DsTangents.at(c),
                                                                            1.0};

            const auto local_space_data = local_space->CreateLocalSpace(reference_element_data, segment);

            const auto edge_dofs_coordinates = local_space->InternalDOFsCoordinates(reference_element_data, local_space_data);

            const Eigen::VectorXd dofs_edge = test.initial_solution(dimension, mesh.Cell1DMarker(c), edge_dofs_coordinates);

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);
                const int global_i = local_dof_i.Global_Index + count_dofs.offsets_DOFs[1];

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    initial_condition.SetValue(global_i, dofs_edge(loc_i));
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
Assembler::Performance_Data_2D Assembler::ComputePerformance_2D(
    const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
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
void Assembler::PostProcessSolution_2D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                       const double &value_time,
                                       const Gedim::MeshMatricesDAO &mesh,
                                       const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                       const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                       const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                       const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                       const Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                                       const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                                       Assembler::PostProcess_Data_2D &result) const
{

    const unsigned int dimension = 2;
    const unsigned int id_domain = 0;

    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        result.cell0Ds_exact[p] = test.exact_solution(dimension, id_domain, mesh.Cell0DCoordinates(p), value_time)[0];

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

        const auto exact_solution_values = test.exact_solution(dimension, id_domain, cell2D_internal_quadrature.Points, value_time);
        const auto exact_derivative_solution_values =
            test.exact_derivative_solution(dimension, id_domain, cell2D_internal_quadrature.Points, value_time);

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
void Assembler::PostProcessSolution_1D(const Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration &config,
                                       const double &value_time,
                                       const Gedim::MeshMatricesDAO &mesh,
                                       const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                       const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                       const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                       const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                       const Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                                       const Polydim::examples::Parabolic_PCC_BulkFace_2D::test::I_Test &test,
                                       Assembler::PostProcess_Data_1D &result) const
{
    const unsigned int dimension = 1;

    const auto local_space =
        Polydim::FEM::PCC::create_FEM_PCC_1D_local_space(FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

    Eigen::Vector3d reference_origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d reference_tangent = Eigen::Vector3d::Zero();
    reference_tangent(0) = 1.0;

    const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(), reference_origin, reference_tangent, 1.0};

    const auto local_space_data = local_space->CreateLocalSpace(reference_element_data, segment);

    // compute vem values
    const auto reference_quadrature =
        Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);
    const Eigen::VectorXd points_curvilinear_coordinates = reference_quadrature.Points.row(0);
    const unsigned int num_ref_quadrature_points = reference_quadrature.Points.cols();

    const auto basis_functions_values =
        local_space->ComputeBasisFunctionsValues(reference_element_data, local_space_data, reference_quadrature.Points);

    const auto basis_functions_derivative_values =
        local_space->ComputeBasisFunctionsDerivativeValues(reference_element_data,
                                                           local_space_data,
                                                           reference_quadrature.Points);

    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        result.cell0Ds_exact[p] = test.exact_solution(dimension, mesh.Cell0DMarker(p), mesh.Cell0DCoordinates(p), value_time)[0];

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
        // map edge internal quadrature points
        const Eigen::Vector3d &edgeStart = mesh_geometric_data.Cell1DsVertices.at(c).col(0);
        const Eigen::Vector3d &edgeTangent = mesh_geometric_data.Cell1DsTangents.at(c);

        Eigen::MatrixXd quadrature_points_2D(3, num_ref_quadrature_points);
        for (unsigned int q = 0; q < num_ref_quadrature_points; q++)
            quadrature_points_2D.col(q) = edgeStart + reference_quadrature.Points(0, q) * edgeTangent;

        const Eigen::VectorXd quadrature_weights_2D = reference_quadrature.Weights * mesh_geometric_data.Cell1DsLengths.at(c);

        const auto exact_solution_values = test.exact_solution(dimension, mesh.Cell1DMarker(c), quadrature_points_2D, value_time);
        const auto exact_derivative_solution_values =
            test.exact_derivative_solution(dimension, mesh.Cell1DMarker(c), quadrature_points_2D, value_time);

        const auto tangent = mesh_geometric_data.Cell1DsTangents[c];
        const Eigen::VectorXd tangent_derivatives_values =
            exact_derivative_solution_values[0] * tangent[0] + exact_derivative_solution_values[1] * tangent[1];

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

        result.cell1Ds_error_L2[c] = quadrature_weights_2D.transpose() * local_error_L2;
        result.cell1Ds_norm_L2[c] = quadrature_weights_2D.transpose() * local_norm_L2;

        const Eigen::VectorXd local_error_H1 =
            (basis_functions_derivative_values[0] * dofs_values - tangent_derivatives_values).array().square();
        const Eigen::VectorXd local_norm_H1 = (basis_functions_derivative_values[0] * dofs_values).array().square();

        result.cell1Ds_error_H1[c] = quadrature_weights_2D.transpose() * local_error_H1;
        result.cell1Ds_norm_H1[c] = quadrature_weights_2D.transpose() * local_norm_H1;

        if (mesh_geometric_data.Cell1DsLengths.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell1DsLengths.at(c);
    }

    result.error_L2 = std::sqrt(result.cell1Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell1Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell1Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell1Ds_norm_H1.sum());
}
} // namespace Parabolic_PCC_BulkFace_2D
} // namespace examples
} // namespace Polydim
