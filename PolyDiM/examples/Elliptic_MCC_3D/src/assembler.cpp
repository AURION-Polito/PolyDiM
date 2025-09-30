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

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_3D
{
// ***************************************************************************
void Assembler::ComputeStrongTerm(const unsigned int &cell3DIndex,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry &polyhedron,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                  const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &velocity_local_space_data_data,
                                  const test::I_Test &test,
                                  Elliptic_MCC_3D_Problem_Data &assembler_data) const
{
    unsigned int offsetQuadraturePoints = 0;
    for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
    {
        const unsigned int cell2D_index = mesh.Cell3DFace(cell3DIndex, f);
        const unsigned int numFaceQuadraturePoints =
            velocity_local_space_data_data.BoundaryQuadrature.FacesQuadrature[f].Weights.size();

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(cell2D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
        {
            offsetQuadraturePoints += numFaceQuadraturePoints;
            continue;
        }

        // compute values of Neumann condition
        const double globalDirection = polyhedron.FacesGlobalNormalDirection[f] ? 1.0 : -1.0;

        const VectorXd neumannValues = test.strong_boundary_condition(
            boundary_info.Marker,
            velocity_local_space_data_data.BoundaryQuadrature.Quadrature.Points.middleCols(offsetQuadraturePoints,
                                                                                           numFaceQuadraturePoints));

        const VectorXd strong_boundary_values =
            (1.0 / polyhedron.FacesMeasure[f]) * globalDirection *
            velocity_local_space_data_data.FacesVanderInternal[f].transpose() *
            velocity_local_space_data_data.BoundaryQuadrature.FacesQuadrature[f].Weights.asDiagonal() * neumannValues;

        const auto local_dofs = dofs_data.CellsDOFs.at(2).at(cell2D_index);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); loc_i++)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                assembler_data.solutionNeumann.SetValue(local_dof_i.Global_Index, strong_boundary_values(loc_i));
            }
            break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }

        offsetQuadraturePoints += numFaceQuadraturePoints;
    }
}
// ***************************************************************************
void Assembler::ComputeWeakTerm(const unsigned int cell3DIndex,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry &polyhedron,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &velocity_local_space_data_data,
                                const test::I_Test &test,
                                Elliptic_MCC_3D_Problem_Data &assembler_data) const
{
    unsigned int offsetQuadraturePoints = 0;
    for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
    {
        const unsigned int cell2D_index = mesh.Cell3DFace(cell3DIndex, f);
        const unsigned int numFaceQuadraturePoints =
            velocity_local_space_data_data.BoundaryQuadrature.FacesQuadrature[f].Weights.size();

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(cell2D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
        {
            offsetQuadraturePoints += numFaceQuadraturePoints;
            continue;
        }

        const double globalDirection = polyhedron.FacesGlobalNormalDirection[f] ? 1.0 : -1.0;

        // compute values of Neumann condition
        const VectorXd dirichletValues = test.weak_boundary_condition(
            boundary_info.Marker,
            velocity_local_space_data_data.BoundaryQuadrature.Quadrature.Points.middleCols(offsetQuadraturePoints,
                                                                                           numFaceQuadraturePoints));

        const VectorXd weak_boundary_values =
            -globalDirection * velocity_local_space_data_data.VanderBasisFunctionValuesOnFace[f].transpose() *
            velocity_local_space_data_data.BoundaryQuadrature.FacesQuadrature[f].Weights.asDiagonal() * dirichletValues;

        const auto local_dofs = dofs_data.CellsDOFs.at(2).at(cell2D_index);

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
        offsetQuadraturePoints += numFaceQuadraturePoints;
    }
}
// ***************************************************************************
Assembler::Elliptic_MCC_3D_Problem_Data Assembler::Assemble(
    const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::MCC::VEM_MCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Polydim::VEM::MCC::I_VEM_MCC_3D_Velocity_LocalSpace &vem_velocity_space,
    const Polydim::VEM::MCC::I_VEM_MCC_3D_Pressure_LocalSpace &vem_pressure_space,
    const Polydim::examples::Elliptic_MCC_3D::test::I_Test &test) const
{
    Elliptic_MCC_3D_Problem_Data result;
    result.globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    result.neumannMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    result.rightHandSide.SetSize(count_dofs.num_total_dofs);
    result.solution.SetSize(count_dofs.num_total_dofs);
    result.solutionNeumann.SetSize(count_dofs.num_total_strong);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry polyhedron = {
            config.GeometricTolerance1D(),
            config.GeometricTolerance2D(),
            config.GeometricTolerance3D(),
            mesh_geometric_data.Cell3DsVertices.at(c),
            mesh_geometric_data.Cell3DsCentroids.at(c),
            mesh_geometric_data.Cell3DsVolumes.at(c),
            mesh_geometric_data.Cell3DsDiameters.at(c),
            mesh_geometric_data.Cell3DsTetrahedronPoints.at(c),
            mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c),
            mesh_geometric_data.Cell3DsFacesTranslations.at(c),
            mesh_geometric_data.Cell3DsFacesNormals.at(c),
            mesh_geometric_data.Cell3DsFacesNormalDirections.at(c),
            mesh_geometric_data.Cell3DsFacesNormalGlobalDirection.at(c),
            mesh_geometric_data.Cell3DsFacesAreas.at(c),
            mesh_geometric_data.Cell3DsFaces2DCentroids.at(c),
            mesh_geometric_data.Cell3DsFacesDiameters.at(c),
            mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c)};

        const auto velocity_local_space_data = vem_velocity_space.CreateLocalSpace(velocity_reference_element_data, polyhedron);
        const auto pressure_local_space_data = vem_pressure_space.CreateLocalSpace(pressure_reference_element_data, polyhedron);

        const auto velocity_basis_functions_values =
            vem_velocity_space.ComputeBasisFunctionsValues(velocity_local_space_data, VEM::MCC::ProjectionTypes::Pi0k);
        const auto velocity_basis_functions_divergence_values =
            vem_velocity_space.ComputeBasisFunctionsDivergenceValues(velocity_local_space_data);
        const auto pressure_basis_functions_values = vem_pressure_space.ComputeBasisFunctionsValues(pressure_local_space_data);

        const auto reaction_term_values = test.reaction_term(velocity_local_space_data.InternalQuadrature.Points);
        const auto advection_term_values = test.mixed_advection_term(velocity_local_space_data.InternalQuadrature.Points);
        const auto diffusion_term_values = test.inverse_diffusion_term(velocity_local_space_data.InternalQuadrature.Points);
        const auto source_term_values = test.source_term(velocity_local_space_data.InternalQuadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                           velocity_basis_functions_values,
                                                           velocity_local_space_data.InternalQuadrature.Weights);

        double kmax = 0.0;
        for (const auto &diffusion_term : diffusion_term_values)
        {
            const double max_k = diffusion_term.cwiseAbs().maxCoeff();
            kmax = kmax < max_k ? max_k : kmax;
        }

        local_A += kmax * vem_velocity_space.ComputeDofiDofiStabilizationMatrix(velocity_local_space_data,
                                                                                VEM::MCC::ProjectionTypes::Pi0k);

        const auto local_M = equation.ComputeCellReactionMatrix(reaction_term_values,
                                                                pressure_basis_functions_values,
                                                                velocity_local_space_data.InternalQuadrature.Weights);

        const auto local_T = equation.ComputeCellAdvectionMatrix(advection_term_values,
                                                                 pressure_basis_functions_values,
                                                                 velocity_basis_functions_values,
                                                                 velocity_local_space_data.InternalQuadrature.Weights);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        velocity_local_space_data.InternalQuadrature.Weights.asDiagonal() *
                                        velocity_basis_functions_divergence_values;

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               pressure_basis_functions_values,
                                                               velocity_local_space_data.InternalQuadrature.Weights);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<3>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[3].at(c).size();

        Eigen::MatrixXd elemental_matrix = MatrixXd::Zero(local_count_dofs.num_total_dofs, local_count_dofs.num_total_dofs);
        Eigen::VectorXd elemental_rhs = VectorXd::Zero(local_count_dofs.num_total_dofs);
        elemental_matrix << local_A, -(local_B + local_T).transpose(), local_B, local_M;
        elemental_rhs << VectorXd::Zero(local_count_dofs.num_total_dofs - num_local_dofs_pressure), local_rhs;

        assert(velocity_local_space_data.NumBasisFunctions == local_count_dofs.num_total_dofs - num_local_dofs_pressure);

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data = {
            {std::cref(dofs_data[0]), std::cref(dofs_data[1])},
            local_count_dofs.offsets_DOFs,
            count_dofs.offsets_DOFs,
            count_dofs.offsets_Strongs};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<3>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          elemental_matrix,
                                                                                          elemental_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.neumannMatrixA,
                                                                                          result.rightHandSide);

        ComputeWeakTerm(c, mesh, polyhedron, mesh_dofs_info[0], dofs_data[0], velocity_reference_element_data, velocity_local_space_data, test, result);

        ComputeStrongTerm(c, mesh, polyhedron, mesh_dofs_info[0], dofs_data[0], velocity_reference_element_data, velocity_local_space_data, test, result);
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
Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(
    const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
    const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::MCC::I_VEM_MCC_3D_Velocity_LocalSpace &vem_velocity_space) const
{
    Assembler::VEM_Performance_Result result;
    result.Cell3DsPerformance.resize(mesh.Cell3DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry polyhedron = {
            config.GeometricTolerance1D(),
            config.GeometricTolerance2D(),
            config.GeometricTolerance3D(),
            mesh_geometric_data.Cell3DsVertices.at(c),
            mesh_geometric_data.Cell3DsCentroids.at(c),
            mesh_geometric_data.Cell3DsVolumes.at(c),
            mesh_geometric_data.Cell3DsDiameters.at(c),
            mesh_geometric_data.Cell3DsTetrahedronPoints.at(c),
            mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c),
            mesh_geometric_data.Cell3DsFacesTranslations.at(c),
            mesh_geometric_data.Cell3DsFacesNormals.at(c),
            mesh_geometric_data.Cell3DsFacesNormalDirections.at(c),
            mesh_geometric_data.Cell3DsFacesNormalGlobalDirection.at(c),
            mesh_geometric_data.Cell3DsFacesAreas.at(c),
            mesh_geometric_data.Cell3DsFaces2DCentroids.at(c),
            mesh_geometric_data.Cell3DsFacesDiameters.at(c),
            mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c)};

        const auto velocity_local_space_data = vem_velocity_space.CreateLocalSpace(velocity_reference_element_data, polyhedron);

        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        result.Cell3DsPerformance[c].Analysis = performanceAnalysis.Compute(Polydim::Utilities::Monomials_3D(),
                                                                            velocity_reference_element_data.MonomialsKp1,
                                                                            vem_velocity_space,
                                                                            velocity_local_space_data);

        result.Cell3DsPerformance[c].NumInternalQuadraturePoints = velocity_local_space_data.InternalQuadrature.Weights.size();
        result.Cell3DsPerformance[c].NumBoundaryQuadraturePoints =
            velocity_local_space_data.BoundaryQuadrature.Quadrature.Weights.size();
    }

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(
    const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
    const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
    const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::MCC::VEM_MCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Polydim::VEM::MCC::I_VEM_MCC_3D_Velocity_LocalSpace &vem_velocity_space,
    const Polydim::VEM::MCC::I_VEM_MCC_3D_Pressure_LocalSpace &vem_pressure_space,
    const Elliptic_MCC_3D_Problem_Data &assembler_data,
    const Polydim::examples::Elliptic_MCC_3D::test::I_Test &test) const
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

    result.cell3Ds_numeric_pressure.resize(mesh.Cell3DTotalNumber());
    result.cell3Ds_exact_pressure.resize(mesh.Cell3DTotalNumber());

    result.cell3Ds_error_L2_pressure.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_super_error_L2_pressure.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_norm_L2_pressure.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_error_L2_velocity.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_norm_L2_velocity.setZero(mesh.Cell3DTotalNumber());
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry polyhedron = {
            config.GeometricTolerance1D(),
            config.GeometricTolerance2D(),
            config.GeometricTolerance3D(),
            mesh_geometric_data.Cell3DsVertices.at(c),
            mesh_geometric_data.Cell3DsCentroids.at(c),
            mesh_geometric_data.Cell3DsVolumes.at(c),
            mesh_geometric_data.Cell3DsDiameters.at(c),
            mesh_geometric_data.Cell3DsTetrahedronPoints.at(c),
            mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c),
            mesh_geometric_data.Cell3DsFacesTranslations.at(c),
            mesh_geometric_data.Cell3DsFacesNormals.at(c),
            mesh_geometric_data.Cell3DsFacesNormalDirections.at(c),
            mesh_geometric_data.Cell3DsFacesNormalGlobalDirection.at(c),
            mesh_geometric_data.Cell3DsFacesAreas.at(c),
            mesh_geometric_data.Cell3DsFaces2DCentroids.at(c),
            mesh_geometric_data.Cell3DsFacesDiameters.at(c),
            mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c)};

        const auto velocity_local_space_data = vem_velocity_space.CreateLocalSpace(velocity_reference_element_data, polyhedron);
        const auto pressure_local_space_data = vem_pressure_space.CreateLocalSpace(pressure_reference_element_data, polyhedron);

        const auto velocity_basis_functions_values =
            vem_velocity_space.ComputeBasisFunctionsValues(velocity_local_space_data, VEM::MCC::ProjectionTypes::Pi0k);
        const auto pressure_basis_functions_values = vem_pressure_space.ComputeBasisFunctionsValues(pressure_local_space_data);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<3>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[3].at(c).size();

        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<3>(c,
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

        const auto exact_pressure_values = test.exact_pressure(velocity_local_space_data.InternalQuadrature.Points);
        const auto exact_velocity_values = test.exact_velocity(velocity_local_space_data.InternalQuadrature.Points);

        result.cell3Ds_numeric_pressure[c] =
            ((1.0 / polyhedron.Measure) * (pressure_basis_functions_values * pressure_dofs_values).transpose() *
             pressure_local_space_data.InternalQuadrature.Weights);

        result.cell3Ds_exact_pressure[c] = ((1.0 / polyhedron.Measure) * (exact_pressure_values).transpose() *
                                            pressure_local_space_data.InternalQuadrature.Weights);

        // Interpolate Exact Solution
        const VectorXd rightHandSide = pressure_local_space_data.VanderInternal.transpose() *
                                       pressure_local_space_data.InternalQuadrature.Weights.asDiagonal() * exact_pressure_values;
        const VectorXd coeffPolynomial = pressure_local_space_data.Hmatrix.llt().solve(rightHandSide);

        const Eigen::VectorXd local_error_L2_pressure =
            (pressure_basis_functions_values * pressure_dofs_values - exact_pressure_values).array().square();
        const Eigen::VectorXd local_super_error_L2_pressure =
            (pressure_basis_functions_values * pressure_dofs_values - pressure_local_space_data.VanderInternal * coeffPolynomial)
                .array()
                .square();
        const Eigen::VectorXd local_norm_L2_pressure = (pressure_basis_functions_values * pressure_dofs_values).array().square();

        result.cell3Ds_super_error_L2_pressure[c] =
            velocity_local_space_data.InternalQuadrature.Weights.transpose() * local_super_error_L2_pressure;
        result.cell3Ds_error_L2_pressure[c] = velocity_local_space_data.InternalQuadrature.Weights.transpose() * local_error_L2_pressure;
        result.cell3Ds_norm_L2_pressure[c] = velocity_local_space_data.InternalQuadrature.Weights.transpose() * local_norm_L2_pressure;

        const Eigen::VectorXd local_error_L2_velocity =
            (velocity_basis_functions_values[0] * velocity_dofs_values - exact_velocity_values[0]).array().square() +
            (velocity_basis_functions_values[1] * velocity_dofs_values - exact_velocity_values[1]).array().square();
        const Eigen::VectorXd local_norm_L2_velocity =
            (velocity_basis_functions_values[0] * velocity_dofs_values).array().square() +
            (velocity_basis_functions_values[1] * velocity_dofs_values).array().square();

        result.cell3Ds_error_L2_velocity[c] = velocity_local_space_data.InternalQuadrature.Weights.transpose() * local_error_L2_velocity;
        result.cell3Ds_norm_L2_velocity[c] = velocity_local_space_data.InternalQuadrature.Weights.transpose() * local_norm_L2_velocity;

        if (mesh_geometric_data.Cell3DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell3DsDiameters.at(c);
    }

    result.error_L2_pressure = std::sqrt(result.cell3Ds_error_L2_pressure.sum());
    result.super_error_L2_pressure = std::sqrt(result.cell3Ds_super_error_L2_pressure.sum());
    result.norm_L2_pressure = std::sqrt(result.cell3Ds_norm_L2_pressure.sum());
    result.error_L2_velocity = std::sqrt(result.cell3Ds_error_L2_velocity.sum());
    result.norm_L2_velocity = std::sqrt(result.cell3Ds_norm_L2_velocity.sum());

    return result;
}
// ***************************************************************************
} // namespace Elliptic_MCC_3D
} // namespace examples
} // namespace Polydim
