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

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
//***************************************************************************
void Assembler::ComputeStrongTerm(const unsigned int &cell3DIndex,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                  const PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                  const test::I_Test &test,
                                  Elliptic_PCC_3D_Problem_Data &assembler_data) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int p = 0; p < mesh.Cell3DNumberVertices(cell3DIndex); ++p)
    {
        const unsigned int cell0D_index = mesh.Cell3DVertex(cell3DIndex, p);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
            continue;

        const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, coordinates);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
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
    for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(cell3DIndex); ++e)
    {
        const unsigned int cell1D_index = mesh.Cell3DEdge(cell3DIndex, e);

        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
            local_dofs.size() == 0)
            continue;

        const auto edge_dofs_coordinates =
            Polydim::PDETools::LocalSpace_PCC_3D::EdgeDofsCoordinates(reference_element_data, local_space_data, e);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, edge_dofs_coordinates);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    // Assemble strong boundary condition on Cell2Ds
    unsigned int quadraturePointOffset = 0;
    for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
    {
        const unsigned int cell2D_index = mesh.Cell3DFace(cell3DIndex, f);

        const auto local_dofs = dofs_data.CellsDOFs.at(2).at(cell2D_index);
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(cell2D_index);

        const auto face_dofs_coordinates =
            Polydim::PDETools::LocalSpace_PCC_3D::FaceDofsCoordinates(reference_element_data, local_space_data, f, quadraturePointOffset);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
            local_dofs.size() == 0)
            continue;

        const Eigen::VectorXd dirichletValues =
            test.strong_boundary_condition(boundary_info.Marker, face_dofs_coordinates.Points);

        const Eigen::VectorXd strong_boundary_values =
            Polydim::PDETools::LocalSpace_PCC_3D::FaceDofs(reference_element_data, local_space_data, f, dirichletValues, face_dofs_coordinates);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
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
void Assembler::ComputeWeakTerm(const unsigned int cell3DIndex,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                const PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                const Polydim::examples::Elliptic_PCC_3D::test::I_Test &test,
                                Elliptic_PCC_3D_Problem_Data &assembler_data) const
{
    // Assemble strong boundary condition on Cell2Ds
    unsigned int quadraturePointOffset = 0;
    for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
    {
        const unsigned int cell2D_index = mesh.Cell3DFace(cell3DIndex, f);
        const auto face_quadrature =
            Polydim::PDETools::LocalSpace_PCC_3D::FaceQuadrature(reference_element_data, local_space_data, f, quadraturePointOffset);

        const auto &boundary_info_2D = mesh_dofs_info.CellsBoundaryInfo.at(2).at(cell2D_index);
        if (boundary_info_2D.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
            continue;

        const Eigen::VectorXd neumannValues = test.weak_boundary_condition(boundary_info_2D.Marker, face_quadrature.Points);

        const auto face_basis_function_values =
            Polydim::PDETools::LocalSpace_PCC_3D::BasisFunctionsValuesOnFace(f,
                                                                             reference_element_data,
                                                                             local_space_data,
                                                                             face_quadrature.Points);

        const Eigen::VectorXd weak_boundary_values =
            face_basis_function_values.transpose() * face_quadrature.Weights.asDiagonal() * neumannValues;

        unsigned int offsetPi0km1 = 0;
        for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(cell2D_index); v++)
        {
            const unsigned int cell0D_index = mesh.Cell2DVertex(cell2D_index, v);

            const auto local_dofs_0D = dofs_data.CellsDOFs.at(0).at(cell0D_index);
            const unsigned int numCell0DLocals = local_dofs_0D.size();
            for (unsigned int loc_i = 0; loc_i < numCell0DLocals; ++loc_i)
            {
                const auto &local_dof_i = local_dofs_0D.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index, weak_boundary_values[loc_i + offsetPi0km1]);
                }
                break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    continue;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            offsetPi0km1 += numCell0DLocals;
        }

        for (unsigned int e = 0; e < mesh.Cell2DNumberEdges(cell2D_index); e++)
        {
            const unsigned int cell1D_index = mesh.Cell2DEdge(cell2D_index, e);

            const auto local_dofs_1D = dofs_data.CellsDOFs.at(1).at(cell1D_index);
            const unsigned int numCell1DLocals = local_dofs_1D.size();
            for (unsigned int loc_i = 0; loc_i < numCell1DLocals; ++loc_i)
            {
                const auto &local_dof_i = local_dofs_1D.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index, weak_boundary_values[loc_i + offsetPi0km1]);
                }
                break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    continue;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            offsetPi0km1 += numCell1DLocals;
        }

        const auto local_dofs_2D = dofs_data.CellsDOFs.at(2).at(cell2D_index);
        for (unsigned int loc_i = 0; loc_i < local_dofs_2D.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs_2D.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index, weak_boundary_values[loc_i + offsetPi0km1]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }
}
// ***************************************************************************
typename Assembler::Elliptic_PCC_3D_Problem_Data Assembler::Assemble(
    const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
    const Polydim::examples::Elliptic_PCC_3D::test::I_Test &test) const
{
    Elliptic_PCC_3D_Problem_Data result;

    result.globalMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
    result.dirichletMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberStrongs);
    result.rightHandSide.SetSize(dofs_data.NumberDOFs);
    result.solution.SetSize(dofs_data.NumberDOFs);
    result.solutionDirichlet.SetSize(dofs_data.NumberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
        // TODO("Remove comments")
        // std::cout << "C3D " << c << " ";
        // std::cout << "C0D " << mesh.Cell3DVertex(c, 0) << ", ";
        // std::cout << mesh.Cell3DVertex(c, 1) << ", ";
        // std::cout << mesh.Cell3DVertex(c, 2) << ", ";
        // std::cout << mesh.Cell3DVertex(c, 3) << " ";
        // std::cout << "C1D " << mesh.Cell3DEdge(c, 0) << ", ";
        // std::cout << mesh.Cell3DEdge(c, 1) << ", ";
        // std::cout << mesh.Cell3DEdge(c, 2) << ", ";
        // std::cout << mesh.Cell3DEdge(c, 3) << ", ";
        // std::cout << mesh.Cell3DEdge(c, 4) << ", ";
        // std::cout << mesh.Cell3DEdge(c, 5) << " ";
        // std::cout << "C2D " << mesh.Cell3DFace(c, 0) << ", ";
        // std::cout << mesh.Cell3DFace(c, 1) << ", ";
        // std::cout << mesh.Cell3DFace(c, 2) << ", ";
        // std::cout << mesh.Cell3DFace(c, 3) << std::endl;

        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_3D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             config.GeometricTolerance3D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const auto basis_functions_values =
            Polydim::PDETools::LocalSpace_PCC_3D::BasisFunctionsValues(reference_element_data, local_space_data);

        const auto basis_functions_derivative_values =
            Polydim::PDETools::LocalSpace_PCC_3D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto cell3D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_3D::InternalQuadrature(reference_element_data, local_space_data);

        const auto diffusion_term_values = test.diffusion_term(cell3D_internal_quadrature.Points);
        const auto source_term_values = test.source_term(cell3D_internal_quadrature.Points);

        const auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                 basis_functions_derivative_values,
                                                                 cell3D_internal_quadrature.Weights);

        const Eigen::MatrixXd local_stab_A =
            diffusion_term_values.cwiseAbs().maxCoeff() *
            Polydim::PDETools::LocalSpace_PCC_3D::StabilizationMatrix(reference_element_data, local_space_data);

        const auto local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell3D_internal_quadrature.Weights);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[3].at(c);

        assert(Polydim::PDETools::LocalSpace_PCC_3D::Size(reference_element_data, local_space_data) == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {0}, {0}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<3>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A + local_stab_A,
                                                                                          local_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.dirichletMatrixA,
                                                                                          result.rightHandSide);

        ComputeWeakTerm(c, mesh, mesh_dofs_info, dofs_data, reference_element_data, local_space_data, test, result);

        ComputeStrongTerm(c, mesh, mesh_dofs_info, dofs_data, reference_element_data, local_space_data, test, result);
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
Assembler::Performance_Data Assembler::ComputePerformance(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
                                                          const Gedim::MeshMatricesDAO &mesh,
                                                          const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                                          const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data) const
{
    Assembler::Performance_Data result;
    result.Cell3DsPerformance.resize(mesh.Cell3DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_3D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             config.GeometricTolerance3D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        result.Cell3DsPerformance[c] =
            Polydim::PDETools::LocalSpace_PCC_3D::ComputePerformance(reference_element_data, local_space_data);
    }

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Polydim::examples::Elliptic_PCC_3D::Program_configuration &config,
                                                           const Gedim::MeshMatricesDAO &mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                           const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                                           const Elliptic_PCC_3D_Problem_Data &assembler_data,
                                                           const Polydim::examples::Elliptic_PCC_3D::test::I_Test &test) const
{
    PostProcess_Data result;

    result.residual_norm = 0.0;
    if (dofs_data.NumberDOFs > 0)
    {
        Gedim::Eigen_Array<> residual;
        residual.SetSize(dofs_data.NumberDOFs);
        residual.SumMultiplication(assembler_data.globalMatrixA, assembler_data.solution);
        residual -= assembler_data.rightHandSide;

        result.residual_norm = residual.Norm();
    }

    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        result.cell0Ds_exact[p] = test.exact_solution(mesh.Cell0DCoordinates(p))[0];

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                result.cell0Ds_numeric[p] = assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index);
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                result.cell0Ds_numeric[p] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    result.cell3Ds_error_L2.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_norm_L2.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_error_H1.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_norm_H1.setZero(mesh.Cell3DTotalNumber());
    result.error_L2 = 0.0;
    result.norm_L2 = 0.0;
    result.error_H1 = 0.0;
    result.norm_H1 = 0.0;
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_3D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             config.GeometricTolerance3D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const auto basis_functions_values =
            Polydim::PDETools::LocalSpace_PCC_3D::BasisFunctionsValues(reference_element_data,
                                                                       local_space_data,
                                                                       Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto basis_functions_derivative_values =
            Polydim::PDETools::LocalSpace_PCC_3D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto cell3D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_3D::InternalQuadrature(reference_element_data, local_space_data);

        const auto exact_solution_values = test.exact_solution(cell3D_internal_quadrature.Points);
        const auto exact_derivative_solution_values = test.exact_derivative_solution(cell3D_internal_quadrature.Points);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<3>(c, dofs_data);
        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<3>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                {0},
                                                                                {0},
                                                                                assembler_data.solution,
                                                                                assembler_data.solutionDirichlet);

        const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values - exact_solution_values).array().square();
        const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();

        result.cell3Ds_error_L2[c] = cell3D_internal_quadrature.Weights.transpose() * local_error_L2;
        result.cell3Ds_norm_L2[c] = cell3D_internal_quadrature.Weights.transpose() * local_norm_L2;

        const Eigen::VectorXd local_error_H1 =
            (basis_functions_derivative_values[0] * dofs_values - exact_derivative_solution_values[0]).array().square() +
            (basis_functions_derivative_values[1] * dofs_values - exact_derivative_solution_values[1]).array().square() +
            (basis_functions_derivative_values[2] * dofs_values - exact_derivative_solution_values[2]).array().square();
        const Eigen::VectorXd local_norm_H1 = (basis_functions_derivative_values[0] * dofs_values).array().square() +
                                              (basis_functions_derivative_values[1] * dofs_values).array().square() +
                                              (basis_functions_derivative_values[2] * dofs_values).array().square();

        result.cell3Ds_error_H1[c] = cell3D_internal_quadrature.Weights.transpose() * local_error_H1;
        result.cell3Ds_norm_H1[c] = cell3D_internal_quadrature.Weights.transpose() * local_norm_H1;

        if (mesh_geometric_data.Cell3DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell3DsDiameters.at(c);
    }

    result.error_L2 = std::sqrt(result.cell3Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell3Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell3Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell3Ds_norm_H1.sum());

    return result;
}
} // namespace Elliptic_PCC_3D
} // namespace examples
} // namespace Polydim
