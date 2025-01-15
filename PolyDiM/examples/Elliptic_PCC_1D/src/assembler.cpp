#include "assembler.hpp"

#include "Assembler_Utilities.hpp"
#include "EllipticEquation.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_1D
{
//***************************************************************************
void Assembler::ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                  const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test,
                                  Elliptic_PCC_1D_Problem_Data &assembler_data) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(p);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
            continue;

        const auto coordinates = mesh.Cell0DCoordinates(p);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, coordinates);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

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
//***************************************************************************
void Assembler::ComputeWeakTerm(const unsigned int cell1DIndex,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test,
                                Elliptic_PCC_1D_Problem_Data &assembler_data) const
{
    for (unsigned int p = 0; p < 2; ++p)
    {
        const unsigned int cell0D_index = mesh.Cell1DVertex(cell1DIndex, p);
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
            continue;

        const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);
        const double direction = p == 0 ? -1.0 : 1.0;

        const Eigen::VectorXd weak_boundary_values = direction * test.weak_boundary_condition(boundary_info.Marker, coordinates);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

        assert(local_dofs.size() == weak_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                continue;
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index, weak_boundary_values[loc_i]);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }
}
// ***************************************************************************
Assembler::Elliptic_PCC_1D_Problem_Data Assembler::Assemble(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                                                            const Gedim::MeshMatricesDAO &mesh,
                                                            const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                                            const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                                            const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                            const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                            const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace &local_space,
                                                            const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test) const
{
    Elliptic_PCC_1D_Problem_Data result;

    result.globalMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
    result.dirichletMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberStrongs);
    result.rightHandSide.SetSize(dofs_data.NumberDOFs);
    result.solution.SetSize(dofs_data.NumberDOFs);
    result.solutionDirichlet.SetSize(dofs_data.NumberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const Polydim::FEM::PCC::FEM_PCC_1D_Segment_Geometry segment = {config.GeometricTolerance1D(),
                                                                        mesh_geometric_data.Cell1DsVertices.at(c).col(0),
                                                                        mesh_geometric_data.Cell1DsTangents.at(c),
                                                                        mesh_geometric_data.Cell1DsLengths.at(c)};

        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, segment);

        const auto basis_functions_values = local_space.ComputeBasisFunctionsValues(reference_element_data, local_space_data);

        const auto basis_functions_derivative_values =
            local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto diffusion_term_values = test.diffusion_term(local_space_data.InternalQuadrature.Points);
        const auto source_term_values = test.source_term(local_space_data.InternalQuadrature.Points);

        const auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                 basis_functions_derivative_values,
                                                                 local_space_data.InternalQuadrature.Weights);

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               basis_functions_values,
                                                               local_space_data.InternalQuadrature.Weights);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[1].at(c);

        assert(local_space_data.NumberOfBasisFunctions == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {0}, {0}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<1>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A,
                                                                                          local_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.dirichletMatrixA,
                                                                                          result.rightHandSide);

        ComputeWeakTerm(c, mesh, mesh_geometric_data, mesh_dofs_info, dofs_data, reference_element_data, test, result);
    }

    ComputeStrongTerm(mesh, mesh_geometric_data, mesh_dofs_info, dofs_data, reference_element_data, test, result);

    result.rightHandSide.Create();
    result.solutionDirichlet.Create();
    result.globalMatrixA.Create();
    result.dirichletMatrixA.Create();

    if (dofs_data.NumberStrongs > 0)
        result.rightHandSide.SubtractionMultiplication(result.dirichletMatrixA, result.solutionDirichlet);

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                                                           const Gedim::MeshMatricesDAO &mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                           const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data,
                                                           const Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace &local_space,
                                                           const Elliptic_PCC_1D_Problem_Data &assembler_data,
                                                           const Polydim::examples::Elliptic_PCC_1D::test::I_Test &test) const
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

        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, segment);

        const auto basis_functions_values = local_space.ComputeBasisFunctionsValues(reference_element_data, local_space_data);

        const auto basis_functions_derivative_values =
            local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto exact_solution_values = test.exact_solution(local_space_data.InternalQuadrature.Points);
        const auto exact_derivative_solution_values = test.exact_derivative_solution(local_space_data.InternalQuadrature.Points);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<1>(c, {std::cref(dofs_data)});
        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<1>(c,
                                                                                {std::cref(dofs_data)},
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                {0},
                                                                                {0},
                                                                                assembler_data.solution,
                                                                                assembler_data.solutionDirichlet);

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

    return result;
}
} // namespace Elliptic_PCC_1D
} // namespace examples
} // namespace Polydim
