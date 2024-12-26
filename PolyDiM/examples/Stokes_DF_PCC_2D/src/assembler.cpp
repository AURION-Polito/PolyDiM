#include "assembler.hpp"

#include "Assembler_Utilities.hpp"
#include "EllipticEquation.hpp"
#include "VEM_DF_PCC_2D_Velocity_LocalSpace_Creator.hpp"
#include "program_configuration.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_2D
{
// ***************************************************************************
void Assembler::ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                  const std::vector<size_t> &offsetStrongs,
                                  const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                  const test::I_Test &test,
                                  Stokes_DF_PCC_2D_Problem_Data &assembler_data) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int h = 0; h < 2; h++)
    {
        for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
        {
            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(0).at(p);

            if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
                continue;

            const auto coordinates = mesh.Cell0DCoordinates(p);

            const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, coordinates)[h];

            const auto local_dofs = dofs_data[h].CellsDOFs.at(0).at(p);

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
        const auto &referenceSegmentInternalPoints = reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints;
        const unsigned int numReferenceSegmentInternalPoints = referenceSegmentInternalPoints.cols();

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(e);

            if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
                continue;

            const auto cell1D_origin = mesh.Cell1DOriginCoordinates(e);
            const auto cell1D_end = mesh.Cell1DEndCoordinates(e);
            const auto cell1D_tangent = cell1D_end - cell1D_origin;

            Eigen::MatrixXd coordinates = Eigen::MatrixXd::Zero(3, numReferenceSegmentInternalPoints);
            for (unsigned int r = 0; r < numReferenceSegmentInternalPoints; r++)
                coordinates.col(r) << cell1D_origin + referenceSegmentInternalPoints(0, r) * cell1D_tangent;

            const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, coordinates)[h];

            const auto local_dofs = dofs_data[h].CellsDOFs.at(1).at(e);

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
Assembler::Stokes_DF_PCC_2D_Problem_Data Assembler::Assemble(
    const Polydim::examples::Stokes_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Polydim::examples::Stokes_DF_PCC_2D::test::I_Test &test) const
{
    const unsigned int numDOFHandler = mesh_dofs_info.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<size_t> offsetDOFs = {0,
                                      dofs_data[0].NumberDOFs,
                                      dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs,
                                      dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs + dofs_data[2].NumberDOFs};
    std::vector<size_t> offsetStrongs = {0,
                                         dofs_data[0].NumberStrongs,
                                         dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs,
                                         dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs + dofs_data[2].NumberStrongs};

    for (unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }

    Stokes_DF_PCC_2D_Problem_Data result;
    result.globalMatrixA.SetSize(numberDOFs + 1, numberDOFs + 1, Gedim::ISparseArray::SparseArrayTypes::None);
    result.dirichletMatrixA.SetSize(numberDOFs + 1, numberStrongs);
    result.rightHandSide.SetSize(numberDOFs + 1);
    result.solution.SetSize(numberDOFs + 1);
    result.solutionDirichlet.SetSize(numberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry polygon = {config.GeometricTolerance1D(),
                                                                              config.GeometricTolerance2D(),
                                                                              mesh_geometric_data.Cell2DsVertices.at(c),
                                                                              mesh_geometric_data.Cell2DsCentroids.at(c),
                                                                              mesh_geometric_data.Cell2DsAreas.at(c),
                                                                              mesh_geometric_data.Cell2DsDiameters.at(c),
                                                                              mesh_geometric_data.Cell2DsTriangulations.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeLengths.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeDirections.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeTangents.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeNormals.at(c)};

        const auto vem_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_local_space(config.VemType());
        const auto local_space = vem_local_space->CreateLocalSpace(velocity_reference_element_data, polygon);

        const auto velocity_basis_functions_values =
            vem_local_space->ComputeBasisFunctionsValues(local_space, Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

        const auto velocity_basis_functions_derivatives_values =
            vem_local_space->ComputeBasisFunctionsDerivativeValues(local_space, Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);
        const auto velocity_basis_functions_divergence_values = vem_local_space->ComputeBasisFunctionsDivergenceValues(local_space);

        const MatrixXd pressure_basis_functions_values =
            vem_local_space->ComputePolynomialsValues(local_space).leftCols(local_space.Nkm1);

        const auto diffusion_term_values = test.diffusion_term(local_space.InternalQuadrature.Points);
        const auto source_term_values = test.source_term(local_space.InternalQuadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                           velocity_basis_functions_derivatives_values,
                                                           local_space.InternalQuadrature.Weights);

        double kmax = diffusion_term_values.cwiseAbs().maxCoeff();
        local_A += kmax * local_space.StabMatrix;

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        local_space.InternalQuadrature.Weights.asDiagonal() * velocity_basis_functions_divergence_values;

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               velocity_basis_functions_values,
                                                               local_space.InternalQuadrature.Weights);

        const std::vector<size_t> offset_global_dofs = {
            0lu,
            dofs_data[0].CellsGlobalDOFs[2].at(c).size(),
            dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size(),
            dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                dofs_data[2].CellsGlobalDOFs[2].at(c).size()};
        const std::vector<std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData::GlobalCell_DOF>> global_dofs = {
            dofs_data[0].CellsGlobalDOFs[2].at(c),
            dofs_data[1].CellsGlobalDOFs[2].at(c),
            dofs_data[2].CellsGlobalDOFs[2].at(c),
            dofs_data[3].CellsGlobalDOFs[2].at(c)};

        Eigen::MatrixXd elemental_matrix =
            MatrixXd::Zero(global_dofs[0].size() + global_dofs[1].size() + global_dofs[2].size() + global_dofs[3].size(),
                           global_dofs[0].size() + global_dofs[1].size() + global_dofs[2].size() + global_dofs[3].size());
        Eigen::VectorXd elemental_rhs =
            VectorXd::Zero(global_dofs[0].size() + global_dofs[1].size() + global_dofs[2].size() + global_dofs[3].size());
        elemental_matrix << local_A, local_B.transpose(), local_B, MatrixXd::Zero(global_dofs[3].size(), global_dofs[3].size());

        elemental_rhs << local_rhs, VectorXd::Zero(global_dofs[3].size());

        assert(local_space.NumBasisFunctions == global_dofs[0].size() + global_dofs[1].size() + global_dofs[2].size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data = {
            {std::cref(dofs_data[0]), std::cref(dofs_data[1]), std::cref(dofs_data[2]), std::cref(dofs_data[3])},
            {0lu,
             dofs_data[0].CellsGlobalDOFs[2].at(c).size(),
             dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size(),
             dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                 dofs_data[2].CellsGlobalDOFs[2].at(c).size()},
            offsetDOFs,
            offsetStrongs};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          elemental_matrix,
                                                                                          elemental_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.dirichletMatrixA,
                                                                                          result.rightHandSide);

        // Compute mean values
        const VectorXd mean_value_pressure =
            pressure_basis_functions_values.transpose() * local_space.InternalQuadrature.Weights;

        // Mean value condition
        const unsigned int h1 = 3;
        const unsigned int num_global_offset_lagrange =
            dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs + dofs_data[2].NumberDOFs + dofs_data[3].NumberDOFs;
        for (unsigned int loc_i = 0; loc_i < global_dofs[h1].size(); loc_i++)
        {
            const auto global_dof_i = global_dofs[h1].at(loc_i);
            const auto local_dof_i =
                dofs_data[h1].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);
            const unsigned int global_index_i = local_dof_i.Global_Index + offsetDOFs[h1];

            result.globalMatrixA.Triplet(global_index_i, num_global_offset_lagrange, mean_value_pressure(loc_i));

            result.globalMatrixA.Triplet(num_global_offset_lagrange, global_index_i, mean_value_pressure(loc_i));
        }
    }

    ComputeStrongTerm(mesh, mesh_geometric_data, mesh_dofs_info, dofs_data, offsetStrongs, velocity_reference_element_data, test, result);

    result.rightHandSide.Create();
    result.solutionDirichlet.Create();
    result.globalMatrixA.Create();
    result.dirichletMatrixA.Create();

    if (numberStrongs > 0)
        result.rightHandSide.SubtractionMultiplication(result.dirichletMatrixA, result.solutionDirichlet);

    return result;
}

// ***************************************************************************
Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(
    const Polydim::examples::Stokes_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data) const
{
    Assembler::VEM_Performance_Result result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry polygon = {config.GeometricTolerance1D(),
                                                                              config.GeometricTolerance2D(),
                                                                              mesh_geometric_data.Cell2DsVertices.at(c),
                                                                              mesh_geometric_data.Cell2DsCentroids.at(c),
                                                                              mesh_geometric_data.Cell2DsAreas.at(c),
                                                                              mesh_geometric_data.Cell2DsDiameters.at(c),
                                                                              mesh_geometric_data.Cell2DsTriangulations.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeLengths.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeDirections.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeTangents.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeNormals.at(c)};

        const auto vem_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_local_space(config.VemType());

        const auto local_space = vem_local_space->CreateLocalSpace(velocity_reference_element_data, polygon);

        Polydim::VEM::DF_PCC::VEM_DF_PCC_PerformanceAnalysis performanceAnalysis;

        result.Cell2DsPerformance[c].Analysis = performanceAnalysis.Compute(polygon.Measure,
                                                                            polygon.Diameter,
                                                                            Polydim::VEM::Monomials::VEM_Monomials_2D(),
                                                                            velocity_reference_element_data.Monomials,
                                                                            vem_local_space,
                                                                            local_space);

        result.Cell2DsPerformance[c].NumInternalQuadraturePoints = local_space.InternalQuadrature.Weights.size();
        result.Cell2DsPerformance[c].NumBoundaryQuadraturePoints = local_space.BoundaryQuadrature.Quadrature.Weights.size();
    }

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(
    const Polydim::examples::Stokes_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Stokes_DF_PCC_2D_Problem_Data &assembler_data,
    const Polydim::examples::Stokes_DF_PCC_2D::test::I_Test &test) const
{
    const unsigned int numDOFHandler = dofs_data.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<unsigned int> offsetDOFs = {0,
                                            dofs_data[0].NumberDOFs,
                                            dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs,
                                            dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs + dofs_data[2].NumberDOFs};
    std::vector<unsigned int> offsetStrongs = {0,
                                               dofs_data[0].NumberStrongs,
                                               dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs,
                                               dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs + dofs_data[2].NumberStrongs};
    for (unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }

    PostProcess_Data result;

    result.residual_norm = 0.0;
    if (numberDOFs > 0)
    {
        Gedim::Eigen_Array<> residual;
        residual.SetSize(numberDOFs + 1);
        residual.SumMultiplication(assembler_data.globalMatrixA, assembler_data.solution);
        residual -= assembler_data.rightHandSide;

        result.residual_norm = residual.Norm();
    }

    result.cell0Ds_numeric_pressure.resize(mesh.Cell2DTotalNumber());
    result.cell0Ds_exact_pressure.resize(mesh.Cell2DTotalNumber());

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
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry polygon = {config.GeometricTolerance1D(),
                                                                              config.GeometricTolerance2D(),
                                                                              mesh_geometric_data.Cell2DsVertices.at(c),
                                                                              mesh_geometric_data.Cell2DsCentroids.at(c),
                                                                              mesh_geometric_data.Cell2DsAreas.at(c),
                                                                              mesh_geometric_data.Cell2DsDiameters.at(c),
                                                                              mesh_geometric_data.Cell2DsTriangulations.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeLengths.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeDirections.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeTangents.at(c),
                                                                              mesh_geometric_data.Cell2DsEdgeNormals.at(c)};

        const auto vem_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_local_space(config.VemType());

        const auto local_space = vem_local_space->CreateLocalSpace(velocity_reference_element_data, polygon);

        const auto velocity_basis_functions_derivatives_values =
            vem_local_space->ComputeBasisFunctionsDerivativeValues(local_space, Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);
        const Eigen::MatrixXd pressure_basis_functions_values =
            vem_local_space->ComputePolynomialsValues(local_space).leftCols(local_space.Nkm1);

        const std::vector<size_t> offset_local_dofs = {
            0lu,
            dofs_data[0].CellsGlobalDOFs[2].at(c).size(),
            dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size(),
            dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                dofs_data[2].CellsGlobalDOFs[2].at(c).size()};
        const std::vector<std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData::GlobalCell_DOF>> global_dofs = {
            dofs_data[0].CellsGlobalDOFs[2].at(c),
            dofs_data[1].CellsGlobalDOFs[2].at(c),
            dofs_data[2].CellsGlobalDOFs[2].at(c),
            dofs_data[3].CellsGlobalDOFs[2].at(c)};

        const unsigned int num_velocity_dof = dofs_data[0].CellsGlobalDOFs[2].at(c).size() +
                                              dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                                              dofs_data[2].CellsGlobalDOFs[2].at(c).size();

        Eigen::VectorXd velocity_dofs_values = Eigen::VectorXd::Zero(num_velocity_dof);
        for (unsigned int h = 0; h < numDOFHandler - 1; h++)
        {
            const auto &global_dofs_velocity = dofs_data[h].CellsGlobalDOFs[2].at(c);
            for (unsigned int loc_i = 0; loc_i < global_dofs_velocity.size(); ++loc_i)
            {
                const auto &global_dof_i = global_dofs_velocity.at(loc_i);
                const auto &local_dof_i =
                    dofs_data[h].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    velocity_dofs_values[loc_i + offset_local_dofs[h]] =
                        assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index + offsetStrongs[h]);
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    velocity_dofs_values[loc_i + offset_local_dofs[h]] =
                        assembler_data.solution.GetValue(local_dof_i.Global_Index + offsetDOFs[h]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }

        const auto &pressure_global_dofs = dofs_data[3].CellsGlobalDOFs[2].at(c);
        Eigen::VectorXd pressure_dofs_values = Eigen::VectorXd::Zero(pressure_global_dofs.size());
        for (unsigned int loc_i = 0; loc_i < pressure_global_dofs.size(); ++loc_i)
        {
            const auto &global_dof_i = pressure_global_dofs.at(loc_i);
            const auto &local_dof_i =
                dofs_data[3].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                pressure_dofs_values[loc_i] = assembler_data.solution.GetValue(local_dof_i.Global_Index + offsetDOFs[3]);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }

        const VectorXd numeric_pressure_values = pressure_basis_functions_values * pressure_dofs_values;
        const VectorXd exact_pressure_values = test.exact_pressure(local_space.InternalQuadrature.Points);
        const auto exact_velocity_derivatives_values = test.exact_derivatives_velocity(local_space.InternalQuadrature.Points);

        const Eigen::VectorXd local_error_L2_pressure = (numeric_pressure_values - exact_pressure_values).array().square();
        const Eigen::VectorXd local_norm_L2_pressure = (numeric_pressure_values).array().square();

        result.cell2Ds_error_L2_pressure[c] = local_space.InternalQuadrature.Weights.transpose() * local_error_L2_pressure;
        result.cell2Ds_norm_L2_pressure[c] = local_space.InternalQuadrature.Weights.transpose() * local_norm_L2_pressure;

        const unsigned int numQuadraturePoints = local_space.InternalQuadrature.Points.cols();
        Eigen::VectorXd local_error_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        Eigen::VectorXd local_norm_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        for (unsigned int d1 = 0; d1 < local_space.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < local_space.Dimension; d2++)
            {
                local_error_H1_velocity.array() +=
                    (velocity_basis_functions_derivatives_values[local_space.Dimension * d1 + d2] * velocity_dofs_values -
                     exact_velocity_derivatives_values[3 * d1 + d2])
                        .array()
                        .square();

                local_norm_H1_velocity.array() +=
                    (velocity_basis_functions_derivatives_values[local_space.Dimension * d1 + d2] * velocity_dofs_values)
                        .array()
                        .square();
            }
        }
        result.cell2Ds_error_H1_velocity[c] = local_space.InternalQuadrature.Weights.transpose() * local_error_H1_velocity;
        result.cell2Ds_norm_H1_velocity[c] = local_space.InternalQuadrature.Weights.transpose() * local_norm_H1_velocity;

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
} // namespace Stokes_DF_PCC_2D
} // namespace examples
} // namespace Polydim
