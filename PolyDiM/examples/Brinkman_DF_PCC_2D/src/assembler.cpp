#include "assembler.hpp"

#include "EllipticEquation.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "VEM_DF_PCC_PerformanceAnalysis.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{
// ***************************************************************************
void Assembler::ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                  const std::vector<size_t> &offsetStrongs,
                                  const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
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
        const auto &referenceSegmentInternalPoints = velocity_reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints;
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
void Assembler::ComputeWeakTerm(const unsigned int cell2DIndex,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_local_space,
                                const Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test &test,
                                Stokes_DF_PCC_2D_Problem_Data &assembler_data) const
{
    const unsigned numVertices = polygon.Vertices.cols();
    for (unsigned int h = 0; h < reference_element_data.Dimension; h++)
    {
        for (unsigned int ed = 0; ed < numVertices; ed++)
        {
            const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex, ed);

            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(cell1D_index);

            if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
                continue;

            // compute vem values
            const auto weakReferenceSegment =
                Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

            const Eigen::VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0);

            const auto weak_basis_function_values =
                vem_local_space.ComputeValuesOnEdge(reference_element_data, pointsCurvilinearCoordinates);

            // map edge internal quadrature points
            const Eigen::Vector3d &edgeStart = polygon.EdgesDirection[ed] ? polygon.Vertices.col(ed)
                                                                          : polygon.Vertices.col((ed + 1) % numVertices);

            const Eigen::Vector3d &edgeTangent = polygon.EdgesTangent.col(ed);
            const double direction = polygon.EdgesDirection[ed] ? 1.0 : -1.0;

            const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegment.Points.cols();
            Eigen::MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
            for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
            {
                weakQuadraturePoints.col(q) = edgeStart + direction * weakReferenceSegment.Points(0, q) * edgeTangent;
            }
            const double absMapDeterminant = std::abs(polygon.EdgesLength[ed]);
            const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights * absMapDeterminant;

            const auto neumannValues = test.weak_boundary_condition(boundary_info.Marker, weakQuadraturePoints);

            // compute values of Neumann condition
            const Eigen::VectorXd neumannContributions =
                -weak_basis_function_values.transpose() * weakQuadratureWeights.asDiagonal() * neumannValues[h];

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
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_pressure_local_space,
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

        const auto pressure_local_space = vem_pressure_local_space.CreateLocalSpace(pressure_reference_element_data, polygon);
        const auto velocity_local_space = vem_velocity_local_space.CreateLocalSpace(velocity_reference_element_data, polygon);

        const auto velocity_basis_functions_values =
            vem_velocity_local_space.ComputeBasisFunctionsValues(velocity_local_space, Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

        const auto velocity_basis_functions_derivatives_values =
            vem_velocity_local_space.ComputeBasisFunctionsDerivativeValues(velocity_local_space,
                                                                           Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);
        const auto velocity_basis_functions_divergence_values =
            vem_velocity_local_space.ComputeBasisFunctionsDivergenceValues(velocity_local_space);

        const MatrixXd pressure_basis_functions_values = vem_pressure_local_space.ComputeBasisFunctionsValues(pressure_local_space);

        const auto fluid_viscosity_values = test.fluid_viscosity(velocity_local_space.InternalQuadrature.Points);
        const auto inverse_diffusion_term_values = test.inverse_diffusion_term(velocity_local_space.InternalQuadrature.Points);
        const auto source_term_values = test.source_term(velocity_local_space.InternalQuadrature.Points);
        const auto divergence_term_values = test.divergence_term(velocity_local_space.InternalQuadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(fluid_viscosity_values,
                                                           velocity_basis_functions_derivatives_values,
                                                           velocity_local_space.InternalQuadrature.Weights);

        const double mu_max = fluid_viscosity_values.cwiseAbs().maxCoeff();
        local_A += mu_max * vem_velocity_local_space.ComputeDofiDofiStabilizationMatrix(velocity_local_space,
                                                                                        Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

        local_A += equation.ComputeCellDiffusionMatrix(inverse_diffusion_term_values,
                                                       velocity_basis_functions_values,
                                                       velocity_local_space.InternalQuadrature.Weights);

        double k_max = 0.0;
        for (const auto &inv_diffusion_term : inverse_diffusion_term_values)
        {
            const double max_k = inv_diffusion_term.cwiseAbs().maxCoeff();
            k_max = k_max < max_k ? max_k : k_max;
        }
        local_A += k_max * vem_velocity_local_space.ComputeDofiDofiStabilizationMatrix(velocity_local_space,
                                                                                       Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        velocity_local_space.InternalQuadrature.Weights.asDiagonal() *
                                        velocity_basis_functions_divergence_values;

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               velocity_basis_functions_values,
                                                               velocity_local_space.InternalQuadrature.Weights);

        const auto local_div = equation.ComputeCellForcingTerm(divergence_term_values,
                                                               pressure_basis_functions_values,
                                                               velocity_local_space.InternalQuadrature.Weights);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[3].CellsGlobalDOFs[2].at(c).size();

        Eigen::MatrixXd elemental_matrix = MatrixXd::Zero(local_count_dofs.num_total_dofs, local_count_dofs.num_total_dofs);
        Eigen::VectorXd elemental_rhs = VectorXd::Zero(local_count_dofs.num_total_dofs);

        elemental_matrix << local_A, -local_B.transpose(), local_B, MatrixXd::Zero(num_local_dofs_pressure, num_local_dofs_pressure);

        elemental_rhs << local_rhs, local_div;

        assert(velocity_local_space.NumBasisFunctions == local_count_dofs.num_total_dofs - num_local_dofs_pressure);

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
            const VectorXd mean_value_pressure =
                pressure_basis_functions_values.transpose() * velocity_local_space.InternalQuadrature.Weights;

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

        ComputeWeakTerm(c, mesh, polygon, mesh_dofs_info, dofs_data, count_dofs, velocity_reference_element_data, vem_velocity_local_space, test, result);
    }

    ComputeStrongTerm(mesh, mesh_geometric_data, mesh_dofs_info, dofs_data, count_dofs.offsets_Strongs, velocity_reference_element_data, test, result);

    result.rightHandSide.Create();
    result.solutionDirichlet.Create();
    result.globalMatrixA.Create();
    result.dirichletMatrixA.Create();

    if (count_dofs.num_total_strong > 0)
        result.rightHandSide.SubtractionMultiplication(result.dirichletMatrixA, result.solutionDirichlet);

    return result;
}
// ***************************************************************************
Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(
    const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space) const
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

        const auto velocity_local_space = vem_velocity_local_space.CreateLocalSpace(velocity_reference_element_data, polygon);

        Polydim::VEM::DF_PCC::VEM_DF_PCC_PerformanceAnalysis performanceAnalysis;

        const auto Analysis = performanceAnalysis.Compute(Polydim::VEM::Utilities::VEM_Monomials_2D(),
                                                          velocity_reference_element_data.Monomials,
                                                          vem_velocity_local_space,
                                                          velocity_local_space);

        result.Cell2DsPerformance[c].maxErrorGBD = *std::max_element(Analysis.ErrorGBD.begin(), Analysis.ErrorGBD.end());
        result.Cell2DsPerformance[c].maxErrorHCD = *std::max_element(Analysis.ErrorHCD.begin(), Analysis.ErrorHCD.end());
        result.Cell2DsPerformance[c].maxErrorPiNabla =
            *std::max_element(Analysis.ErrorPiNabla.begin(), Analysis.ErrorPiNabla.end());
        result.Cell2DsPerformance[c].maxErrorPi0k = *std::max_element(Analysis.ErrorPi0k.begin(), Analysis.ErrorPi0k.end());
        result.Cell2DsPerformance[c].ErrorStabilization = Analysis.ErrorStabilization;
        result.Cell2DsPerformance[c].maxPiNablaConditioning =
            *std::max_element(Analysis.PiNablaConditioning.begin(), Analysis.PiNablaConditioning.end());
        result.Cell2DsPerformance[c].maxPi0kConditioning =
            *std::max_element(Analysis.Pi0kConditioning.begin(), Analysis.Pi0kConditioning.end());
        result.Cell2DsPerformance[c].NumInternalQuadraturePoints = velocity_local_space.InternalQuadrature.Weights.size();
        result.Cell2DsPerformance[c].NumBoundaryQuadraturePoints =
            velocity_local_space.BoundaryQuadrature.Quadrature.Weights.size();
    }

    return result;
}
// ***************************************************************************
Assembler::DiscrepancyErrors_Data Assembler::ComputeDiscrepancyErrors(
    const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const std::vector<PDETools::DOFs::DOFsManager::DOFsData> &full_dofs_data,
    const Polydim::PDETools::Assembler_Utilities::count_dofs_data &full_count_dofs,
    const vector<PDETools::DOFs::DOFsManager::DOFsData> &reduced_dofs_data,
    const Polydim::PDETools::Assembler_Utilities::count_dofs_data &reduced_count_dofs,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &full_velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &full_pressure_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reduced_velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reduced_pressure_reference_element_data,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_full_velocity_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_full_pressure_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_reduced_velocity_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_reduced_pressure_local_space,
    const Stokes_DF_PCC_2D_Problem_Data &full_assembler_data,
    const Stokes_DF_PCC_2D_Problem_Data &reduced_assembler_data) const
{

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

        const auto full_pressure_local_space =
            vem_full_pressure_local_space.CreateLocalSpace(full_pressure_reference_element_data, polygon);

        const auto full_velocity_local_space =
            vem_full_velocity_local_space.CreateLocalSpace(full_velocity_reference_element_data, polygon);

        const auto full_velocity_basis_functions_derivatives_values =
            vem_full_velocity_local_space.ComputeBasisFunctionsDerivativeValues(full_velocity_local_space,
                                                                                Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);
        const Eigen::MatrixXd full_pressure_basis_functions_values =
            vem_full_pressure_local_space.ComputeBasisFunctionsValues(full_pressure_local_space);

        const auto reduced_pressure_local_space =
            vem_reduced_pressure_local_space.CreateLocalSpace(reduced_pressure_reference_element_data, polygon);

        const auto reduced_velocity_local_space =
            vem_reduced_velocity_local_space.CreateLocalSpace(reduced_velocity_reference_element_data, polygon);

        const auto reduced_velocity_basis_functions_derivatives_values =
            vem_reduced_velocity_local_space.ComputeBasisFunctionsDerivativeValues(reduced_velocity_local_space,
                                                                                   Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);
        const Eigen::MatrixXd reduced_pressure_basis_functions_values =
            vem_reduced_pressure_local_space.ComputeBasisFunctionsValues(reduced_pressure_local_space);

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
            ((1.0 / polygon.Measure) * full_numeric_pressure_values.transpose() *
             full_pressure_local_space.InternalQuadrature.Weights) *
            reduced_pressure_basis_functions_values;

        const Eigen::VectorXd local_error_L2_pressure =
            (full_numeric_projected_pressure_values - reduced_numeric_pressure_values).array().square();
        const Eigen::VectorXd local_norm_L2_pressure = (full_numeric_projected_pressure_values).array().square();

        result.cell2Ds_discrepancy_error_L2_pressure[c] =
            full_velocity_local_space.InternalQuadrature.Weights.transpose() * local_error_L2_pressure;
        result.cell2Ds_full_norm_L2_pressure[c] = full_velocity_local_space.InternalQuadrature.Weights.transpose() * local_norm_L2_pressure;

        const unsigned int numQuadraturePoints = full_velocity_local_space.InternalQuadrature.Points.cols();
        Eigen::VectorXd local_error_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        Eigen::VectorXd local_norm_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        for (unsigned int d1 = 0; d1 < full_velocity_local_space.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < full_velocity_local_space.Dimension; d2++)
            {
                local_error_H1_velocity.array() +=
                    (full_velocity_basis_functions_derivatives_values[full_velocity_local_space.Dimension * d1 + d2] * full_velocity_dofs_values -
                     reduced_velocity_basis_functions_derivatives_values[reduced_velocity_local_space.Dimension * d1 + d2] * reduced_velocity_dofs_values)
                        .array()
                        .square();

                local_norm_H1_velocity.array() +=
                    (full_velocity_basis_functions_derivatives_values[full_velocity_local_space.Dimension * d1 + d2] * full_velocity_dofs_values)
                        .array()
                        .square();
            }
        }
        result.cell2Ds_discrepancy_error_H1_velocity[c] =
            full_velocity_local_space.InternalQuadrature.Weights.transpose() * local_error_H1_velocity;
        result.cell2Ds_full_norm_H1_velocity[c] = full_velocity_local_space.InternalQuadrature.Weights.transpose() * local_norm_H1_velocity;
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
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_pressure_local_space,
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

    for (unsigned int d = 0; d < velocity_reference_element_data.Dimension; d++)
    {
        result.cell0Ds_numeric_velocity[d].resize(mesh.Cell0DTotalNumber());
        result.cell0Ds_exact_velocity[d].resize(mesh.Cell0DTotalNumber());
    }

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        const auto exact_velocity = test.exact_velocity(mesh.Cell0DCoordinates(p));

        for (unsigned int d = 0; d < velocity_reference_element_data.Dimension; d++)
            result.cell0Ds_exact_velocity[d](p) = exact_velocity[d](0);

        for (unsigned int d = 0; d < velocity_reference_element_data.Dimension; d++)
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
    for(unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        for(unsigned int v = 0; v < mesh.Cell2DNumberVertices(c); v++)
        {
            result.repeated_connectivity[c][v] = num_repetead_points;
            num_repetead_points++;
        }
    }
    result.repeated_vertices_coordinates = Eigen::MatrixXd::Zero(3, num_repetead_points);
    result.cell0Ds_exact_pressure = Eigen::VectorXd::Zero(num_repetead_points);
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


    num_repetead_points = 0;
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

        const auto pressure_local_space = vem_pressure_local_space.CreateLocalSpace(pressure_reference_element_data, polygon);
        const auto velocity_local_space = vem_velocity_local_space.CreateLocalSpace(velocity_reference_element_data, polygon);

        const auto velocity_basis_functions_derivatives_values =
            vem_velocity_local_space.ComputeBasisFunctionsDerivativeValues(velocity_local_space,
                                                                           Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);
        const Eigen::MatrixXd pressure_basis_functions_values =
            vem_pressure_local_space.ComputeBasisFunctionsValues(pressure_local_space);

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

        const unsigned int num_cell_vertices =  mesh.Cell2DNumberVertices(c);
        const Eigen::MatrixXd pressure_basis_functions_vertices_values =
            vem_pressure_local_space.ComputeBasisFunctionsValues(pressure_reference_element_data, pressure_local_space, mesh_geometric_data.Cell2DsVertices.at(c));
        result.repeated_vertices_coordinates.middleCols(num_repetead_points, num_cell_vertices) = mesh_geometric_data.Cell2DsVertices.at(c);
        result.cell0Ds_numeric_pressure.segment(num_repetead_points, num_cell_vertices) = pressure_basis_functions_vertices_values * pressure_dofs_values;
        result.cell0Ds_exact_pressure.segment(num_repetead_points, num_cell_vertices) = test.exact_pressure(mesh_geometric_data.Cell2DsVertices.at(c));
        num_repetead_points += num_cell_vertices;

        const VectorXd numeric_pressure_values = pressure_basis_functions_values * pressure_dofs_values;
        const VectorXd exact_pressure_values = test.exact_pressure(velocity_local_space.InternalQuadrature.Points);
        const auto exact_velocity_derivatives_values =
            test.exact_derivatives_velocity(velocity_local_space.InternalQuadrature.Points);

        const Eigen::VectorXd local_error_L2_pressure = (numeric_pressure_values - exact_pressure_values).array().square();
        const Eigen::VectorXd local_norm_L2_pressure = (numeric_pressure_values).array().square();

        result.cell2Ds_error_L2_pressure[c] = velocity_local_space.InternalQuadrature.Weights.transpose() * local_error_L2_pressure;
        result.cell2Ds_norm_L2_pressure[c] = velocity_local_space.InternalQuadrature.Weights.transpose() * local_norm_L2_pressure;

        const unsigned int numQuadraturePoints = velocity_local_space.InternalQuadrature.Points.cols();
        Eigen::VectorXd local_error_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        Eigen::VectorXd local_norm_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
        for (unsigned int d1 = 0; d1 < velocity_local_space.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < velocity_local_space.Dimension; d2++)
            {
                local_error_H1_velocity.array() +=
                    (velocity_basis_functions_derivatives_values[velocity_local_space.Dimension * d1 + d2] * velocity_dofs_values -
                     exact_velocity_derivatives_values[3 * d1 + d2])
                        .array()
                        .square();

                local_norm_H1_velocity.array() +=
                    (velocity_basis_functions_derivatives_values[velocity_local_space.Dimension * d1 + d2] * velocity_dofs_values)
                        .array()
                        .square();
            }
        }
        result.cell2Ds_error_H1_velocity[c] = velocity_local_space.InternalQuadrature.Weights.transpose() * local_error_H1_velocity;
        result.cell2Ds_norm_H1_velocity[c] = velocity_local_space.InternalQuadrature.Weights.transpose() * local_norm_H1_velocity;

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
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim
