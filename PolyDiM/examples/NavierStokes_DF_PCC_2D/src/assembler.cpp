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
namespace NavierStokes_DF_PCC_2D
{
// ***************************************************************************
void Assembler::ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                  const std::vector<size_t> &offsetStrongs,
                                  const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                  const test::I_Test &test,
                                  NavierStokes_DF_PCC_2D_Problem_Data &assembler_data) const
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
                                const Polydim::examples::NavierStokes_DF_PCC_2D::test::I_Test &test,
                                NavierStokes_DF_PCC_2D_Problem_Data &assembler_data) const
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
                weak_basis_function_values.transpose() * weakQuadratureWeights.asDiagonal() * neumannValues[h];

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
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_pressure_local_space,
    const Polydim::examples::NavierStokes_DF_PCC_2D::test::I_Test &test) const
{
    NavierStokes_DF_PCC_2D_Problem_Data result;
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
        const auto source_term_values = test.source_term(velocity_local_space.InternalQuadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(fluid_viscosity_values,
                                                           velocity_basis_functions_derivatives_values,
                                                           velocity_local_space.InternalQuadrature.Weights);

        double mu_max = fluid_viscosity_values.cwiseAbs().maxCoeff();
        local_A += mu_max * vem_velocity_local_space.ComputeDofiDofiStabilizationMatrix(velocity_local_space,
                                                                                        Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        velocity_local_space.InternalQuadrature.Weights.asDiagonal() *
                                        velocity_basis_functions_divergence_values;

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               velocity_basis_functions_values,
                                                               velocity_local_space.InternalQuadrature.Weights);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[3].CellsGlobalDOFs[2].at(c).size();

        Eigen::MatrixXd elemental_matrix = MatrixXd::Zero(local_count_dofs.num_total_dofs, local_count_dofs.num_total_dofs);
        Eigen::VectorXd elemental_rhs = VectorXd::Zero(local_count_dofs.num_total_dofs);

        elemental_matrix << local_A, local_B.transpose(), local_B, MatrixXd::Zero(num_local_dofs_pressure, num_local_dofs_pressure);
        elemental_rhs << local_rhs, VectorXd::Zero(num_local_dofs_pressure);

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
void Assembler::AssembleNavierStokes(const Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration &config,
                                     const Gedim::MeshMatricesDAO &mesh,
                                     const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                     const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                                     const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                     const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
                                     const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space,
                                     const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_pressure_local_space,
                                     NavierStokes_DF_PCC_2D_Problem_Data &result)
{

    result.globalMatrixC.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    result.dirichletMatrixC.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    result.rightHandSideC.SetSize(count_dofs.num_total_dofs);

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

        const auto pressure_local_space = vem_pressure_local_space.CreateLocalSpace(pressure_reference_element_data, polygon);
        const auto velocity_local_space = vem_velocity_local_space.CreateLocalSpace(velocity_reference_element_data, polygon);

        const auto velocity_basis_functions_values =
            vem_velocity_local_space.ComputeBasisFunctionsValues(velocity_local_space, Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

        const auto velocity_basis_functions_derivatives_values =
            vem_velocity_local_space.ComputeBasisFunctionsDerivativeValues(velocity_local_space,
                                                                           Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[3].CellsGlobalDOFs[2].at(c).size();
        const unsigned int num_local_dofs_velocity = local_count_dofs.num_total_dofs - num_local_dofs_pressure;
        MatrixXd cellMatrixC = MatrixXd::Zero(num_local_dofs_velocity, num_local_dofs_velocity);
        VectorXd cellRightHandSideC = VectorXd::Zero(num_local_dofs_velocity);

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
        switch (config.ConvectiveForm())
        {
        case Program_configuration::ConvectiveFormType::None:
            break;
        case Program_configuration::ConvectiveFormType::Conv: {
            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellMatrixC +=
                        velocity_basis_functions_values[d1].transpose() *
                        velocity_local_space.InternalQuadrature.Weights.cwiseProduct(previous_iteration_values[d2]).asDiagonal() *
                        velocity_basis_functions_derivatives_values[2 * d1 + d2];
            }

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellMatrixC += velocity_basis_functions_values[d1].transpose() *
                                   velocity_local_space.InternalQuadrature.Weights
                                       .cwiseProduct(previous_iteration_derivatives_values[2 * d1 + d2])
                                       .asDiagonal() *
                                   velocity_basis_functions_values[d2];
            }

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellRightHandSideC +=
                        velocity_basis_functions_values[d1].transpose() *
                        velocity_local_space.InternalQuadrature.Weights.asDiagonal() *
                        previous_iteration_values[d2].cwiseProduct(previous_iteration_derivatives_values[2 * d1 + d2]);
            }
        }
        break;
        case Program_configuration::ConvectiveFormType::Skew: {
            MatrixXd cellCMatrix1 = MatrixXd::Zero(num_local_dofs_velocity, num_local_dofs_velocity);
            MatrixXd cellCMatrix2 = MatrixXd::Zero(num_local_dofs_velocity, num_local_dofs_velocity);

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix1 +=
                        velocity_basis_functions_values[d1].transpose() *
                        velocity_local_space.InternalQuadrature.Weights.cwiseProduct(previous_iteration_values[d2]).asDiagonal() *
                        velocity_basis_functions_derivatives_values[2 * d1 + d2];
            }

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix1 += velocity_basis_functions_values[d1].transpose() *
                                    velocity_local_space.InternalQuadrature.Weights
                                        .cwiseProduct(previous_iteration_derivatives_values[2 * d1 + d2])
                                        .asDiagonal() *
                                    velocity_basis_functions_values[d2];
            }

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix2 +=
                        velocity_basis_functions_derivatives_values[2 * d1 + d2].transpose() *
                        velocity_local_space.InternalQuadrature.Weights.cwiseProduct(previous_iteration_values[d2]).asDiagonal() *
                        velocity_basis_functions_values[d1];
            }

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellCMatrix2 +=
                        velocity_basis_functions_derivatives_values[2 * d1 + d2].transpose() *
                        velocity_local_space.InternalQuadrature.Weights.cwiseProduct(previous_iteration_values[d1]).asDiagonal() *
                        velocity_basis_functions_values[d2];
            }

            cellMatrixC = 0.5 * (cellCMatrix1 - cellCMatrix2);

            VectorXd cellRightHandSide1 = VectorXd::Zero(num_local_dofs_velocity);
            VectorXd cellRightHandSide2 = VectorXd::Zero(num_local_dofs_velocity);

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellRightHandSide1 +=
                        velocity_basis_functions_values[d1].transpose() *
                        velocity_local_space.InternalQuadrature.Weights.asDiagonal() *
                        previous_iteration_values[d2].cwiseProduct(previous_iteration_derivatives_values[2 * d1 + d2]);
            }

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                    cellRightHandSide2 += velocity_basis_functions_derivatives_values[2 * d1 + d2].transpose() *
                                          velocity_local_space.InternalQuadrature.Weights.asDiagonal() *
                                          previous_iteration_values[d2].cwiseProduct(previous_iteration_values[d1]);
            }

            cellRightHandSideC = 0.5 * (cellRightHandSide1 - cellRightHandSide2);
        }
        break;
        default:
            throw runtime_error("Not valid convective form type");
        }

        assert(velocity_local_space.NumBasisFunctions == num_local_dofs_velocity);

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
Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(
    const Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration &config,
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

        const auto Analysis = performanceAnalysis.Compute(Polydim::VEM::Monomials::VEM_Monomials_2D(),
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

        result.Cell2DsPerformance[c].NumInternalQuadraturePoints = velocity_local_space.InternalQuadrature.Weights.size();
        result.Cell2DsPerformance[c].NumBoundaryQuadraturePoints =
            velocity_local_space.BoundaryQuadrature.Quadrature.Weights.size();
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
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace &vem_velocity_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace &vem_pressure_local_space,
    const NavierStokes_DF_PCC_2D_Problem_Data &assembler_data,
    const Polydim::examples::NavierStokes_DF_PCC_2D::test::I_Test &test) const
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

        const Eigen::VectorXd velocity_dofs_values =
            dofs_values.segment(0, local_count_dofs.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd pressure_dofs_values =
            dofs_values.segment(local_count_dofs.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

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
} // namespace NavierStokes_DF_PCC_2D
} // namespace examples
} // namespace Polydim
