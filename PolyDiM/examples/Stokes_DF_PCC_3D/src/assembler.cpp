#include "assembler.hpp"

#include "Assembler_Utilities.hpp"
#include "EllipticEquation.hpp"
#include "program_configuration.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_3D
{
// ***************************************************************************
void Assembler::ComputeStrongTerm(const Gedim::MeshMatricesDAO &mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                                  const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                                  const std::vector<size_t> &offsetStrongs,
                                  const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
                                  const test::I_Test &test,
                                  Stokes_DF_PCC_3D_Problem_Data &assembler_data) const
{
    //    // Assemble strong boundary condition on Cell0Ds
    //    for (unsigned int h = 0; h < 2; h++)
    //    {
    //        for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
    //        {
    //            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(0).at(p);

    //            if (boundary_info.Type !=
    //            Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
    //                continue;

    //            const auto coordinates = mesh.Cell0DCoordinates(p);

    //            const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker,
    //            coordinates)[h];

    //            const auto local_dofs = dofs_data[h].CellsDOFs.at(0).at(p);

    //            assert(local_dofs.size() == strong_boundary_values.size());

    //            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
    //            {
    //                const auto &local_dof_i = local_dofs.at(loc_i);

    //                switch (local_dof_i.Type)
    //                {
    //                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
    //                    assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index + offsetStrongs[h],
    //                                                              strong_boundary_values[loc_i]);
    //                }
    //                break;
    //                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
    //                    continue;
    //                default:
    //                    throw std::runtime_error("Unknown DOF Type");
    //                }
    //            }
    //        }

    //        // Assemble strong boundary condition on Cell1Ds
    //        const auto &referenceSegmentInternalPoints =
    //        velocity_reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints; const unsigned int
    //        numReferenceSegmentInternalPoints = referenceSegmentInternalPoints.cols();

    //        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    //        {
    //            const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(e);

    //            if (boundary_info.Type !=
    //            Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
    //                continue;

    //            const auto cell1D_origin = mesh.Cell1DOriginCoordinates(e);
    //            const auto cell1D_end = mesh.Cell1DEndCoordinates(e);
    //            const auto cell1D_tangent = cell1D_end - cell1D_origin;

    //            Eigen::MatrixXd coordinates = Eigen::MatrixXd::Zero(3, numReferenceSegmentInternalPoints);
    //            for (unsigned int r = 0; r < numReferenceSegmentInternalPoints; r++)
    //                coordinates.col(r) << cell1D_origin + referenceSegmentInternalPoints(0, r) * cell1D_tangent;

    //            const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker,
    //            coordinates)[h];

    //            const auto local_dofs = dofs_data[h].CellsDOFs.at(1).at(e);

    //            assert(local_dofs.size() == strong_boundary_values.size());

    //            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
    //            {
    //                const auto &local_dof_i = local_dofs.at(loc_i);

    //                switch (local_dof_i.Type)
    //                {
    //                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
    //                    assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index + offsetStrongs[h],
    //                                                              strong_boundary_values[loc_i]);
    //                }
    //                break;
    //                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
    //                    continue;
    //                default:
    //                    throw std::runtime_error("Unknown DOF Type");
    //                }
    //            }
    //        }
    //    }
}
// ***************************************************************************
Assembler::Stokes_DF_PCC_3D_Problem_Data Assembler::Assemble(
    const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
    const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
    const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &velocity_reference_element_data_2D,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data_3D,
    const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data_2D,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Velocity_LocalSpace &vem_velocity_local_space,
    const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Pressure_LocalSpace &vem_pressure_local_space,
    const Polydim::examples::Stokes_DF_PCC_3D::test::I_Test &test) const
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

    Stokes_DF_PCC_3D_Problem_Data result;
    result.globalMatrixA.SetSize(numberDOFs + 1, numberDOFs + 1, Gedim::ISparseArray::SparseArrayTypes::None);
    result.dirichletMatrixA.SetSize(numberDOFs + 1, numberStrongs);
    result.rightHandSide.SetSize(numberDOFs + 1);
    result.solution.SetSize(numberDOFs + 1);
    result.solutionDirichlet.SetSize(numberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        // Assemble equation elements
        for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
        {

            const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(c).size();
            std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> polygonalFaces;
            for (unsigned int f = 0; f < numFaces; f++)
            {
                polygonalFaces.push_back({config.GeometricTolerance1D(),
                                          config.GeometricTolerance2D(),
                                          mesh_geometric_data.Cell3DsFaces2DVertices.at(c)[f],
                                          mesh_geometric_data.Cell3DsFaces2DCentroids.at(c)[f],
                                          mesh_geometric_data.Cell3DsFacesAreas.at(c)[f],
                                          mesh_geometric_data.Cell3DsFacesDiameters.at(c)[f],
                                          mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c)[f],
                                          mesh_geometric_data.Cell3DsFacesEdgeLengths.at(c)[f],
                                          mesh_geometric_data.Cell3DsFacesEdgeDirections.at(c)[f],
                                          mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(c)[f],
                                          mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(c)[f]});
            }

            const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polyhedron_Geometry polyhedron = {
                config.GeometricTolerance1D(),
                config.GeometricTolerance2D(),
                config.GeometricTolerance3D(),
                mesh_geometric_data.Cell3DsVertices.at(c),
                mesh_geometric_data.Cell3DsEdges.at(c),
                mesh_geometric_data.Cell3DsFaces.at(c),
                mesh_geometric_data.Cell3DsCentroids.at(c),
                mesh_geometric_data.Cell3DsVolumes.at(c),
                mesh_geometric_data.Cell3DsDiameters.at(c),
                mesh_geometric_data.Cell3DsTetrahedronPoints.at(c),
                mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c),
                mesh_geometric_data.Cell3DsFacesTranslations.at(c),
                mesh_geometric_data.Cell3DsFacesNormals.at(c),
                mesh_geometric_data.Cell3DsFacesNormalDirections.at(c),
                mesh_geometric_data.Cell3DsFacesNormalGlobalDirection.at(c),
                mesh_geometric_data.Cell3DsFacesTangents.at(c),
                mesh_geometric_data.Cell3DsFacesTangentsGlobalDirection.at(c),
                mesh_geometric_data.Cell3DsEdgeDirections.at(c),
                mesh_geometric_data.Cell3DsEdgeTangents.at(c)};

            const auto pressure_local_space = vem_pressure_local_space.CreateLocalSpace(pressure_reference_element_data, polygon);
            const auto velocity_local_space = vem_velocity_local_space.CreateLocalSpace(velocity_reference_element_data, polygon);

            const auto velocity_basis_functions_values =
                vem_velocity_local_space.ComputeBasisFunctionsValues(velocity_local_space,
                                                                     Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

            const auto velocity_basis_functions_derivatives_values =
                vem_velocity_local_space.ComputeBasisFunctionsDerivativeValues(velocity_local_space,
                                                                               Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);
            const auto velocity_basis_functions_divergence_values =
                vem_velocity_local_space.ComputeBasisFunctionsDivergenceValues(velocity_local_space);

            const MatrixXd pressure_basis_functions_values =
                vem_pressure_local_space.ComputeBasisFunctionsValues(pressure_local_space);

            const auto fluid_viscosity_values = test.fluid_viscosity(velocity_local_space.InternalQuadrature.Points);
            const auto source_term_values = test.source_term(velocity_local_space.InternalQuadrature.Points);

            auto local_A = equation.ComputeCellDiffusionMatrix(fluid_viscosity_values,
                                                               velocity_basis_functions_derivatives_values,
                                                               velocity_local_space.InternalQuadrature.Weights);

            double mu_max = fluid_viscosity_values.cwiseAbs().maxCoeff();
            local_A += mu_max * vem_velocity_local_space.ComputeDofiDofiStabilizationMatrix(
                                    velocity_local_space,
                                    Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

            const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                            velocity_local_space.InternalQuadrature.Weights.asDiagonal() *
                                            velocity_basis_functions_divergence_values;

            const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                                   velocity_basis_functions_values,
                                                                   velocity_local_space.InternalQuadrature.Weights);

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
            elemental_matrix << local_A, local_B.transpose(), local_B,
                MatrixXd::Zero(global_dofs[3].size(), global_dofs[3].size());

            elemental_rhs << local_rhs, VectorXd::Zero(global_dofs[3].size());

            assert(velocity_local_space.NumBasisFunctions ==
                   global_dofs[0].size() + global_dofs[1].size() + global_dofs[2].size());

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
                pressure_basis_functions_values.transpose() * velocity_local_space.InternalQuadrature.Weights;

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
        const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Velocity_LocalSpace &vem_velocity_local_space) const
    {
        Assembler::VEM_Performance_Result result;
        result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

        // Assemble equation elements
        for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
        {
            const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polygon_Geometry polygon = {
                config.GeometricTolerance1D(),
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

            result.Cell2DsPerformance[c].Analysis = performanceAnalysis.Compute(Polydim::VEM::Monomials::VEM_Monomials_3D(),
                                                                                velocity_reference_element_data.Monomials,
                                                                                vem_velocity_local_space,
                                                                                velocity_local_space);

            result.Cell2DsPerformance[c].NumInternalQuadraturePoints = velocity_local_space.InternalQuadrature.Weights.size();
            result.Cell2DsPerformance[c].NumBoundaryQuadraturePoints =
                velocity_local_space.BoundaryQuadrature.Quadrature.Weights.size();
        }

        return result;
    }
    // ***************************************************************************
    Assembler::DiscrepancyErrors_Data Assembler::ComputeDiscrepancyErrors(
        const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &full_dofs_data,
        const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &reduced_dofs_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &full_velocity_reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &full_pressure_reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reduced_velocity_reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reduced_pressure_reference_element_data,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Velocity_LocalSpace &vem_full_velocity_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Pressure_LocalSpace &vem_full_pressure_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Velocity_LocalSpace &vem_reduced_velocity_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Pressure_LocalSpace &vem_reduced_pressure_local_space,
        const Stokes_DF_PCC_3D_Problem_Data &full_assembler_data,
        const Stokes_DF_PCC_3D_Problem_Data &reduced_assembler_data) const
    {

        const unsigned int full_numDOFHandler = full_dofs_data.size();

        unsigned int full_numberDOFs = 0;
        unsigned int full_numberStrongs = 0;
        std::vector<unsigned int> full_offsetDOFs = {0,
                                                     full_dofs_data[0].NumberDOFs,
                                                     full_dofs_data[0].NumberDOFs + full_dofs_data[1].NumberDOFs,
                                                     full_dofs_data[0].NumberDOFs + full_dofs_data[1].NumberDOFs +
                                                         full_dofs_data[2].NumberDOFs};
        std::vector<unsigned int> full_offsetStrongs = {0,
                                                        full_dofs_data[0].NumberStrongs,
                                                        full_dofs_data[0].NumberStrongs + full_dofs_data[1].NumberStrongs,
                                                        full_dofs_data[0].NumberStrongs + full_dofs_data[1].NumberStrongs +
                                                            full_dofs_data[2].NumberStrongs};
        for (unsigned int i = 0; i < full_numDOFHandler; i++)
        {
            full_numberDOFs += full_dofs_data[i].NumberDOFs;
            full_numberStrongs += full_dofs_data[i].NumberStrongs;
        }

        const unsigned int reduced_numDOFHandler = reduced_dofs_data.size();

        unsigned int reduced_numberDOFs = 0;
        unsigned int reduced_numberStrongs = 0;
        std::vector<unsigned int> reduced_offsetDOFs = {0,
                                                        reduced_dofs_data[0].NumberDOFs,
                                                        reduced_dofs_data[0].NumberDOFs + reduced_dofs_data[1].NumberDOFs,
                                                        reduced_dofs_data[0].NumberDOFs + reduced_dofs_data[1].NumberDOFs +
                                                            reduced_dofs_data[2].NumberDOFs};
        std::vector<unsigned int> reduced_offsetStrongs = {
            0,
            reduced_dofs_data[0].NumberStrongs,
            reduced_dofs_data[0].NumberStrongs + reduced_dofs_data[1].NumberStrongs,
            reduced_dofs_data[0].NumberStrongs + reduced_dofs_data[1].NumberStrongs + reduced_dofs_data[2].NumberStrongs};
        for (unsigned int i = 0; i < reduced_numDOFHandler; i++)
        {
            reduced_numberDOFs += reduced_dofs_data[i].NumberDOFs;
            reduced_numberStrongs += reduced_dofs_data[i].NumberStrongs;
        }

        DiscrepancyErrors_Data result;

        result.residual_norm = 0.0;
        if (full_numberDOFs > 0)
        {
            Gedim::Eigen_Array<> residual;
            residual.SetSize(full_numberDOFs + 1);
            residual.SumMultiplication(full_assembler_data.globalMatrixA, full_assembler_data.solution);
            residual -= full_assembler_data.rightHandSide;

            result.residual_norm = residual.Norm();
        }

        result.reduced_residual_norm = 0.0;
        if (reduced_numberDOFs > 0)
        {
            Gedim::Eigen_Array<> residual;
            residual.SetSize(reduced_numberDOFs + 1);
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
            const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polygon_Geometry polygon = {
                config.GeometricTolerance1D(),
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

            const std::vector<size_t> full_offset_local_dofs = {
                0lu,
                full_dofs_data[0].CellsGlobalDOFs[2].at(c).size(),
                full_dofs_data[0].CellsGlobalDOFs[2].at(c).size() + full_dofs_data[1].CellsGlobalDOFs[2].at(c).size(),
                full_dofs_data[0].CellsGlobalDOFs[2].at(c).size() + full_dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                    full_dofs_data[2].CellsGlobalDOFs[2].at(c).size()};

            const unsigned int num_full_velocity_dof = full_dofs_data[0].CellsGlobalDOFs[2].at(c).size() +
                                                       full_dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                                                       full_dofs_data[2].CellsGlobalDOFs[2].at(c).size();

            Eigen::VectorXd full_velocity_dofs_values = Eigen::VectorXd::Zero(num_full_velocity_dof);
            for (unsigned int h = 0; h < full_numDOFHandler - 1; h++)
            {
                const auto &global_dofs_velocity = full_dofs_data[h].CellsGlobalDOFs[2].at(c);
                for (unsigned int loc_i = 0; loc_i < global_dofs_velocity.size(); ++loc_i)
                {
                    const auto &global_dof_i = global_dofs_velocity.at(loc_i);
                    const auto &local_dof_i = full_dofs_data[h]
                                                  .CellsDOFs.at(global_dof_i.Dimension)
                                                  .at(global_dof_i.CellIndex)
                                                  .at(global_dof_i.DOFIndex);

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                        full_velocity_dofs_values[loc_i + full_offset_local_dofs[h]] =
                            full_assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index + full_offsetStrongs[h]);
                        break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                        full_velocity_dofs_values[loc_i + full_offset_local_dofs[h]] =
                            full_assembler_data.solution.GetValue(local_dof_i.Global_Index + full_offsetDOFs[h]);
                        break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }

            const auto &full_pressure_global_dofs = full_dofs_data[3].CellsGlobalDOFs[2].at(c);
            Eigen::VectorXd full_pressure_dofs_values = Eigen::VectorXd::Zero(full_pressure_global_dofs.size());
            for (unsigned int loc_i = 0; loc_i < full_pressure_global_dofs.size(); ++loc_i)
            {
                const auto &global_dof_i = full_pressure_global_dofs.at(loc_i);
                const auto &local_dof_i =
                    full_dofs_data[3].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    full_pressure_dofs_values[loc_i] =
                        full_assembler_data.solution.GetValue(local_dof_i.Global_Index + full_offsetDOFs[3]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            const std::vector<size_t> reduced_offset_local_dofs = {
                0lu,
                reduced_dofs_data[0].CellsGlobalDOFs[2].at(c).size(),
                reduced_dofs_data[0].CellsGlobalDOFs[2].at(c).size() + reduced_dofs_data[1].CellsGlobalDOFs[2].at(c).size(),
                reduced_dofs_data[0].CellsGlobalDOFs[2].at(c).size() + reduced_dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                    reduced_dofs_data[2].CellsGlobalDOFs[2].at(c).size()};

            const unsigned int num_reduced_velocity_dof = reduced_dofs_data[0].CellsGlobalDOFs[2].at(c).size() +
                                                          reduced_dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                                                          reduced_dofs_data[2].CellsGlobalDOFs[2].at(c).size();

            Eigen::VectorXd reduced_velocity_dofs_values = Eigen::VectorXd::Zero(num_reduced_velocity_dof);
            for (unsigned int h = 0; h < reduced_numDOFHandler - 1; h++)
            {
                const auto &global_dofs_velocity = reduced_dofs_data[h].CellsGlobalDOFs[2].at(c);
                for (unsigned int loc_i = 0; loc_i < global_dofs_velocity.size(); ++loc_i)
                {
                    const auto &global_dof_i = global_dofs_velocity.at(loc_i);
                    const auto &local_dof_i = reduced_dofs_data[h]
                                                  .CellsDOFs.at(global_dof_i.Dimension)
                                                  .at(global_dof_i.CellIndex)
                                                  .at(global_dof_i.DOFIndex);

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                        reduced_velocity_dofs_values[loc_i + reduced_offset_local_dofs[h]] =
                            reduced_assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index + reduced_offsetStrongs[h]);
                        break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                        reduced_velocity_dofs_values[loc_i + reduced_offset_local_dofs[h]] =
                            reduced_assembler_data.solution.GetValue(local_dof_i.Global_Index + reduced_offsetDOFs[h]);
                        break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }

            const auto &reduced_pressure_global_dofs = reduced_dofs_data[3].CellsGlobalDOFs[2].at(c);
            Eigen::VectorXd reduced_pressure_dofs_values = Eigen::VectorXd::Zero(reduced_pressure_global_dofs.size());
            for (unsigned int loc_i = 0; loc_i < reduced_pressure_global_dofs.size(); ++loc_i)
            {
                const auto &global_dof_i = reduced_pressure_global_dofs.at(loc_i);
                const auto &local_dof_i =
                    reduced_dofs_data[3].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    reduced_pressure_dofs_values[loc_i] =
                        reduced_assembler_data.solution.GetValue(local_dof_i.Global_Index + reduced_offsetDOFs[3]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

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
            result.cell2Ds_full_norm_L2_pressure[c] =
                full_velocity_local_space.InternalQuadrature.Weights.transpose() * local_norm_L2_pressure;

            const unsigned int numQuadraturePoints = full_velocity_local_space.InternalQuadrature.Points.cols();
            Eigen::VectorXd local_error_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
            Eigen::VectorXd local_norm_H1_velocity = Eigen::VectorXd::Zero(numQuadraturePoints);
            for (unsigned int d1 = 0; d1 < full_velocity_local_space.Dimension; d1++)
            {
                for (unsigned int d2 = 0; d2 < full_velocity_local_space.Dimension; d2++)
                {
                    local_error_H1_velocity.array() +=
                        (full_velocity_basis_functions_derivatives_values[full_velocity_local_space.Dimension * d1 + d2] * full_velocity_dofs_values -
                         reduced_velocity_basis_functions_derivatives_values[reduced_velocity_local_space.Dimension * d1 + d2] *
                             reduced_velocity_dofs_values)
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
            result.cell2Ds_full_norm_H1_velocity[c] =
                full_velocity_local_space.InternalQuadrature.Weights.transpose() * local_norm_H1_velocity;
        }

        result.discrepancy_error_L2_pressure = std::sqrt(result.cell2Ds_discrepancy_error_L2_pressure.sum());
        result.discrepancy_error_H1_velocity = std::sqrt(result.cell2Ds_discrepancy_error_H1_velocity.sum());
        result.full_norm_H1_velocity = std::sqrt(result.cell2Ds_full_norm_H1_velocity.sum());
        result.full_norm_L2_pressure = std::sqrt(result.cell2Ds_full_norm_L2_pressure.sum());

        result.pressure_dofs_ratio = ((double)reduced_dofs_data[3].NumberDOFs) / full_dofs_data[3].NumberDOFs;
        result.velocity_dofs_ratio = ((double)reduced_offsetDOFs[3]) / full_offsetDOFs[3];

        return result;
    }
    // ***************************************************************************
    Assembler::PostProcess_Data Assembler::PostProcessSolution(
        const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
        const Gedim::MeshMatricesDAO &mesh,
        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
        const vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &velocity_reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &pressure_reference_element_data,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Velocity_LocalSpace &vem_velocity_local_space,
        const Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Pressure_LocalSpace &vem_pressure_local_space,
        const Stokes_DF_PCC_3D_Problem_Data &assembler_data,
        const Polydim::examples::Stokes_DF_PCC_3D::test::I_Test &test) const
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
                                                   dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs +
                                                       dofs_data[2].NumberStrongs};
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
                        result.cell0Ds_exact_velocity[d][p] =
                            assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index + offsetStrongs[d]);
                        break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                        result.cell0Ds_exact_velocity[d][p] =
                            assembler_data.solution.GetValue(local_dof_i.Global_Index + offsetDOFs[d]);
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
            const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polygon_Geometry polygon = {
                config.GeometricTolerance1D(),
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

            const std::vector<size_t> offset_local_dofs = {
                0lu,
                dofs_data[0].CellsGlobalDOFs[2].at(c).size(),
                dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size(),
                dofs_data[0].CellsGlobalDOFs[2].at(c).size() + dofs_data[1].CellsGlobalDOFs[2].at(c).size() +
                    dofs_data[2].CellsGlobalDOFs[2].at(c).size()};

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
} // namespace Stokes_DF_PCC_3D
} // namespace examples
} // namespace Polydim
