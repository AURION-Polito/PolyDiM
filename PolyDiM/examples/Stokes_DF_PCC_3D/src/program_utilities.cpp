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

#include "program_utilities.hpp"

#include "Eigen_BiCGSTABSolver.hpp"
#include "Eigen_LUSolver.hpp"
#include "VTKUtilities.hpp"
#include <numbers>

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_3D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Stokes_DF_PCC_3D::test::I_Test> create_test(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Stokes_DF_PCC_3D::test::Test_Types::Patch_Test:
        return std::make_unique<Polydim::examples::Stokes_DF_PCC_3D::test::Patch_Test>();
    case Polydim::examples::Stokes_DF_PCC_3D::test::Test_Types::Stokes:
        return std::make_unique<Polydim::examples::Stokes_DF_PCC_3D::test::Stokes>();
    case Polydim::examples::Stokes_DF_PCC_3D::test::Test_Types::Stokes_Benchmark_1:
        return std::make_unique<Polydim::examples::Stokes_DF_PCC_3D::test::Stokes_Benchmark_1>();
    case Polydim::examples::Stokes_DF_PCC_3D::test::Test_Types::Stokes_Benchmark_2:
        return std::make_unique<Polydim::examples::Stokes_DF_PCC_3D::test::Stokes_Benchmark_2>();
    default:
        throw std::runtime_error("Test type " + std::to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void solve_stokes(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                  Polydim::examples::Stokes_DF_PCC_3D::Assembler::Stokes_DF_PCC_3D_Problem_Data &assembler_data)
{
    switch (config.SolverType())
    {
    case Program_configuration::Solver_Types::EigenLU: {
        Gedim::Output::PrintGenericMessage("Factorize...", true);
        Gedim::Profiler::StartTime("Factorize");

        Gedim::Eigen_LUSolver solver;
        solver.Initialize(assembler_data.globalMatrixA);

        Gedim::Profiler::StopTime("Factorize");
        Gedim::Output::PrintStatusProgram("Factorize");

        solver.Solve(assembler_data.rightHandSide, assembler_data.solution);
    }
    break;
    case Program_configuration::Solver_Types::EigenBICGSTAB: {
        Gedim::Eigen_BiCGSTABSolver<Eigen::VectorXd, Eigen::SparseMatrix<double>, Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>> solver;

        solver.Initialize(assembler_data.globalMatrixA,
                          assembler_data.rightHandSide,
                          assembler_data.solution,
                          {config.MaxNumberIterations(), config.RelResidualTolerance()});

        const auto solver_data = solver.Solve();

        std::cerr << "EigenBICGSTAB Iterative solver: " << solver_data.Iterations << " " << solver_data.Residual << std::endl;
    }
    break;
    default:
        throw std::runtime_error("not valid solver");
        break;
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D &domain,
                        Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    switch (config.MeshGenerator())
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Tetrahedral:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Polyhedral:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Cubic: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_3D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxVolume(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::CsvImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::VtkImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::OVMImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_3D(meshUtilities, config.MeshGenerator(), config.MeshImportFilePath(), mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator " + std::to_string((unsigned int)config.MeshGenerator()) + " not supported");
    }
}
// ***************************************************************************
Gedim::MeshUtilities::MeshGeometricData3D create_domain_mesh_geometric_properties(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                                                                                  Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_3D_geometry_data(geometryUtilities, mesh);
}
// ***************************************************************************
void export_solution(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                     const Polydim::examples::Stokes_DF_PCC_3D::Assembler::Stokes_DF_PCC_3D_Problem_Data &assembler_data,
                     const Polydim::examples::Stokes_DF_PCC_3D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder)
{
    const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    {
        const char separator = ';';
        std::cout << "ProgramType" << separator;
        std::cout << "VemType" << separator;
        std::cout << "VemOrder" << separator;
        std::cout << "Cell3Ds" << separator;
        std::cout << "Dofs" << separator;
        std::cout << "Strongs" << separator;
        std::cout << "h" << separator;
        std::cout << "errorH1Velocity" << separator;
        std::cout << "errorL2Pressure" << separator;
        std::cout << "normH1Velocity" << separator;
        std::cout << "normL2Pressure" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << VEM_ID << separator;
        std::cout << std::scientific << config.VemOrder() << separator;
        std::cout << std::scientific << mesh.Cell3DTotalNumber() << separator;
        std::cout << std::scientific << count_dofs.num_total_dofs << separator;
        std::cout << std::scientific << count_dofs.num_total_strong << separator;
        std::cout << std::scientific << post_process_data.mesh_size << separator;
        std::cout << std::scientific << post_process_data.error_H1_velocity << separator;
        std::cout << std::scientific << post_process_data.error_L2_pressure << separator;
        std::cout << std::scientific << post_process_data.norm_H1_velocity << separator;
        std::cout << std::scientific << post_process_data.norm_L2_pressure << separator;
        std::cout << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        std::cout << std::scientific << post_process_data.residual_norm << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);

        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "VemType" << separator;
            errorFile << "VemOrder" << separator;
            errorFile << "Cell3Ds" << separator;
            errorFile << "Dofs" << separator;
            errorFile << "Strongs" << separator;
            errorFile << "h" << separator;
            errorFile << "errorH1Velocity" << separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "normH1Velocity" << separator;
            errorFile << "normL2Pressure" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << VEM_ID << separator;
        errorFile << std::scientific << config.VemOrder() << separator;
        errorFile << std::scientific << mesh.Cell3DTotalNumber() << separator;
        errorFile << std::scientific << count_dofs.num_total_dofs << separator;
        errorFile << std::scientific << count_dofs.num_total_strong << separator;
        errorFile << std::scientific << post_process_data.mesh_size << separator;
        errorFile << std::scientific << post_process_data.error_H1_velocity << separator;
        errorFile << std::scientific << post_process_data.error_L2_pressure << separator;
        errorFile << std::scientific << post_process_data.norm_H1_velocity << separator;
        errorFile << std::scientific << post_process_data.norm_L2_pressure << separator;
        errorFile << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << std::scientific << post_process_data.residual_norm << std::endl;

        errorFile.close();
    }

    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolyhedrons(mesh.Cell0DsCoordinates(),
                                    mesh.Cell3DsFacesVertices(),
                                    {{"Numeric Velocity - X",
                                      Gedim::VTPProperty::Formats::Points,
                                      static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[0].size()),
                                      post_process_data.cell0Ds_numeric_velocity[0].data()},
                                     {"Numeric Velocity - Y",
                                      Gedim::VTPProperty::Formats::Points,
                                      static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[1].size()),
                                      post_process_data.cell0Ds_numeric_velocity[1].data()},
                                     {"Numeric Velocity - Z",
                                      Gedim::VTPProperty::Formats::Points,
                                      static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[2].size()),
                                      post_process_data.cell0Ds_numeric_velocity[2].data()},
                                     {"Exact Velocity - X",
                                      Gedim::VTPProperty::Formats::Points,
                                      static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[0].size()),
                                      post_process_data.cell0Ds_exact_velocity[0].data()},
                                     {"Exact Velocity - Y",
                                      Gedim::VTPProperty::Formats::Points,
                                      static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[1].size()),
                                      post_process_data.cell0Ds_exact_velocity[1].data()},
                                     {"Exact Velocity - Z",
                                      Gedim::VTPProperty::Formats::Points,
                                      static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[2].size()),
                                      post_process_data.cell0Ds_exact_velocity[2].data()},
                                     {"ErrorL2Pressure",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(post_process_data.cell3Ds_error_L2_pressure.size()),
                                      post_process_data.cell3Ds_error_L2_pressure.data()},
                                     {"ErrorH1Velocity",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(post_process_data.cell3Ds_error_H1_velocity.size()),
                                      post_process_data.cell3Ds_error_H1_velocity.data()}});

            exporter.Export(exportVtuFolder + "/Solution_" + std::to_string(TEST_ID) + "_" + std::to_string(VEM_ID) +
                            +"_" + std::to_string(config.VemOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_performance(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                        const Assembler::VEM_Performance_Result &performance_data,
                        const std::string &exportFolder)
{
    {
        const char separator = ',';
        std::ofstream exporter;
        const unsigned int Method_ID = static_cast<unsigned int>(config.VemType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
        exporter.open(exportFolder + "/Cell3Ds_MethodPerformance_" + std::to_string(TEST_ID) + "_" +
                      std::to_string(Method_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv");
        exporter.precision(16);

        if (exporter.fail())
            throw std::runtime_error("Error on mesh cell2Ds file");

        exporter << "Cell3D_Index" << separator;
        exporter << "NumQuadPoints_Boundary" << separator;
        exporter << "NumQuadPoints_Internal" << separator;
        exporter << "PiNabla_Cond" << separator;
        exporter << "Pi0k_Cond" << separator;
        exporter << "PiNabla_Error" << separator;
        exporter << "Pi0k_Error" << separator;
        exporter << "GBD_Error" << separator;
        exporter << "HCD_Error" << separator;
        exporter << "Stab_Error" << std::endl;

        for (unsigned int v = 0; v < performance_data.Cell3DsPerformance.size(); v++)
        {
            const auto &cell3DPerformance = performance_data.Cell3DsPerformance[v];

            exporter << std::scientific << v << separator;
            exporter << std::scientific << cell3DPerformance.NumBoundaryQuadraturePoints << separator;
            exporter << std::scientific << cell3DPerformance.NumInternalQuadraturePoints << separator;
            exporter << std::scientific << cell3DPerformance.maxPiNablaConditioning << separator;
            exporter << std::scientific << cell3DPerformance.maxPi0kConditioning << separator;
            exporter << std::scientific << cell3DPerformance.maxErrorPiNabla << separator;
            exporter << std::scientific << cell3DPerformance.maxErrorPi0k << separator;
            exporter << std::scientific << cell3DPerformance.maxErrorGBD << separator;
            exporter << std::scientific << cell3DPerformance.maxErrorHCD << separator;
            exporter << std::scientific << cell3DPerformance.ErrorStabilization << std::endl;
        }

        exporter.close();
    }
}
// ***************************************************************************
void export_velocity_dofs(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &vem_velocity_reference_element_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Polydim::examples::Stokes_DF_PCC_3D::Assembler::Stokes_DF_PCC_3D_Problem_Data &assembler_data,
                          const Polydim::examples::Stokes_DF_PCC_3D::Assembler::PostProcess_Data &post_process_data,
                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    std::list<Eigen::Vector3d> dofs_coordinate;
    std::list<double> dof_dimension_values;
    std::list<double> dof_cell_index_values;
    std::list<std::array<double, 3>> solution_values;
    std::list<std::array<double, 3>> rhs_values;
    std::list<std::array<double, 3>> dof_global_index_values;
    std::list<std::array<double, 3>> dof_type_values;
    std::list<std::array<double, 3>> dof_boundary_type_values;
    std::list<std::array<double, 3>> dof_boundary_marker_values;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
    {
        for (unsigned int loc_i = 0; loc_i < dofs_data[0].CellsDOFs[0].at(c).size(); ++loc_i)
        {
            dofs_coordinate.push_back(mesh.Cell0DCoordinates(c));
            dof_dimension_values.push_back(0);
            dof_cell_index_values.push_back(c);

            std::array<double, 3> sol;
            std::array<double, 3> rhs;
            std::array<double, 3> dof_global;
            std::array<double, 3> dof_type;
            std::array<double, 3> dof_boundary_type;
            std::array<double, 3> dof_boundary_marker;

            for (unsigned int h = 0; h < 3; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(0).at(c);

                const auto &local_dofs = dofs_data[h].CellsDOFs[0].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h];
                    sol[h] = assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h];
                    sol[h] = assembler_data.solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    rhs[h] = assembler_data.rightHandSide.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const std::vector<double> local_edge_coordinates =
            geometryUtilities.EquispaceCoordinates(dofs_data[0].CellsDOFs[1].at(c).size(), 0.0, 1.0, false);
        const Eigen::Vector3d edge_origin = mesh.Cell1DOriginCoordinates(c);
        const Eigen::Vector3d edge_tangent = mesh.Cell1DEndCoordinates(c) - edge_origin;

        for (unsigned int loc_i = 0; loc_i < dofs_data[0].CellsDOFs[1].at(c).size(); ++loc_i)
        {
            std::array<double, 3> sol;
            std::array<double, 3> rhs;
            std::array<double, 3> dof_global;
            std::array<double, 3> dof_type;
            std::array<double, 3> dof_boundary_type;
            std::array<double, 3> dof_boundary_marker;

            dofs_coordinate.push_back(edge_origin + local_edge_coordinates[loc_i] * edge_tangent);
            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(1);

            for (unsigned int h = 0; h < 3; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(c);
                const auto &local_dof = dofs_data[h].CellsDOFs[1].at(c).at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h];
                    sol[h] = assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h];
                    sol[h] = assembler_data.solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    rhs[h] = assembler_data.rightHandSide.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const unsigned int num_loc_dofs = dofs_data[3].CellsDOFs[2].at(c).size();

        if (num_loc_dofs == 0)
            continue;

        const std::vector<double> local_polygon_coordinates =
            geometryUtilities.EquispaceCoordinates(num_loc_dofs + 1, 0.0, 1.0, true);

        unsigned int neigh = mesh.Cell2DNeighbourCell3D(c, 0);
        if (neigh == std::numeric_limits<unsigned int>::max())
            neigh = mesh.Cell2DNeighbourCell3D(c, 1);

        const auto local_index_face = mesh.Cell3DFindFace(neigh, c);

        const Eigen::Vector3d polygon_centroid = mesh_geometric_data.Cell3DsFaces2DCentroids.at(neigh)[local_index_face];
        const auto polygonCentroidEdgesDistance = geometryUtilities.PolygonCentroidEdgesDistance(
            mesh_geometric_data.Cell3DsFaces2DVertices.at(neigh)[local_index_face],
            polygon_centroid,
            mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(neigh)[local_index_face]);

        const double circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            std::array<double, 3> sol;
            std::array<double, 3> rhs;
            std::array<double, 3> dof_global;
            std::array<double, 3> dof_type;
            std::array<double, 3> dof_boundary_type;
            std::array<double, 3> dof_boundary_marker;

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(2);

            if (num_loc_dofs > 1)
                dofs_coordinate.push_back(geometryUtilities.RotatePointsFrom2DTo3D(
                    polygon_centroid +
                        circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                          sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                          0.0),
                    mesh_geometric_data.Cell3DsFacesRotationMatrices.at(neigh)[local_index_face],
                    mesh_geometric_data.Cell3DsFacesTranslations.at(neigh)[local_index_face]));
            else
                dofs_coordinate.push_back(geometryUtilities.RotatePointsFrom2DTo3D(
                    polygon_centroid,
                    mesh_geometric_data.Cell3DsFacesRotationMatrices.at(neigh)[local_index_face],
                    mesh_geometric_data.Cell3DsFacesTranslations.at(neigh)[local_index_face]));

            for (unsigned int h = 0; h < 3; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h + 3].CellsBoundaryInfo.at(2).at(c);
                const auto &local_dof = dofs_data[h + 3].CellsDOFs[2].at(c).at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h + 3];
                    sol[h] = assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h + 3]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h + 3];
                    sol[h] = assembler_data.solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 3]);
                    rhs[h] = assembler_data.rightHandSide.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 3]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
        const unsigned int num_loc_dofs = dofs_data[6].CellsDOFs[3].at(c).size();

        if (num_loc_dofs == 0)
            break;

        const auto local_polyhedron_coordinates = geometryUtilities.fibonacci_sphere(num_loc_dofs);
        const Eigen::Vector3d polyhedron_centroid = mesh_geometric_data.Cell3DsCentroids.at(c);
        const auto polyhedron_centroid_faces_distance =
            geometryUtilities.PolyhedronCentroidFacesDistance(polyhedron_centroid,
                                                              mesh_geometric_data.Cell3DsFacesNormals.at(c),
                                                              mesh_geometric_data.Cell3DsFaces3DVertices.at(c));
        const double polyhedron_in_radius = geometryUtilities.PolyhedronInRadius(polyhedron_centroid_faces_distance);

        const double sphere_diameter = 0.5 * polyhedron_in_radius;

        const auto &local_dofs = dofs_data[6].CellsDOFs[3].at(c);

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(3);
            dofs_coordinate.push_back(polyhedron_centroid + sphere_diameter * local_polyhedron_coordinates.col(loc_i));

            std::array<double, 3> sol = {std::nan(""), std::nan(""), std::nan("")};
            std::array<double, 3> rhs = {std::nan(""), std::nan(""), std::nan("")};
            std::array<double, 3> dof_global = {std::nan(""), std::nan(""), std::nan("")};
            std::array<double, 3> dof_type = {std::nan(""), std::nan(""), std::nan("")};
            std::array<double, 3> dof_boundary_type = {std::nan(""), std::nan(""), std::nan("")};
            std::array<double, 3> dof_boundary_marker = {std::nan(""), std::nan(""), std::nan("")};

            for (unsigned int h = 0; h < 1; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h + 6].CellsBoundaryInfo.at(3).at(c);
                const auto &local_dofs = dofs_data[h + 6].CellsDOFs[3].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h + 6];
                    sol[h] = assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h + 6]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h + 6];
                    sol[h] = assembler_data.solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 6]);
                    rhs[h] = assembler_data.rightHandSide.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 6]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    {
        const unsigned int n = dofs_coordinate.size();
        Eigen::MatrixXd coordinates(3, dofs_coordinate.size());
        unsigned int c = 0;
        for (const auto &dof_coordinate : dofs_coordinate)
            coordinates.col(c++) << dof_coordinate;
        const auto dof_dimension_values_data = std::vector<double>(dof_dimension_values.begin(), dof_dimension_values.end());
        const auto dof_cell_index_values_data = std::vector<double>(dof_cell_index_values.begin(), dof_cell_index_values.end());
        Eigen::VectorXd rhs_values_data(3 * n);
        c = 0;
        for (const auto &v : rhs_values)
        {
            for (unsigned int d = 0; d < 3; d++)
                rhs_values_data[3 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd solution_values_data(3 * n);
        c = 0;
        for (const auto &v : solution_values)
        {
            for (unsigned int d = 0; d < 3; d++)
                solution_values_data[3 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_global_index_values_data(3 * n);
        c = 0;
        for (const auto &v : dof_global_index_values)
        {
            for (unsigned int d = 0; d < 3; d++)
                dof_global_index_values_data[3 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_type_values_data(3 * n);
        c = 0;
        for (const auto &v : dof_type_values)
        {
            for (unsigned int d = 0; d < 3; d++)
                dof_type_values_data[3 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_boundary_type_values_data(3 * n);
        c = 0;
        for (const auto &v : dof_boundary_type_values)
        {
            for (unsigned int d = 0; d < 3; d++)
                dof_boundary_type_values_data[3 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_boundary_marker_values_data(3 * n);
        c = 0;
        for (const auto &v : dof_boundary_marker_values)
        {
            for (unsigned int d = 0; d < 3; d++)
                dof_boundary_marker_values_data[3 * c + d] = v[d];
            c++;
        }

        Gedim::VTKUtilities exporter;
        exporter.AddPoints(coordinates,
                           {{"cell_dimension",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_dimension_values_data.size()),
                             dof_dimension_values_data.data()},
                            {"cell_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_cell_index_values_data.size()),
                             dof_cell_index_values_data.data()},
                            {"boundary_type",
                             Gedim::VTPProperty::Formats::PointsArray,
                             static_cast<unsigned int>(dof_boundary_type_values_data.size()),
                             dof_boundary_type_values_data.data()},
                            {"boundary_marker",
                             Gedim::VTPProperty::Formats::PointsArray,
                             static_cast<unsigned int>(dof_boundary_marker_values_data.size()),
                             dof_boundary_marker_values_data.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::PointsArray,
                             static_cast<unsigned int>(dof_global_index_values_data.size()),
                             dof_global_index_values_data.data()},
                            {"dof_type",
                             Gedim::VTPProperty::Formats::PointsArray,
                             static_cast<unsigned int>(dof_type_values_data.size()),
                             dof_type_values_data.data()},
                            {"rhs",
                             Gedim::VTPProperty::Formats::PointsArray,
                             static_cast<unsigned int>(rhs_values_data.size()),
                             rhs_values_data.data()},
                            {"solution",
                             Gedim::VTPProperty::Formats::PointsArray,
                             static_cast<unsigned int>(solution_values_data.size()),
                             solution_values_data.data()}});

        const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

        exporter.Export(exportVtuFolder + "/dofs_" + std::to_string(TEST_ID) + "_" + std::to_string(VEM_ID) + +"_" +
                        std::to_string(config.VemOrder()) + ".vtu");
    }
}
// ***************************************************************************
void export_discrepancy_errors(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Polydim::examples::Stokes_DF_PCC_3D::Assembler::DiscrepancyErrors_Data &discrepancy_errors_data,
                               const std::string &exportSolutionFolder,
                               const std::string &exportVtuFolder)
{
    const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    {
        const char separator = ';';
        std::cout << "ProgramType" << separator;
        std::cout << "VemType" << separator;
        std::cout << "VemOrder" << separator;
        std::cout << "Cell3Ds" << separator;
        std::cout << "VelocityDofsRatio" << separator;
        std::cout << "PressureDofsRatio" << separator;
        std::cout << "discrepancyErrorH1Velocity" << separator;
        std::cout << "discrepancyErrorL2Pressure" << separator;
        std::cout << "normH1FULLVelocity" << separator;
        std::cout << "normL2FULLProjectedPressure" << separator;
        std::cout << "reducedResidual" << separator;
        std::cout << "fullResidual" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << VEM_ID << separator;
        std::cout << std::scientific << config.VemOrder() << separator;
        std::cout << std::scientific << mesh.Cell3DTotalNumber() << separator;
        std::cout << std::scientific << discrepancy_errors_data.velocity_dofs_ratio << separator;
        std::cout << std::scientific << discrepancy_errors_data.pressure_dofs_ratio << separator;
        std::cout << std::scientific << discrepancy_errors_data.discrepancy_error_H1_velocity << separator;
        std::cout << std::scientific << discrepancy_errors_data.discrepancy_error_L2_pressure << separator;
        std::cout << std::scientific << discrepancy_errors_data.full_norm_H1_velocity << separator;
        std::cout << std::scientific << discrepancy_errors_data.full_norm_L2_pressure << separator;
        std::cout << std::scientific << discrepancy_errors_data.reduced_residual_norm << separator;
        std::cout << std::scientific << discrepancy_errors_data.residual_norm << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/DiscrepancyErrors_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);

        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "VemType" << separator;
            errorFile << "VemOrder" << separator;
            errorFile << "Cell3Ds" << separator;
            errorFile << "VelocityDofsRatio" << separator;
            errorFile << "PressureDofsRatio" << separator;
            errorFile << "discrepancyErrorH1Velocity" << separator;
            errorFile << "discrepancyErrorL2Pressure" << separator;
            errorFile << "normH1FULLVelocity" << separator;
            errorFile << "normL2FULLProjectedPressure" << separator;
            errorFile << "reducedResidual" << separator;
            errorFile << "fullResidual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << VEM_ID << separator;
        errorFile << std::scientific << config.VemOrder() << separator;
        errorFile << std::scientific << mesh.Cell3DTotalNumber() << separator;
        errorFile << std::scientific << discrepancy_errors_data.velocity_dofs_ratio << separator;
        errorFile << std::scientific << discrepancy_errors_data.pressure_dofs_ratio << separator;
        errorFile << std::scientific << discrepancy_errors_data.discrepancy_error_H1_velocity << separator;
        errorFile << std::scientific << discrepancy_errors_data.discrepancy_error_L2_pressure << separator;
        errorFile << std::scientific << discrepancy_errors_data.full_norm_H1_velocity << separator;
        errorFile << std::scientific << discrepancy_errors_data.full_norm_L2_pressure << separator;
        errorFile << std::scientific << discrepancy_errors_data.reduced_residual_norm << separator;
        errorFile << std::scientific << discrepancy_errors_data.residual_norm << std::endl;

        errorFile.close();
    }

    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(
                mesh.Cell0DsCoordinates(),
                mesh.Cell2DsVertices(),
                {{"ErrorL2Pressure",
                  Gedim::VTPProperty::Formats::Cells,
                  static_cast<unsigned int>(discrepancy_errors_data.cell3Ds_discrepancy_error_L2_pressure.size()),
                  discrepancy_errors_data.cell3Ds_discrepancy_error_L2_pressure.data()},
                 {"ErrorH1Velocity",
                  Gedim::VTPProperty::Formats::Cells,
                  static_cast<unsigned int>(discrepancy_errors_data.cell3Ds_discrepancy_error_H1_velocity.size()),
                  discrepancy_errors_data.cell3Ds_discrepancy_error_H1_velocity.data()}});

            exporter.Export(exportVtuFolder + "/DiscrepancyErrors_" + std::to_string(TEST_ID) + "_" +
                            std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Stokes_DF_PCC_3D
} // namespace examples
} // namespace Polydim
