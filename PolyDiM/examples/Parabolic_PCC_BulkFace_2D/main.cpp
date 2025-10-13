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

#include "Eigen_LUSolver.hpp"
#include "FEM_PCC_1D_Creator.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "VTKUtilities.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

int main(int argc, char **argv)
{
    Polydim::examples::Parabolic_PCC_BulkFace_2D::Program_configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini", false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    /// Create folders
    const std::string exportFolder = config.ExportFolder();
    Gedim::Output::CreateFolder(exportFolder);

    const std::string exportCsvFolder = exportFolder + "/Mesh";
    Gedim::Output::CreateFolder(exportCsvFolder);
    const std::string exportVtuFolder = exportFolder + "/Paraview";
    Gedim::Output::CreateFolder(exportVtuFolder);
    const std::string exportSolutionFolder = exportFolder + "/Solution";
    Gedim::Output::CreateFolder(exportSolutionFolder);

    const std::string logFolder = exportFolder + "/Log";

    /// Set Profiler
    Gedim::Profiler::ActivateProfiler = true;

    /// Set Log folder
    Gedim::Output::CreateFolder(logFolder);
    Gedim::LogFile::LogFolder = logFolder;

    /// Export Configuration of the following Run
    Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini", false);

    /// Set problem
    Gedim::Output::PrintGenericMessage("SetProblem...", true);
    Gedim::Profiler::StartTime("SetProblem");

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    const auto test = Polydim::examples::Parabolic_PCC_BulkFace_2D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(domain.spatial_domain.vertices);
        vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }

    Gedim::Profiler::StopTime("SetProblem");
    Gedim::Output::PrintStatusProgram("SetProblem");

    /// Create domain mesh
    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData_2D;
    Gedim::MeshMatricesDAO mesh_2D(meshData_2D);

    Polydim::examples::Parabolic_PCC_BulkFace_2D::program_utilities::create_domain_mesh_2D(config, domain.spatial_domain, mesh_2D);

    Gedim::MeshMatrices meshData_1D;
    Gedim::MeshMatricesDAO mesh_1D(meshData_1D);

    std::vector<unsigned int> cell0DsFilter;
    std::vector<unsigned int> cell1DsFilter;

    for (unsigned int c = 0; c < mesh_2D.Cell0DTotalNumber(); c++)
        if (mesh_2D.Cell0DMarker(c) != 0)
            cell0DsFilter.push_back(c);

    for (unsigned int c = 0; c < mesh_2D.Cell1DTotalNumber(); c++)
        if (mesh_2D.Cell1DMarker(c) != 0)
            cell1DsFilter.push_back(c);

    Gedim::MeshUtilities::ExtractMeshData extract_data =
        meshUtilities.ExtractMesh1D(cell0DsFilter, cell1DsFilter, mesh_2D, mesh_1D);

    const auto time_steps =
        Polydim::examples::Parabolic_PCC_BulkFace_2D::program_utilities::create_time_steps(config, domain.time_domain);
    const double delta_time = config.TimeStep();

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    // Export the domain mesh
    {
        meshUtilities.ExportMeshToVTU(mesh_2D, exportVtuFolder, "Domain_Mesh_2D");
        meshUtilities.ExportMeshToVTU(mesh_1D, exportVtuFolder, "Domain_Mesh_1D");
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto mesh_geometric_data_2D =
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometryUtilities, meshUtilities, mesh_2D);

    const auto mesh_geometric_data_1D =
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_1D_geometry_data(geometryUtilities, meshUtilities, mesh_1D);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    const auto reference_element_data_2D =
        Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    const auto reference_element_1D = Polydim::FEM::PCC::create_FEM_PCC_1D_reference_element(
        Polydim::FEM::PCC::FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace);

    Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data reference_element_data_1D;

    if (config.MethodType() == Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes::FEM_PCC)
        reference_element_data_1D =
            reference_element_1D->Create(config.MethodOrder(), Polydim::FEM::PCC::FEM_PCC_1D_Types::Equispaced);
    else
        reference_element_data_1D =
            reference_element_1D->Create(config.MethodOrder(), Polydim::FEM::PCC::FEM_PCC_1D_Types::GaussLobatto);

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data_2D(mesh_2D);
    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data_1D(mesh_1D);

    Polydim::PDETools::DOFs::DOFsManager dofManager;

    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> mesh_dofs_info(2);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(2);

    mesh_dofs_info[0] = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(reference_element_data_2D, mesh_2D, boundary_info);
    dofs_data[0] = dofManager.CreateDOFs_2D(mesh_dofs_info[0], mesh_connectivity_data_2D);

    mesh_dofs_info[1] = dofManager.Create_Constant_DOFsInfo_1D(
        mesh_connectivity_data_1D,
        {{reference_element_data_1D.NumDofs0D, reference_element_data_1D.NumDofs1D, 0, 0}, boundary_info});

    dofs_data[1] = dofManager.CreateDOFs_1D(mesh_dofs_info[1], mesh_connectivity_data_1D);

    const auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(dofs_data);

    Gedim::Output::PrintGenericMessage("2D Discrete Space with " + std::to_string(dofs_data[0].NumberDOFs) +
                                           " DOFs and " + std::to_string(dofs_data[0].NumberStrongs) + " STRONGs",
                                       true);

    Gedim::Output::PrintGenericMessage("1D Discrete Space with " + std::to_string(dofs_data[1].NumberDOFs) +
                                           " DOFs and " + std::to_string(dofs_data[1].NumberStrongs) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem Discrete Type " +
                                           std::to_string(static_cast<unsigned int>(config.MethodType())) + "...",
                                       true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Parabolic_PCC_BulkFace_2D::Assembler assembler;
    Polydim::examples::Parabolic_PCC_BulkFace_2D::Assembler::Elliptic_PCC_BF_2D_Problem_Data assembler_data;

    assembler.ComputeInitialCondition(config,
                                      mesh_2D,
                                      mesh_geometric_data_2D,
                                      mesh_1D,
                                      mesh_geometric_data_1D,
                                      dofs_data,
                                      count_dofs,
                                      reference_element_data_2D,
                                      reference_element_data_1D,
                                      *test,
                                      assembler_data.initial_solution);

    assembler.AssembleMatrix(config,
                             mesh_2D,
                             mesh_geometric_data_2D,
                             mesh_1D,
                             mesh_geometric_data_1D,
                             extract_data,
                             mesh_dofs_info,
                             dofs_data,
                             count_dofs,
                             reference_element_data_2D,
                             reference_element_data_1D,
                             *test,
                             assembler_data.globalMatrixA,
                             assembler_data.globalMatrixM);

    assembler_data.globalMatrixA *= delta_time;
    assembler_data.globalMatrixA += assembler_data.globalMatrixM;

    Gedim::Profiler::StopTime("AssembleSystem");
    Gedim::Output::PrintStatusProgram("AssembleSystem");

    Gedim::Output::PrintGenericMessage("Factorize...", true);
    Gedim::Profiler::StartTime("Factorize");

    Gedim::Eigen_LUSolver solver;
    solver.Initialize(assembler_data.globalMatrixA);

    Gedim::Profiler::StopTime("Factorize");
    Gedim::Output::PrintStatusProgram("Factorize");

    for (unsigned int t = 1; t < time_steps.size(); t++)
    {
        const double value_time = time_steps[t];

        assembler.AssembleRhs(config,
                              value_time,
                              mesh_2D,
                              mesh_geometric_data_2D,
                              mesh_1D,
                              mesh_geometric_data_1D,
                              extract_data,
                              mesh_dofs_info,
                              dofs_data,
                              count_dofs,
                              reference_element_data_2D,
                              reference_element_data_1D,
                              *test,
                              assembler_data.rightHandSide);

        assembler_data.rightHandSide *= delta_time;
        assembler_data.rightHandSide.SumMultiplication(assembler_data.globalMatrixM, assembler_data.initial_solution);

        Gedim::Output::PrintGenericMessage("Solve...", true);
        Gedim::Profiler::StartTime("Solve");

        assembler_data.solution.SetSize(count_dofs.num_total_dofs);
        solver.Solve(assembler_data.rightHandSide, assembler_data.solution);

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");

        Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
        Gedim::Profiler::StartTime("ComputeErrors");

        auto post_process_data = assembler.PostProcessSolution(config,
                                                               value_time,
                                                               mesh_2D,
                                                               mesh_geometric_data_2D,
                                                               mesh_1D,
                                                               mesh_geometric_data_1D,
                                                               mesh_dofs_info,
                                                               dofs_data,
                                                               count_dofs,
                                                               reference_element_data_2D,
                                                               reference_element_data_1D,
                                                               *test,
                                                               assembler_data);

        Gedim::Profiler::StopTime("ComputeErrors");
        Gedim::Output::PrintStatusProgram("ComputeErrors");

        Gedim::Output::PrintGenericMessage("ExportSolution...", true);
        Gedim::Profiler::StartTime("ExportSolution");

        Polydim::examples::Parabolic_PCC_BulkFace_2D::program_utilities::export_solution(config,
                                                                                         value_time,
                                                                                         mesh_2D,
                                                                                         mesh_1D,
                                                                                         dofs_data,
                                                                                         assembler_data,
                                                                                         post_process_data,
                                                                                         exportSolutionFolder,
                                                                                         exportVtuFolder);

        Gedim::Profiler::StopTime("ExportSolution");
        Gedim::Output::PrintStatusProgram("ExportSolution");

        assembler_data.initial_solution.Copy(assembler_data.solution);
    }

    Gedim::Output::PrintGenericMessage("ComputeMethodPerformance...", true);
    Gedim::Profiler::StartTime("ComputeMethodPerformance");

    if (config.ComputeMethodPerformance())
    {
        const auto performance = assembler.ComputePerformance_2D(config, mesh_2D, mesh_geometric_data_2D, reference_element_data_2D);

        Polydim::examples::Parabolic_PCC_BulkFace_2D::program_utilities::export_performance_2D(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeMethodPerformance");
    Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");

    return 0;
}
