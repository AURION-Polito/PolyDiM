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

#include "DOFsManager.hpp"
#include "Eigen_LUSolver.hpp"
#include "LocalSpace_MCC_2D.hpp"
#include "MeshDAOExporterToCsv.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

int main(int argc, char **argv)
{
    Polydim::examples::Elliptic_MCC_2D::Program_configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini", false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    /// Create folders
    const std::string exportFolder = config.ExportFolder();
    Gedim::Output::CreateFolder(exportFolder);

    const std::string exportCsvFolder = exportFolder + "/Csv";
    Gedim::Output::CreateFolder(exportCsvFolder);
    const std::string exportVtuFolder = exportFolder + "/Paraview";
    Gedim::Output::CreateFolder(exportVtuFolder);

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

    const auto test = Polydim::examples::Elliptic_MCC_2D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
    if (config.ExportFormat()[1])
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(domain.vertices);
        vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
        Gedim::Output::PrintGenericMessage(Gedim::Output::MagentaColor + "Domain is exported in: " + exportVtuFolder +
                                               "/Domain.vtu" + Gedim::Output::EndColor,
                                           true);
    }

    Gedim::Profiler::StopTime("SetProblem");
    Gedim::Output::PrintStatusProgram("SetProblem");

    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Elliptic_MCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    // Export the domain mesh
    if (config.ExportFormat()[1])
    {
        Gedim::MeshUtilities meshUtilities;
        meshUtilities.ExportMeshToVTU(mesh, exportVtuFolder, "Domain_Mesh");
        Gedim::Output::PrintGenericMessage(Gedim::Output::MagentaColor + "Mesh is exported in: " + exportVtuFolder +
                                               Gedim::Output::EndColor,
                                           true);
    }

    if (config.ExportFormat()[0])
    {
        const std::string exportMeshFolder = exportCsvFolder + "/Mesh";
        Gedim::Output::CreateFolder(exportMeshFolder);

        const Gedim::MeshFromCsvUtilities csv_utilities;
        Gedim::MeshFromCsvUtilities::Configuration csv_configuration;
        csv_configuration.Folder = exportMeshFolder;
        Gedim::MeshDAOExporterToCsv exporter_to_csv(csv_utilities);
        exporter_to_csv.Export(csv_configuration, mesh);
        Gedim::Output::PrintGenericMessage(Gedim::Output::MagentaColor + "Mesh is exported in: " + exportMeshFolder +
                                               Gedim::Output::EndColor,
                                           true);
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto meshGeometricData =
        Polydim::examples::Elliptic_MCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space

    Gedim::Output::PrintGenericMessage("CreateSpace...", true);
    Gedim::Profiler::StartTime("CreateSpace");

    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_MCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);

    const auto reference_element_num_dofs = Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElementNumDOFs(reference_element_data);

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(2);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(2);

    for (unsigned int h = 0; h < 2; h++)
    {
        meshDOFsInfo[h] =
            dofManager.Create_Constant_DOFsInfo_2D(mesh_connectivity_data, {reference_element_num_dofs[h], boundary_info});

        dofs_data[h] = dofManager.CreateDOFs_2D(meshDOFsInfo[h], mesh_connectivity_data);
    }

    const auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(dofs_data);

    Gedim::Profiler::StopTime("CreateSpace");
    Gedim::Output::PrintStatusProgram("CreateSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Elliptic_MCC_2D::Assembler assembler;
    auto assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, count_dofs, reference_element_data, *test);

    Gedim::Profiler::StopTime("AssembleSystem");
    Gedim::Output::PrintStatusProgram("AssembleSystem");

    if (count_dofs.num_total_dofs > 0)
    {
        Gedim::Output::PrintGenericMessage("Factorize...", true);
        Gedim::Profiler::StartTime("Factorize");

        Gedim::Eigen_LUSolver solver;
        solver.Initialize(assembler_data.globalMatrixA);

        Gedim::Profiler::StopTime("Factorize");
        Gedim::Output::PrintStatusProgram("Factorize");

        Gedim::Output::PrintGenericMessage("Solve...", true);
        Gedim::Profiler::StartTime("Solve");

        solver.Solve(assembler_data.rightHandSide, assembler_data.solution);

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");
    }

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    auto post_process_data =
        assembler.PostProcessSolution(config, mesh, meshGeometricData, dofs_data, count_dofs, reference_element_data, assembler_data, *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolutionAndErrors...", true);
    Gedim::Profiler::StartTime("ExportSolutionAndErrors");

    Polydim::examples::Elliptic_MCC_2D::program_utilities::export_solution(config, mesh, dofs_data, assembler_data, post_process_data, exportCsvFolder, exportVtuFolder);
    Polydim::examples::Elliptic_MCC_2D::program_utilities::export_velocity_dofs(config,
                                                                                mesh,
                                                                                meshGeometricData,
                                                                                meshDOFsInfo[0],
                                                                                dofs_data[0],
                                                                                assembler_data,
                                                                                post_process_data,
                                                                                exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolutionAndErrors");
    Gedim::Output::PrintStatusProgram("ExportSolutionAndErrors");

    if (config.ComputeMethodPerformance())
    {
        Gedim::Output::PrintGenericMessage("ComputeMethodPerformance...", true);
        Gedim::Profiler::StartTime("ComputeMethodPerformance");

        const auto performance_data = assembler.ComputePerformance(config, mesh, meshGeometricData, reference_element_data);
        Polydim::examples::Elliptic_MCC_2D::program_utilities::export_performance(config, performance_data, exportCsvFolder);

        Gedim::Profiler::StopTime("ComputeMethodPerformance");
        Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");
    }

    return 0;
}
