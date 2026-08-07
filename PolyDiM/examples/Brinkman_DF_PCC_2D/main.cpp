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

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "Eigen_LUSolver.hpp"
#include "LocalSpace_DF_PCC_2D.hpp"
#include "MeshDAOExporterToCsv.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"

int main(int argc, char **argv)
{
    Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration config;

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

    const auto test = Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
    if (config.ExportFormat()[1])
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(domain.vertices);
        vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }

    Gedim::Profiler::StopTime("SetProblem");
    Gedim::Output::PrintStatusProgram("SetProblem");

    /// Create domain mesh
    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

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
        Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);

    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_DF_PCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    Polydim::PDETools::DOFs::DOFsManager dofManager;

    const auto mesh_dofs_info = Polydim::PDETools::LocalSpace_DF_PCC_2D::SetMeshDOFsInfo(reference_element_data, mesh, boundary_info);
    const unsigned int num_mesh_dofs_info = mesh_dofs_info.size();
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(num_mesh_dofs_info);

    for (unsigned int i = 0; i < num_mesh_dofs_info; i++)
        dofs_data[i] = dofManager.CreateDOFs_2D(mesh_dofs_info[i], mesh_connectivity_data);

    auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(dofs_data);

    if (count_dofs.num_total_boundary_dofs == 0)
        count_dofs.num_total_dofs += 1; // lagrange

    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace...", true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Brinkman_DF_PCC_2D::Assembler assembler;
    auto assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, mesh_dofs_info, dofs_data, count_dofs, reference_element_data, *test);

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

    Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::export_solution(config,
                                                                              mesh,
                                                                              dofs_data,
                                                                              count_dofs,
                                                                              assembler_data,
                                                                              post_process_data,
                                                                              exportCsvFolder,
                                                                              exportVtuFolder);

    Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::export_velocity_dofs(config,
                                                                                   mesh,
                                                                                   meshGeometricData,
                                                                                   mesh_dofs_info,
                                                                                   reference_element_data,
                                                                                   dofs_data,
                                                                                   count_dofs,
                                                                                   assembler_data,
                                                                                   exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolutionAndErrors");
    Gedim::Output::PrintStatusProgram("ExportSolutionAndErrors");

    if (config.ComputeMethodPerformance())
    {
        Gedim::Output::PrintGenericMessage("ComputeMethodPerformance...", true);
        Gedim::Profiler::StartTime("ComputeMethodPerformance");

        const auto performance = assembler.ComputeMethodPerformance(config, mesh, meshGeometricData, reference_element_data);
        Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::export_performance(config, performance, exportCsvFolder);

        Gedim::Profiler::StopTime("ComputeMethodPerformance");
        Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");
    }

    if (config.ComputeDiscrepancyError())
    {
        if (config.MethodType() != Polydim::PDETools::LocalSpace_DF_PCC_2D::MethodTypes::VEM_DF_PCC_REDUCED)
            throw std::runtime_error("not valid method type");

        Gedim::Output::PrintGenericMessage("ComputeDiscrepancyErrors...", true);
        Gedim::Profiler::StartTime("ComputeDiscrepancyErrors");

        const auto discrepancy_errors_data =
            assembler.ComputeDiscrepancyErrors(config, mesh, meshGeometricData, dofs_data, count_dofs, reference_element_data, assembler_data, *test);

        Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::export_discrepancy_errors(config,
                                                                                            mesh,
                                                                                            discrepancy_errors_data,
                                                                                            exportCsvFolder,
                                                                                            exportVtuFolder);

        Gedim::Profiler::StopTime("ComputeDiscrepancyErrors");
        Gedim::Output::PrintStatusProgram("ComputeDiscrepancyErrors");
    }

    return 0;
}
