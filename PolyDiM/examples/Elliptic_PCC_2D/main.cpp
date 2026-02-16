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
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "VTKUtilities.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

unsigned int Polydim::examples::Elliptic_PCC_2D::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::Elliptic_PCC_2D::Program_configuration config;

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

    Polydim::examples::Elliptic_PCC_2D::test::Patch_Test::order = config.MethodOrder();

    const auto test = Polydim::examples::Elliptic_PCC_2D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
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

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Elliptic_PCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    // Export the domain mesh
    {
        Gedim::MeshUtilities meshUtilities;
        meshUtilities.ExportMeshToVTU(mesh, exportVtuFolder, "Domain_Mesh");
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    const auto meshGeometricData = Polydim::examples::Elliptic_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(
        config,
        Polydim::PDETools::LocalSpace_PCC_2D::MeshGeometricDataConfigiguration(reference_element_data),
        mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);

    Polydim::PDETools::DOFs::DOFsManager dofManager;

    const auto meshDOFsInfo = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(reference_element_data, mesh, boundary_info);
    const auto dofs_data = dofManager.CreateDOFs_2D(meshDOFsInfo, mesh_connectivity_data);

    Gedim::Output::PrintGenericMessage("Discrete Space with " + std::to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                           std::to_string(dofs_data.NumberStrongs) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem Discrete Type " +
                                           std::to_string(static_cast<unsigned int>(config.MethodType())) + "...",
                                       true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Elliptic_PCC_2D::Assembler assembler;
    auto assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, reference_element_data, *test);

    Gedim::Profiler::StopTime("AssembleSystem");
    Gedim::Output::PrintStatusProgram("AssembleSystem");

    if (dofs_data.NumberDOFs > 0)
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

        {
          using namespace Gedim;
        std::cout.precision(2);
        std::cout<< std::scientific<< "A: "<< assembler_data.globalMatrixA<< std::endl;
        std::cout<< std::scientific<< "A_D: "<< assembler_data.dirichletMatrixA<< std::endl;
        std::cout << std::scientific << "r: " << assembler_data.rightHandSide << std::endl;
        std::cout << std::scientific << "u: " << assembler_data.solution << std::endl;
        std::cout << std::scientific << "u_D: " << assembler_data.solutionDirichlet << std::endl;
        }

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");
    }

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    auto post_process_data =
        assembler.PostProcessSolution(config, mesh, meshGeometricData, dofs_data, reference_element_data, assembler_data, *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_PCC_2D::program_utilities::export_solution(config, mesh, dofs_data, assembler_data, post_process_data, exportSolutionFolder, exportVtuFolder);

    Polydim::examples::Elliptic_PCC_2D::program_utilities::export_dofs(config,
                                                                       mesh,
                                                                       meshGeometricData,
                                                                       meshDOFsInfo,
                                                                       dofs_data,
                                                                       reference_element_data,
                                                                       assembler_data,
                                                                       post_process_data,
                                                                       exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeMethodPerformance...", true);
    Gedim::Profiler::StartTime("ComputeMethodPerformance");

    if (config.ComputeMethodPerformance())
    {
        const auto performance = assembler.ComputePerformance(config, mesh, meshGeometricData, reference_element_data);

        Polydim::examples::Elliptic_PCC_2D::program_utilities::export_performance(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeMethodPerformance");
    Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");

    return 0;
}
