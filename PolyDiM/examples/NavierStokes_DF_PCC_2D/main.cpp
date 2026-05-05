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
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"

int main(int argc, char **argv)
{
    Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration config;

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

    const auto test = Polydim::examples::NavierStokes_DF_PCC_2D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
    if (domain.vertices.size() > 0)
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

    Polydim::examples::NavierStokes_DF_PCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    // Export the domain mesh
    {
        Gedim::MeshUtilities meshUtilities;
        meshUtilities.ExportMeshToVTU(mesh, exportVtuFolder, "Domain_Mesh");
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto meshGeometricData =
        Polydim::examples::NavierStokes_DF_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
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

    Gedim::Output::PrintGenericMessage("Discrete Space with " + std::to_string(count_dofs.num_total_dofs) +
                                           " DOFs and " + std::to_string(count_dofs.num_total_strong) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem Discrete Type " + std::to_string((unsigned int)config.MethodType()) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::NavierStokes_DF_PCC_2D::Assembler assembler;
    auto assembler_data =
        assembler.AssembleStokes(config, mesh, meshGeometricData, mesh_dofs_info, dofs_data, count_dofs, reference_element_data, *test);

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

        solver.Solve(assembler_data.rightHandSide, assembler_data.previousIteration);

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");
    }

    Gedim::Eigen_Array<> residual;
    residual.SetSize(count_dofs.num_total_dofs);
    residual.SumMultiplication(assembler_data.globalMatrixA, assembler_data.previousIteration);
    residual -= assembler_data.rightHandSide;

    double residual_norm;
    const double initial_residual_norm = residual.Norm();
    const double initial_solution_norm = assembler_data.previousIteration.Norm();

    const unsigned int NLMAXNumIterations = config.NLMaxNumberIterations();
    unsigned int num_nl_iterations;
    for (unsigned int l = 0; l < NLMAXNumIterations; l++)
    {
        Gedim::Output::PrintGenericMessage("AssembleNavierStokesTerm Method Type " +
                                               std::to_string((unsigned int)config.MethodType()) + "...",
                                           true);
        Gedim::Profiler::StartTime("AssembleNavierStokesTerm");

        assembler.AssembleNavierStokes(config, mesh, meshGeometricData, mesh_dofs_info, dofs_data, count_dofs, reference_element_data, assembler_data);

        Gedim::Profiler::StopTime("AssembleNavierStokesTerm");
        Gedim::Output::PrintStatusProgram("AssembleNavierStokesTerm");

        if (count_dofs.num_total_dofs > 0)
        {
            Gedim::Output::PrintGenericMessage("Factorize...", true);
            Gedim::Profiler::StartTime("Factorize");

            Gedim::Eigen_LUSolver solver;
            solver.Initialize(assembler_data.globalMatrixC);

            Gedim::Profiler::StopTime("Factorize");
            Gedim::Output::PrintStatusProgram("Factorize");

            Gedim::Output::PrintGenericMessage("Solve...", true);
            Gedim::Profiler::StartTime("Solve");

            assembler_data.solution.Zeros();
            solver.Solve(assembler_data.rightHandSideC, assembler_data.solution);

            Gedim::Profiler::StopTime("Solve");
            Gedim::Output::PrintStatusProgram("Solve");
        }

        residual.Zeros();
        residual.SumMultiplication(assembler_data.globalMatrixC, assembler_data.solution);
        residual -= assembler_data.rightHandSideC;
        residual_norm = residual.Norm();

        Gedim::Eigen_Array<> delta;
        delta.Copy(assembler_data.solution);
        delta -= assembler_data.previousIteration;
        const double delta_norm = delta.Norm();

        assembler_data.previousIteration.Copy(assembler_data.solution);

        std::cout << "NL It: " << l + 1 << ", Residual Norm: " << residual_norm << ", Delta norm: " << delta_norm << std::endl;

        if ((residual_norm <= config.NLRelResidualTolerance() * initial_residual_norm + config.NLAbsResidualTolerance()) &&
            (delta_norm <= config.NLRelChangeInSolutionTolerance() * initial_solution_norm + config.NLAbsChangeInSolutionTolerance()))
        {
            num_nl_iterations = l + 1;
            break;
        }
    }

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    auto post_process_data =
        assembler.PostProcessSolution(config, mesh, meshGeometricData, dofs_data, count_dofs, reference_element_data, assembler_data, residual_norm, *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::NavierStokes_DF_PCC_2D::program_utilities::export_solution(config,
                                                                                  mesh,
                                                                                  dofs_data,
                                                                                  count_dofs,
                                                                                  assembler_data,
                                                                                  post_process_data,
                                                                                  num_nl_iterations,
                                                                                  exportSolutionFolder,
                                                                                  exportVtuFolder);

    Polydim::examples::NavierStokes_DF_PCC_2D::program_utilities::export_velocity_dofs(config,
                                                                                       mesh,
                                                                                       meshGeometricData,
                                                                                       mesh_dofs_info,
                                                                                       reference_element_data,
                                                                                       dofs_data,
                                                                                       count_dofs,
                                                                                       assembler_data,
                                                                                       post_process_data,
                                                                                       exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
    Gedim::Profiler::StartTime("ComputeVEMPerformance");

    if (config.ComputeMethodPerformance())
    {
        const auto performance = assembler.ComputeMethodPerformance(config, mesh, meshGeometricData, reference_element_data);

        Polydim::examples::NavierStokes_DF_PCC_2D::program_utilities::export_performance(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
