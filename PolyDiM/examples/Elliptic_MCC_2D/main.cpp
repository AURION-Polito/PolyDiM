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
#include "MeshMatricesDAO.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

unsigned int Polydim::examples::Elliptic_MCC_2D::test::Patch_Test::order;

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

    Polydim::examples::Elliptic_MCC_2D::test::Patch_Test::order = config.MethodOrder();

    const auto test = Polydim::examples::Elliptic_MCC_2D::program_utilities::create_test(config);

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

    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Elliptic_MCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

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
        Polydim::examples::Elliptic_MCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

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

    Gedim::Output::PrintGenericMessage("VEM Space with " + std::to_string(count_dofs.num_total_dofs) + " DOFs and " +
                                           std::to_string(count_dofs.num_total_strong) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " +
                                           std::to_string(static_cast<unsigned int>(config.MethodType())) + "...",
                                       true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Elliptic_MCC_2D::Assembler assembler;
    auto assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, count_dofs, reference_element_data, *test);

    std::cout << assembler_data.solutionNeumann << std::endl;
    std::cout << assembler_data.rightHandSide << std::endl;

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

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_MCC_2D::program_utilities::export_solution(config, mesh, dofs_data, assembler_data, post_process_data, exportSolutionFolder, exportVtuFolder);
    Polydim::examples::Elliptic_MCC_2D::program_utilities::export_velocity_dofs(config,
                                                                                mesh,
                                                                                meshGeometricData,
                                                                                meshDOFsInfo[0],
                                                                                dofs_data[0],
                                                                                assembler_data,
                                                                                post_process_data,
                                                                                exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
    Gedim::Profiler::StartTime("ComputeVEMPerformance");

    if (config.ComputeMethodPerformance())
    {
        const auto vemPerformance = assembler.ComputePerformance(config, mesh, meshGeometricData, reference_element_data);
        {
            const char separator = ',';
            /// Export Cell2Ds VEM performance
            std::ofstream exporter;

            exporter.open(exportSolutionFolder + "/Cell2Ds_VEMPerformance.csv");
            exporter.precision(16);

            if (exporter.fail())
                throw std::runtime_error("Error on mesh cell2Ds file");

            exporter << "Cell2D_Index" << separator;
            exporter << "NumQuadPoints_Boundary" << separator;
            exporter << "NumQuadPoints_Internal" << separator;
            exporter << "Vmatrix_Cond" << separator;
            exporter << "Hmatrix_Cond" << separator;
            exporter << "Pi0k_Cond" << separator;
            exporter << "Gmatrix_Cond" << separator;
            exporter << "Pi0k_Error" << separator;
            exporter << "GBD_Error" << separator;
            exporter << "Stab_Error" << std::endl;

            for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
            {
                const auto &cell2DPerformance = vemPerformance.Cell2DsPerformance[v].VEM_Performance_Data.Analysis;

                exporter << std::scientific << v << separator;
                exporter << std::scientific << vemPerformance.Cell2DsPerformance[v].VEM_Performance_Data.NumBoundaryQuadraturePoints
                         << separator;
                exporter << std::scientific << vemPerformance.Cell2DsPerformance[v].VEM_Performance_Data.NumInternalQuadraturePoints
                         << separator;
                exporter << std::scientific << cell2DPerformance.VmatrixConditioning << separator;
                exporter << std::scientific << cell2DPerformance.HmatrixConditioning << separator;
                exporter << std::scientific << cell2DPerformance.Pi0kConditioning << separator;
                exporter << std::scientific << cell2DPerformance.GmatrixConditioning << separator;
                exporter << std::scientific << cell2DPerformance.ErrorPi0k << separator;
                exporter << std::scientific << cell2DPerformance.ErrorGBD << separator;
                exporter << std::scientific << cell2DPerformance.ErrorStabilization << std::endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
