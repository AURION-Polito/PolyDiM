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
#include "MeshMatricesDAO.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"

unsigned int Polydim::examples::Stokes_DF_PCC_3D::test::Patch_Test::order;
unsigned int Polydim::examples::Stokes_DF_PCC_3D::test::Stokes_Benchmark_1::order;
unsigned int Polydim::examples::Stokes_DF_PCC_3D::test::Stokes_Benchmark_2::order;

int main(int argc, char **argv)
{
    Polydim::examples::Stokes_DF_PCC_3D::Program_configuration config;

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

    Polydim::examples::Stokes_DF_PCC_3D::test::Patch_Test::order = config.VemOrder();
    Polydim::examples::Stokes_DF_PCC_3D::test::Stokes_Benchmark_1::order = config.VemOrder();
    Polydim::examples::Stokes_DF_PCC_3D::test::Stokes_Benchmark_2::order = config.VemOrder();

    const auto test = Polydim::examples::Stokes_DF_PCC_3D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolyhedron(domain.vertices, domain.edges, domain.faces);
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

    Polydim::examples::Stokes_DF_PCC_3D::program_utilities::create_domain_mesh(config, domain, mesh);

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    // Export the domain mesh
    {
        Gedim::MeshUtilities meshUtilities;
        meshUtilities.ExportMeshToVTU(mesh, exportVtuFolder, "Domain_Mesh");
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto mesh_geometric_data =
        Polydim::examples::Stokes_DF_PCC_3D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + std::to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

    const auto vem_velocity_reference_element_2D =
        Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_velocity_reference_element_2D(config.VemType());
    const auto velocity_reference_element_data_2D = vem_velocity_reference_element_2D->Create(config.VemOrder());
    const auto vem_pressure_reference_element_3D =
        Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_pressure_reference_element_3D(config.VemType());
    const auto pressure_reference_element_data_3D = vem_pressure_reference_element_3D->Create(config.VemOrder());
    const auto vem_velocity_reference_element_3D =
        Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_velocity_reference_element_3D(config.VemType());
    const auto velocity_reference_element_data_3D = vem_velocity_reference_element_3D->Create(config.VemOrder());

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> mesh_dofs_info(8);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(8);

    for (unsigned int i = 0; i < 3; i++)
    {
        mesh_dofs_info[i] = dofManager.Create_Constant_DOFsInfo<3>(
            mesh_connectivity_data,
            {{velocity_reference_element_data_3D.NumDofs0D, velocity_reference_element_data_3D.NumDofs1D, 0, 0}, boundary_info});

        dofs_data[i] = dofManager.CreateDOFs<3>(mesh_dofs_info[i], mesh_connectivity_data);
    }

    for (unsigned int i = 3; i < 6; i++)
    {
        mesh_dofs_info[i] =
            dofManager.Create_Constant_DOFsInfo<3>(mesh_connectivity_data,
                                                   {{0, 0, velocity_reference_element_data_3D.NumDofs2D, 0}, boundary_info});

        dofs_data[i] = dofManager.CreateDOFs<3>(mesh_dofs_info[i], mesh_connectivity_data);
    }

    mesh_dofs_info[6] = dofManager.Create_Constant_DOFsInfo<3>(
        mesh_connectivity_data,
        {{0, 0, 0, velocity_reference_element_data_3D.NumDofs3D_BigOPlus + velocity_reference_element_data_3D.NumDofs3D_Divergence},
         boundary_info});

    dofs_data[6] = dofManager.CreateDOFs<3>(mesh_dofs_info[6], mesh_connectivity_data);

    mesh_dofs_info[7] = dofManager.Create_Constant_DOFsInfo<3>(mesh_connectivity_data,
                                                               {{pressure_reference_element_data_3D.NumDofs0D,
                                                                 pressure_reference_element_data_3D.NumDofs1D,
                                                                 pressure_reference_element_data_3D.NumDofs2D,
                                                                 pressure_reference_element_data_3D.NumDofs3D},
                                                                boundary_info});

    dofs_data[7] = dofManager.CreateDOFs<3>(mesh_dofs_info[7], mesh_connectivity_data);

    auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(dofs_data);
    if (count_dofs.num_total_boundary_dofs == 0)
        count_dofs.num_total_dofs += 1; // lagrange

    Gedim::Output::PrintGenericMessage("VEM Space with " + std::to_string(count_dofs.num_total_dofs) + " DOFs and " +
                                           std::to_string(count_dofs.num_total_strong) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + std::to_string((unsigned int)config.VemType()) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    const auto vem_pressure_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_pressure_local_space_3D(config.VemType());
    const auto vem_velocity_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_velocity_local_space_3D(config.VemType());

    Polydim::examples::Stokes_DF_PCC_3D::Assembler assembler;
    auto assembler_data = assembler.Assemble(config,
                                             mesh,
                                             mesh_geometric_data,
                                             mesh_dofs_info,
                                             dofs_data,
                                             count_dofs,
                                             velocity_reference_element_data_2D,
                                             velocity_reference_element_data_3D,
                                             pressure_reference_element_data_3D,
                                             *vem_velocity_local_space,
                                             *vem_pressure_local_space,
                                             *test);

    Gedim::Profiler::StopTime("AssembleSystem");
    Gedim::Output::PrintStatusProgram("AssembleSystem");

    Gedim::Output::PrintGenericMessage("Solve...", true);
    Gedim::Profiler::StartTime("Solve");

    if (count_dofs.num_total_dofs > 0)
        Polydim::examples::Stokes_DF_PCC_3D::program_utilities::solve_stokes(config, assembler_data);

    Gedim::Profiler::StopTime("Solve");
    Gedim::Output::PrintStatusProgram("Solve");

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    auto post_process_data = assembler.PostProcessSolution(config,
                                                           mesh,
                                                           mesh_geometric_data,
                                                           dofs_data,
                                                           count_dofs,
                                                           velocity_reference_element_data_2D,
                                                           velocity_reference_element_data_3D,
                                                           pressure_reference_element_data_3D,
                                                           *vem_velocity_local_space,
                                                           *vem_pressure_local_space,
                                                           assembler_data,
                                                           *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Stokes_DF_PCC_3D::program_utilities::export_solution(config,
                                                                            mesh,
                                                                            dofs_data,
                                                                            count_dofs,
                                                                            assembler_data,
                                                                            post_process_data,
                                                                            exportSolutionFolder,
                                                                            exportVtuFolder);

    Polydim::examples::Stokes_DF_PCC_3D::program_utilities::export_velocity_dofs(config,
                                                                                 mesh,
                                                                                 mesh_geometric_data,
                                                                                 mesh_dofs_info,
                                                                                 velocity_reference_element_data_3D,
                                                                                 dofs_data,
                                                                                 count_dofs,
                                                                                 assembler_data,
                                                                                 post_process_data,
                                                                                 exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
    Gedim::Profiler::StartTime("ComputeVEMPerformance");

    if (config.ComputeVEMPerformance())
    {
        const auto performance = assembler.ComputeVemPerformance(config,
                                                                 mesh,
                                                                 mesh_geometric_data,
                                                                 velocity_reference_element_data_2D,
                                                                 velocity_reference_element_data_3D,
                                                                 *vem_velocity_local_space);
        Polydim::examples::Stokes_DF_PCC_3D::program_utilities::export_performance(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    if (config.ComputeDiscrepancyError())
    {
        Gedim::Output::PrintGenericMessage("Create Full VEM Space of order " + std::to_string(config.VemOrder()) + " and DOFs...", true);
        Gedim::Profiler::StartTime("CreateFULLVEMSpace");

        const auto vem_full_pressure_reference_element_3D =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_full_pressure_reference_element_3D(config.VemType());
        const auto full_pressure_reference_element_data_3D = vem_full_pressure_reference_element_3D->Create(config.VemOrder());
        const auto vem_full_velocity_reference_element_3D =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_full_velocity_reference_element_3D(config.VemType());
        const auto full_velocity_reference_element_data_3D = vem_full_velocity_reference_element_3D->Create(config.VemOrder());

        std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> full_mesh_dofs_info(8);
        std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> full_dofs_data(8);

        for (unsigned int i = 0; i < 3; i++)
        {
            full_mesh_dofs_info[i] = dofManager.Create_Constant_DOFsInfo<3>(
                mesh_connectivity_data,
                {{full_velocity_reference_element_data_3D.NumDofs0D, full_velocity_reference_element_data_3D.NumDofs1D, 0, 0},
                 boundary_info});

            full_dofs_data[i] = dofManager.CreateDOFs<3>(full_mesh_dofs_info[i], mesh_connectivity_data);
        }

        for (unsigned int i = 3; i < 6; i++)
        {
            full_mesh_dofs_info[i] = dofManager.Create_Constant_DOFsInfo<3>(
                mesh_connectivity_data,
                {{0, 0, full_velocity_reference_element_data_3D.NumDofs2D, 0}, boundary_info});

            full_dofs_data[i] = dofManager.CreateDOFs<3>(full_mesh_dofs_info[i], mesh_connectivity_data);
        }

        full_mesh_dofs_info[6] = dofManager.Create_Constant_DOFsInfo<3>(
            mesh_connectivity_data,
            {{0,
              0,
              0,
              full_velocity_reference_element_data_3D.NumDofs3D_BigOPlus + full_velocity_reference_element_data_3D.NumDofs3D_Divergence},
             boundary_info});

        full_dofs_data[6] = dofManager.CreateDOFs<3>(full_mesh_dofs_info[6], mesh_connectivity_data);

        full_mesh_dofs_info[7] = dofManager.Create_Constant_DOFsInfo<3>(mesh_connectivity_data,
                                                                        {{full_pressure_reference_element_data_3D.NumDofs0D,
                                                                          full_pressure_reference_element_data_3D.NumDofs1D,
                                                                          full_pressure_reference_element_data_3D.NumDofs2D,
                                                                          full_pressure_reference_element_data_3D.NumDofs3D},
                                                                         boundary_info});

        full_dofs_data[7] = dofManager.CreateDOFs<3>(full_mesh_dofs_info[7], mesh_connectivity_data);

        auto full_count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(full_dofs_data);
        if (full_count_dofs.num_total_boundary_dofs == 0)
            full_count_dofs.num_total_dofs += 1; // lagrange

        Gedim::Output::PrintGenericMessage("VEM Space with " + std::to_string(full_count_dofs.num_total_dofs) +
                                               " DOFs and " + std::to_string(full_count_dofs.num_total_strong) + " STRONGs",
                                           true);

        Gedim::Profiler::StopTime("CreateFULLVEMSpace");
        Gedim::Output::PrintStatusProgram("CreateFULLVEMSpace");

        Gedim::Output::PrintGenericMessage("AssembleSystem FULL VEM Type " + std::to_string((unsigned int)config.VemType()) + "...", true);
        Gedim::Profiler::StartTime("AssembleSystem");

        const auto vem_full_pressure_local_space =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_full_pressure_local_space_3D(config.VemType());
        const auto vem_full_velocity_local_space =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_3D_full_velocity_local_space_3D(config.VemType());

        auto full_assembler_data = assembler.Assemble(config,
                                                      mesh,
                                                      mesh_geometric_data,
                                                      full_mesh_dofs_info,
                                                      full_dofs_data,
                                                      full_count_dofs,
                                                      velocity_reference_element_data_2D,
                                                      full_velocity_reference_element_data_3D,
                                                      full_pressure_reference_element_data_3D,
                                                      *vem_full_velocity_local_space,
                                                      *vem_full_pressure_local_space,
                                                      *test);

        Gedim::Profiler::StopTime("AssembleSystem");
        Gedim::Output::PrintStatusProgram("AssembleSystem");

        Gedim::Output::PrintGenericMessage("Solve...", true);
        Gedim::Profiler::StartTime("Solve");

        if (full_count_dofs.num_total_dofs > 0)
            Polydim::examples::Stokes_DF_PCC_3D::program_utilities::solve_stokes(config, full_assembler_data);

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");

        Gedim::Output::PrintGenericMessage("ComputeDiscrepancyErrors...", true);
        Gedim::Profiler::StartTime("ComputeDiscrepancyErrors");

        const auto discrepancy_errors_data = assembler.ComputeDiscrepancyErrors(config,
                                                                                mesh,
                                                                                mesh_geometric_data,
                                                                                full_dofs_data,
                                                                                full_count_dofs,
                                                                                dofs_data,
                                                                                count_dofs,
                                                                                velocity_reference_element_data_2D,
                                                                                full_velocity_reference_element_data_3D,
                                                                                full_pressure_reference_element_data_3D,
                                                                                velocity_reference_element_data_3D,
                                                                                pressure_reference_element_data_3D,
                                                                                *vem_full_velocity_local_space,
                                                                                *vem_full_pressure_local_space,
                                                                                *vem_velocity_local_space,
                                                                                *vem_pressure_local_space,
                                                                                full_assembler_data,
                                                                                assembler_data);

        Polydim::examples::Stokes_DF_PCC_3D::program_utilities::export_discrepancy_errors(config,
                                                                                          mesh,
                                                                                          discrepancy_errors_data,
                                                                                          exportSolutionFolder,
                                                                                          exportVtuFolder);

        Gedim::Profiler::StopTime("ComputeDiscrepancyErrors");
        Gedim::Output::PrintStatusProgram("ComputeDiscrepancyErrors");
    }

    return 0;
}
