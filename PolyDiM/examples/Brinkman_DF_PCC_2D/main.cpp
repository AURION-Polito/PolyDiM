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
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"

unsigned int Polydim::examples::Brinkman_DF_PCC_2D::test::Patch_Test::order;

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

    Polydim::examples::Brinkman_DF_PCC_2D::test::Patch_Test::order = config.VemOrder();

    const auto test = Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::create_test(config);

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

    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

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
        Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + std::to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data{mesh};

    const auto vem_pressure_reference_element =
        Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_pressure_reference_element(config.VemType());
    const auto pressure_reference_element_data = vem_pressure_reference_element->Create(config.VemOrder());
    const auto vem_velocity_reference_element =
        Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_velocity_reference_element(config.VemType());
    const auto velocity_reference_element_data = vem_velocity_reference_element->Create(config.VemOrder());

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(4);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(4);

    for (unsigned int i = 0; i < 2; i++)
    {
        meshDOFsInfo[i] = dofManager.Create_Constant_DOFsInfo<2>(
            mesh_connectivity_data,
            {{velocity_reference_element_data.NumDofs0D, velocity_reference_element_data.NumDofs1D, 0, 0}, boundary_info[i]});

        dofs_data[i] = dofManager.CreateDOFs<2>(meshDOFsInfo[i], mesh_connectivity_data);
    }

    meshDOFsInfo[2] = dofManager.Create_Constant_DOFsInfo<2>(
        mesh_connectivity_data,
        {{0, 0, velocity_reference_element_data.NumDofs2D_BigOPlus + velocity_reference_element_data.NumDofs2D_Divergence, 0},
         boundary_info[2]});

    dofs_data[2] = dofManager.CreateDOFs<2>(meshDOFsInfo[2], mesh_connectivity_data);

    meshDOFsInfo[3] = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                             {{pressure_reference_element_data.NumDofs0D,
                                                               pressure_reference_element_data.NumDofs1D,
                                                               pressure_reference_element_data.NumDofs2D,
                                                               0},
                                                              boundary_info[3]});

    dofs_data[3] = dofManager.CreateDOFs<2>(meshDOFsInfo[3], mesh_connectivity_data);

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

    const auto vem_pressure_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_pressure_local_space(config.VemType());
    const auto vem_velocity_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_velocity_local_space(config.VemType());

    Polydim::examples::Brinkman_DF_PCC_2D::Assembler assembler;
    auto assembler_data = assembler.Assemble(config,
                                             mesh,
                                             meshGeometricData,
                                             meshDOFsInfo,
                                             dofs_data,
                                             count_dofs,
                                             velocity_reference_element_data,
                                             pressure_reference_element_data,
                                             *vem_velocity_local_space,
                                             *vem_pressure_local_space,
                                             *test);

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

    auto post_process_data = assembler.PostProcessSolution(config,
                                                           mesh,
                                                           meshGeometricData,
                                                           dofs_data,
                                                           count_dofs,
                                                           velocity_reference_element_data,
                                                           pressure_reference_element_data,
                                                           *vem_velocity_local_space,
                                                           *vem_pressure_local_space,
                                                           assembler_data,
                                                           *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::export_solution(config,
                                                                              mesh,
                                                                              dofs_data,
                                                                              count_dofs,
                                                                              assembler_data,
                                                                              post_process_data,
                                                                              exportSolutionFolder,
                                                                              exportVtuFolder);

    Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::export_velocity_dofs(config,
                                                                                   mesh,
                                                                                   meshGeometricData,
                                                                                   meshDOFsInfo,
                                                                                   velocity_reference_element_data,
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
        const auto vemPerformance =
            assembler.ComputeVemPerformance(config, mesh, meshGeometricData, velocity_reference_element_data, *vem_velocity_local_space);
        {
            const char separator = ',';
            /// Export Cell2Ds VEM performance
            std::ofstream exporter;

            const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
            const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
            exporter.open(exportSolutionFolder + "/Cell2Ds_VEMPerformance_" + std::to_string(TEST_ID) + "_" +
                          std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv");
            exporter.precision(16);

            if (exporter.fail())
                throw std::runtime_error("Error on mesh cell2Ds file");

            exporter << "Cell2D_Index" << separator;
            exporter << "NumQuadPoints_Boundary" << separator;
            exporter << "NumQuadPoints_Internal" << separator;
            exporter << "max_PiNabla_Cond" << separator;
            exporter << "max_Pi0k_Cond" << separator;
            exporter << "max_PiNabla_Error" << separator;
            exporter << "max_Pi0k_Error" << separator;
            exporter << "max_GBD_Error" << separator;
            exporter << "max_HCD_Error" << separator;
            exporter << "Stab_Error" << std::endl;

            for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
            {
                const auto &cell2DPerformance = vemPerformance.Cell2DsPerformance[v];

                exporter << std::scientific << v << separator;
                exporter << std::scientific << cell2DPerformance.NumBoundaryQuadraturePoints << separator;
                exporter << std::scientific << cell2DPerformance.NumInternalQuadraturePoints << separator;
                exporter << std::scientific << cell2DPerformance.maxPiNablaConditioning << separator;
                exporter << std::scientific << cell2DPerformance.maxPi0kConditioning << separator;
                exporter << std::scientific << cell2DPerformance.maxErrorPiNabla << separator;
                exporter << std::scientific << cell2DPerformance.maxErrorPi0k << separator;
                exporter << std::scientific << cell2DPerformance.maxErrorGBD << separator;
                exporter << std::scientific << cell2DPerformance.maxErrorHCD << separator;
                exporter << std::scientific << cell2DPerformance.ErrorStabilization << std::endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    if (config.ComputeDiscrepancyError())
    {
        Gedim::Output::PrintGenericMessage("Create Full VEM Space of order " + std::to_string(config.VemOrder()) + " and DOFs...", true);
        Gedim::Profiler::StartTime("CreateFULLVEMSpace");

        const auto vem_full_pressure_reference_element =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_full_pressure_reference_element(config.VemType());
        const auto full_pressure_reference_element_data = vem_full_pressure_reference_element->Create(config.VemOrder());
        const auto vem_full_velocity_reference_element =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_full_velocity_reference_element(config.VemType());
        const auto full_velocity_reference_element_data = vem_full_velocity_reference_element->Create(config.VemOrder());

        std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> full_meshDOFsInfo(4);
        std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> full_dofs_data(4);

        for (unsigned int i = 0; i < 2; i++)
        {
            full_meshDOFsInfo[i] = dofManager.Create_Constant_DOFsInfo<2>(
                mesh_connectivity_data,
                {{full_velocity_reference_element_data.NumDofs0D, full_velocity_reference_element_data.NumDofs1D, 0, 0},
                 boundary_info[i]});

            full_dofs_data[i] = dofManager.CreateDOFs<2>(full_meshDOFsInfo[i], mesh_connectivity_data);
        }

        full_meshDOFsInfo[2] = dofManager.Create_Constant_DOFsInfo<2>(
            mesh_connectivity_data,
            {{0, 0, full_velocity_reference_element_data.NumDofs2D_BigOPlus + full_velocity_reference_element_data.NumDofs2D_Divergence, 0},
             boundary_info[2]});

        full_dofs_data[2] = dofManager.CreateDOFs<2>(full_meshDOFsInfo[2], mesh_connectivity_data);

        full_meshDOFsInfo[3] = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                                      {{full_pressure_reference_element_data.NumDofs0D,
                                                                        full_pressure_reference_element_data.NumDofs1D,
                                                                        full_pressure_reference_element_data.NumDofs2D,
                                                                        0},
                                                                       boundary_info[3]});

        full_dofs_data[3] = dofManager.CreateDOFs<2>(full_meshDOFsInfo[3], mesh_connectivity_data);

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
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_full_pressure_local_space(config.VemType());
        const auto vem_full_velocity_local_space =
            Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_full_velocity_local_space(config.VemType());

        auto full_assembler_data = assembler.Assemble(config,
                                                      mesh,
                                                      meshGeometricData,
                                                      full_meshDOFsInfo,
                                                      full_dofs_data,
                                                      full_count_dofs,
                                                      full_velocity_reference_element_data,
                                                      full_pressure_reference_element_data,
                                                      *vem_full_velocity_local_space,
                                                      *vem_full_pressure_local_space,
                                                      *test);

        Gedim::Profiler::StopTime("AssembleSystem");
        Gedim::Output::PrintStatusProgram("AssembleSystem");

        if (full_count_dofs.num_total_dofs > 0)
        {
            Gedim::Output::PrintGenericMessage("Factorize...", true);
            Gedim::Profiler::StartTime("Factorize");

            Gedim::Eigen_LUSolver solver;
            solver.Initialize(full_assembler_data.globalMatrixA);

            Gedim::Profiler::StopTime("Factorize");
            Gedim::Output::PrintStatusProgram("Factorize");

            Gedim::Output::PrintGenericMessage("Solve...", true);
            Gedim::Profiler::StartTime("Solve");

            solver.Solve(full_assembler_data.rightHandSide, full_assembler_data.solution);

            Gedim::Profiler::StopTime("Solve");
            Gedim::Output::PrintStatusProgram("Solve");
        }

        Gedim::Output::PrintGenericMessage("ComputeDiscrepancyErrors...", true);
        Gedim::Profiler::StartTime("ComputeDiscrepancyErrors");

        const auto discrepancy_errors_data = assembler.ComputeDiscrepancyErrors(config,
                                                                                mesh,
                                                                                meshGeometricData,
                                                                                full_dofs_data,
                                                                                full_count_dofs,
                                                                                dofs_data,
                                                                                count_dofs,
                                                                                full_velocity_reference_element_data,
                                                                                full_pressure_reference_element_data,
                                                                                velocity_reference_element_data,
                                                                                pressure_reference_element_data,
                                                                                *vem_full_velocity_local_space,
                                                                                *vem_full_pressure_local_space,
                                                                                *vem_velocity_local_space,
                                                                                *vem_pressure_local_space,
                                                                                full_assembler_data,
                                                                                assembler_data);

        Polydim::examples::Brinkman_DF_PCC_2D::program_utilities::export_discrepancy_errors(config,
                                                                                            mesh,
                                                                                            discrepancy_errors_data,
                                                                                            exportSolutionFolder,
                                                                                            exportVtuFolder);

        Gedim::Profiler::StopTime("ComputeDiscrepancyErrors");
        Gedim::Output::PrintStatusProgram("ComputeDiscrepancyErrors");
    }

    return 0;
}
