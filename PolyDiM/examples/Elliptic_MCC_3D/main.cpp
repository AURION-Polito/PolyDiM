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
#include "I_VEM_MCC_3D_ReferenceElement.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

unsigned int Polydim::examples::Elliptic_MCC_3D::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::Elliptic_MCC_3D::Program_configuration config;

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

    Polydim::examples::Elliptic_MCC_3D::test::Patch_Test::order = config.VemOrder();
    const auto test = Polydim::examples::Elliptic_MCC_3D::program_utilities::create_test(config);

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

    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Elliptic_MCC_3D::program_utilities::create_domain_mesh(config, domain, mesh);

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
        Polydim::examples::Elliptic_MCC_3D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + std::to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

    const auto vem_pressure_reference_element = Polydim::VEM::MCC::create_VEM_MCC_3D_pressure_reference_element(config.VemType());
    const auto pressure_reference_element_data = vem_pressure_reference_element->Create(config.VemOrder());
    const auto vem_velocity_reference_element = Polydim::VEM::MCC::create_VEM_MCC_3D_velocity_reference_element(config.VemType());
    const auto velocity_reference_element_data = vem_velocity_reference_element->Create(config.VemOrder());

    const auto velocity_space = Polydim::VEM::MCC::create_VEM_MCC_3D_velocity_local_space(config.VemType());
    const auto pressure_space = Polydim::VEM::MCC::create_VEM_MCC_3D_pressure_local_space(config.VemType());

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(2);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(2);

    meshDOFsInfo[0] = dofManager.Create_Constant_DOFsInfo<3>(mesh_connectivity_data,
                                                             {{velocity_reference_element_data.NumDofs0D,
                                                               velocity_reference_element_data.NumDofs1D,
                                                               velocity_reference_element_data.NumDofs2D,
                                                               velocity_reference_element_data.NumDofs3D},
                                                              boundary_info});

    dofs_data[0] = dofManager.CreateDOFs<3>(meshDOFsInfo[0], mesh_connectivity_data);

    meshDOFsInfo[1] = dofManager.Create_Constant_DOFsInfo<3>(mesh_connectivity_data,
                                                             {{
                                                                  pressure_reference_element_data.NumDofs0D,
                                                                  pressure_reference_element_data.NumDofs1D,
                                                                  pressure_reference_element_data.NumDofs2D,
                                                                  pressure_reference_element_data.NumDofs3D,
                                                              },
                                                              boundary_info});

    dofs_data[1] = dofManager.CreateDOFs<3>(meshDOFsInfo[1], mesh_connectivity_data);

    const auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(dofs_data);

    Gedim::Output::PrintGenericMessage("VEM Space with " + std::to_string(count_dofs.num_total_dofs) + " DOFs and " +
                                           std::to_string(count_dofs.num_total_strong) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + std::to_string(static_cast<unsigned int>(config.VemType())) + "...",
                                       true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Elliptic_MCC_3D::Assembler assembler;
    auto assembler_data = assembler.Assemble(config,
                                             mesh,
                                             meshGeometricData,
                                             meshDOFsInfo,
                                             dofs_data,
                                             count_dofs,
                                             velocity_reference_element_data,
                                             pressure_reference_element_data,
                                             *velocity_space,
                                             *pressure_space,
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
                                                           *velocity_space,
                                                           *pressure_space,
                                                           assembler_data,
                                                           *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_MCC_3D::program_utilities::export_solution(config, mesh, dofs_data, assembler_data, post_process_data, exportSolutionFolder, exportVtuFolder);

    Polydim::examples::Elliptic_MCC_3D::program_utilities::export_velocity_dofs(config,
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

    if (config.ComputeVEMPerformance())
    {
        const auto vemPerformance =
            assembler.ComputeVemPerformance(config, mesh, meshGeometricData, velocity_reference_element_data, *velocity_space);
        {
            const char separator = ',';
            /// Export Cell3Ds VEM performance
            std::ofstream exporter;

            const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
            const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
            exporter.open(exportSolutionFolder + "/Cell3Ds_VEMPerformance_" + std::to_string(TEST_ID) + "_" +
                          std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv");
            exporter.precision(16);

            if (exporter.fail())
                throw std::runtime_error("Error on mesh Cell3Ds file");

            exporter << "Cell3D_Index" << separator;
            exporter << "NumQuadPoints_Boundary" << separator;
            exporter << "NumQuadPoints_Internal" << separator;
            exporter << "Vmatrix_Cond" << separator;
            exporter << "Hmatrix_Cond" << separator;
            exporter << "Pi0k_Cond" << separator;
            exporter << "Gmatrix_Cond" << separator;
            exporter << "Pi0k_Error" << separator;
            exporter << "GBD_Error" << separator;
            exporter << "Stab_Error" << std::endl;

            for (unsigned int v = 0; v < vemPerformance.Cell3DsPerformance.size(); v++)
            {
                const auto &cell3DPerformance = vemPerformance.Cell3DsPerformance[v].Analysis;

                exporter << std::scientific << v << separator;
                exporter << std::scientific << vemPerformance.Cell3DsPerformance[v].NumBoundaryQuadraturePoints << separator;
                exporter << std::scientific << vemPerformance.Cell3DsPerformance[v].NumInternalQuadraturePoints << separator;
                exporter << std::scientific << cell3DPerformance.VmatrixConditioning << separator;
                exporter << std::scientific << cell3DPerformance.HmatrixConditioning << separator;
                exporter << std::scientific << cell3DPerformance.Pi0kConditioning << separator;
                exporter << std::scientific << cell3DPerformance.GmatrixConditioning << separator;
                exporter << std::scientific << cell3DPerformance.ErrorPi0k << separator;
                exporter << std::scientific << cell3DPerformance.ErrorGBD << separator;
                exporter << std::scientific << cell3DPerformance.ErrorStabilization << std::endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
