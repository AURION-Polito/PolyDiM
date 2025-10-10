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

unsigned int Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::space_order;
unsigned int Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::time_order;

int main(int argc, char **argv)
{
    Polydim::examples::Parabolic_PCC_2D::Program_configuration config;

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

    Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::space_order = config.MethodOrder();
    Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::time_order = config.Theta() == 0.5 ? 2 : 1;

    const auto test = Polydim::examples::Parabolic_PCC_2D::program_utilities::create_test(config);

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

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Parabolic_PCC_2D::program_utilities::create_domain_mesh(config, domain.spatial_domain, mesh);
    const auto time_steps = Polydim::examples::Parabolic_PCC_2D::program_utilities::create_time_steps(config,
                                                                                                      domain.time_domain);
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
        Polydim::examples::Parabolic_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

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

    Polydim::examples::Parabolic_PCC_2D::Assembler assembler;

    const auto static_assembler_data = assembler.StaticAssemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, reference_element_data, *test);
    const auto initial_condition = assembler.ComputeInitalCondition(config, mesh, meshGeometricData, dofs_data, reference_element_data, *test);
    auto initial_assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, reference_element_data,
        *test, static_assembler_data, time_steps[0]);

    initial_assembler_data.solution = initial_condition.initial_condition;
    initial_assembler_data.solutionDirichlet = initial_condition.initial_condition_dirichlet;

    auto u_k = initial_condition.initial_condition;
    auto u_D_k = initial_condition.initial_condition_dirichlet;
    auto f_k = initial_assembler_data.rightHandSide;
    const auto theta = config.Theta();

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    auto initial_post_process_data =
        assembler.PostProcessSolution(config, mesh, meshGeometricData, dofs_data, reference_element_data,
        initial_assembler_data, *test, static_assembler_data.globalMatrixA, f_k, time_steps.at(0));

    Polydim::examples::Parabolic_PCC_2D::program_utilities::export_solution(config,
                                                                            mesh,
                                                                            dofs_data,
                                                                            static_assembler_data.globalMatrixA,
                                                                            initial_post_process_data,
                                                                            0,
                                                                            time_steps.at(0),
                                                                            exportSolutionFolder,
                                                                            exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    for (unsigned int t = 1; t < time_steps.size(); t++)
    {
        const double time_value = time_steps.at(t);
        const double dt = time_steps.at(t) - time_steps.at(t - 1);

        auto Kp1 = static_assembler_data.globalMatrixA;
        Kp1 *= dt * theta;
        Kp1 += static_assembler_data.globalMatrixM;

        auto Kp1_D = static_assembler_data.dirichletMatrixA;
        Kp1_D *= dt * theta;
        Kp1_D += static_assembler_data.dirichletMatrixM;

        auto K = static_assembler_data.globalMatrixA;
        K *= -dt * (1.0 - theta);
        K += static_assembler_data.globalMatrixM;

        auto K_D = static_assembler_data.dirichletMatrixA;
        K_D *= -dt * (1.0 - theta);
        K_D += static_assembler_data.dirichletMatrixM;

        auto assembler_data_kp1 =
            assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, reference_element_data,
            *test, static_assembler_data, time_value);

        auto& u_kp1 = assembler_data_kp1.solution;
        auto& u_D_kp1 = assembler_data_kp1.solutionDirichlet;
        auto& f_kp1 = assembler_data_kp1.rightHandSide;

        auto rhs_f_kp1 = assembler_data_kp1.rightHandSide;
        rhs_f_kp1 *= (dt * theta);
        auto rhs_f_k = f_k;
        rhs_f_k *= dt * (1.0 - theta);
        auto rhs = rhs_f_kp1 +
                   rhs_f_k;
        rhs.SumMultiplication(K,
                              u_k);
        rhs.SumMultiplication(K_D,
                              u_D_k);
        rhs.SubtractionMultiplication(K_D,
                                      u_D_kp1);

        Gedim::Profiler::StopTime("AssembleSystem");
        Gedim::Output::PrintStatusProgram("AssembleSystem");

        if (dofs_data.NumberDOFs > 0)
        {
            Gedim::Output::PrintGenericMessage("Factorize...", true);
            Gedim::Profiler::StartTime("Factorize");

            Gedim::Eigen_LUSolver solver;
            solver.Initialize(Kp1);

            Gedim::Profiler::StopTime("Factorize");
            Gedim::Output::PrintStatusProgram("Factorize");

            Gedim::Output::PrintGenericMessage("Solve...", true);
            Gedim::Profiler::StartTime("Solve");

            solver.Solve(rhs, u_kp1);

            Gedim::Profiler::StopTime("Solve");
            Gedim::Output::PrintStatusProgram("Solve");
        }

        Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
        Gedim::Profiler::StartTime("ComputeErrors");

        auto post_process_data =
            assembler.PostProcessSolution(config, mesh, meshGeometricData, dofs_data, reference_element_data,
            assembler_data_kp1, *test, Kp1, rhs, time_value);

        Gedim::Profiler::StopTime("ComputeErrors");
        Gedim::Output::PrintStatusProgram("ComputeErrors");

        Gedim::Output::PrintGenericMessage("ExportSolution...", true);
        Gedim::Profiler::StartTime("ExportSolution");

        Polydim::examples::Parabolic_PCC_2D::program_utilities::export_solution(config,
                                                                                mesh,
                                                                                dofs_data,
                                                                                Kp1,
                                                                                post_process_data,
                                                                                t,
                                                                                time_value,
                                                                                exportSolutionFolder,
                                                                                exportVtuFolder);

        Polydim::examples::Parabolic_PCC_2D::program_utilities::export_dofs(config,
                                                                            mesh,
                                                                            meshGeometricData,
                                                                            meshDOFsInfo,
                                                                            dofs_data,
                                                                            reference_element_data,
                                                                            assembler_data_kp1,
                                                                            post_process_data,
                                                                            exportVtuFolder);

        Gedim::Profiler::StopTime("ExportSolution");
        Gedim::Output::PrintStatusProgram("ExportSolution");

        u_k = u_kp1;
        u_D_k = u_D_kp1;
        f_k = f_kp1;
    }

    return 0;
}
