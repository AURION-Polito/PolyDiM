#include "Assembler_Utilities.hpp"
#include "Eigen_CholeskySolver.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "VTKUtilities.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

int Polydim::examples::Elastic_PCC_2D::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::Elastic_PCC_2D::Program_configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini", false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    /// Create folders
    const string exportFolder = config.ExportFolder();
    Gedim::Output::CreateFolder(exportFolder);

    const string exportCsvFolder = exportFolder + "/Mesh";
    Gedim::Output::CreateFolder(exportCsvFolder);
    const string exportVtuFolder = exportFolder + "/Paraview";
    Gedim::Output::CreateFolder(exportVtuFolder);
    const string exportSolutionFolder = exportFolder + "/Solution";
    Gedim::Output::CreateFolder(exportSolutionFolder);

    const string logFolder = exportFolder + "/Log";

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

    Polydim::examples::Elastic_PCC_2D::test::Patch_Test::order = config.MethodOrder();

    const auto test = Polydim::examples::Elastic_PCC_2D::program_utilities::create_test(config);

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

    Polydim::examples::Elastic_PCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

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
        Polydim::examples::Elastic_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    const auto reference_element_data =
        Polydim::examples::Elastic_PCC_2D::local_space::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

    const auto reference_element_num_dofs =
        Polydim::examples::Elastic_PCC_2D::local_space::ReferenceElementNumDOFs(reference_element_data);

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    const auto meshDOFsInfo =
        dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data, {reference_element_num_dofs, boundary_info});

    const auto dofs_data = dofManager.CreateDOFs<2>(meshDOFsInfo, mesh_connectivity_data);

    const auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs({std::cref(dofs_data), std::cref(dofs_data)});

    Gedim::Output::PrintGenericMessage("Discrete Space with " + to_string(count_dofs.num_total_dofs) + " DOFs and " +
                                           to_string(count_dofs.num_total_strong) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem Discrete Type " +
                                           to_string(static_cast<unsigned int>(config.MethodType())) + "...",
                                       true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Elastic_PCC_2D::Assembler assembler;
    auto assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, count_dofs, reference_element_data, *test);

    Gedim::Profiler::StopTime("AssembleSystem");
    Gedim::Output::PrintStatusProgram("AssembleSystem");

    if (dofs_data.NumberDOFs > 0)
    {
        Gedim::Output::PrintGenericMessage("Factorize...", true);
        Gedim::Profiler::StartTime("Factorize");

        Gedim::Eigen_CholeskySolver solver;
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

    Polydim::examples::Elastic_PCC_2D::program_utilities::export_solution(config,
                                                                          mesh,
                                                                          dofs_data,
                                                                          count_dofs,
                                                                          assembler_data,
                                                                          post_process_data,
                                                                          exportSolutionFolder,
                                                                          exportVtuFolder);

    Polydim::examples::Elastic_PCC_2D::program_utilities::export_dofs(config,
                                                                      mesh,
                                                                      meshGeometricData,
                                                                      meshDOFsInfo,
                                                                      dofs_data,
                                                                      count_dofs,
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

        Polydim::examples::Elastic_PCC_2D::program_utilities::export_performance(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeMethodPerformance");
    Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");

    return 0;
}
