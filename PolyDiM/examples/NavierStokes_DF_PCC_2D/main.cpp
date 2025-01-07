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

unsigned int Polydim::examples::NavierStokes_DF_PCC_2D::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::NavierStokes_DF_PCC_2D::Program_configuration config;

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

    Polydim::examples::NavierStokes_DF_PCC_2D::test::Patch_Test::order = config.VemOrder();

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

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

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
            {{velocity_reference_element_data.NumDofs0D, velocity_reference_element_data.NumDofs1D, 0, 0}, boundary_info});

        dofs_data[i] = dofManager.CreateDOFs<2>(meshDOFsInfo[i], mesh_connectivity_data);
    }

    meshDOFsInfo[2] = dofManager.Create_Constant_DOFsInfo<2>(
        mesh_connectivity_data,
        {{0, 0, velocity_reference_element_data.NumDofs2D_BigOPlus + velocity_reference_element_data.NumDofs2D_Divergence, 0},
         boundary_info});

    dofs_data[2] = dofManager.CreateDOFs<2>(meshDOFsInfo[2], mesh_connectivity_data);

    meshDOFsInfo[3] = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                             {{pressure_reference_element_data.NumDofs0D,
                                                               pressure_reference_element_data.NumDofs1D,
                                                               pressure_reference_element_data.NumDofs2D,
                                                               0},
                                                              boundary_info});

    dofs_data[3] = dofManager.CreateDOFs<2>(meshDOFsInfo[3], mesh_connectivity_data);

    auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(dofs_data);
    if (count_dofs.num_total_boundary_dofs == 0)
        count_dofs.num_total_dofs += 1; // lagrange

    Gedim::Output::PrintGenericMessage("VEM Space with " + to_string(count_dofs.num_total_dofs) + " DOFs and " +
                                           to_string(count_dofs.num_total_strong) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    Gedim::Output::PrintGenericMessage("AssembleStokesSystem VEM Type " + to_string((unsigned int)config.VemType()) + "...", true);
    Gedim::Profiler::StartTime("AssembleStokesSystem");

    const auto vem_pressure_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_pressure_local_space(config.VemType());
    const auto vem_velocity_local_space = Polydim::VEM::DF_PCC::create_VEM_DF_PCC_2D_velocity_local_space(config.VemType());

    Polydim::examples::NavierStokes_DF_PCC_2D::Assembler assembler;
    auto assembler_data = assembler.AssembleStokes(config,
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

    Gedim::Profiler::StopTime("AssembleStokesSystem");
    Gedim::Output::PrintStatusProgram("AssembleStokesSystem");

    if (count_dofs.num_total_dofs > 0)
    {
        Gedim::Output::PrintGenericMessage("Factorize...", true);
        Gedim::Profiler::StartTime("Factorize");

        Gedim::Eigen_LUSolver solver;
        solver.Initialize(assembler_data.globalMatrixA);

        cout.precision(2);
        cout << scientific << assembler_data.rightHandSide << endl;

        Gedim::Profiler::StopTime("Factorize");
        Gedim::Output::PrintStatusProgram("Factorize");

        Gedim::Output::PrintGenericMessage("Solve...", true);
        Gedim::Profiler::StartTime("Solve");

        solver.Solve(assembler_data.rightHandSide, assembler_data.solution);

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");
    }

    auto post_process_data_2 = assembler.PostProcessSolution(config,
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

    Polydim::examples::NavierStokes_DF_PCC_2D::program_utilities::export_solution(config,
                                                                                  mesh,
                                                                                  dofs_data,
                                                                                  count_dofs,
                                                                                  assembler_data,
                                                                                  post_process_data_2,
                                                                                  0,
                                                                                  exportSolutionFolder,
                                                                                  exportVtuFolder);

    Gedim::Eigen_Array<> residual;
    residual.SetSize(count_dofs.num_total_dofs);
    residual.SumMultiplication(assembler_data.globalMatrixA, assembler_data.previousIteration);
    residual -= assembler_data.rightHandSide;

    const double initial_residual_norm = residual.Norm();
    const double initial_solution_norm = assembler_data.previousIteration.Norm();

    const unsigned int NLMAXNumIterations = config.NLMaxNumberIterations();
    unsigned int num_nl_iterations;
    for (unsigned int l = 0; l < NLMAXNumIterations; l++)
    {
        Gedim::Output::PrintGenericMessage("AssembleNavierStokesTerm VEM Type " + to_string((unsigned int)config.VemType()) + "...", true);
        Gedim::Profiler::StartTime("AssembleNavierStokesTerm");

        assembler.AssembleNavierStokes(config,
                                       mesh,
                                       meshGeometricData,
                                       meshDOFsInfo,
                                       dofs_data,
                                       count_dofs,
                                       velocity_reference_element_data,
                                       pressure_reference_element_data,
                                       *vem_velocity_local_space,
                                       *vem_pressure_local_space,
                                       assembler_data);

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

            solver.Solve(assembler_data.rightHandSideC, assembler_data.solution);

            Gedim::Profiler::StopTime("Solve");
            Gedim::Output::PrintStatusProgram("Solve");
        }

        residual.Zeros();
        residual.SumMultiplication(assembler_data.globalMatrixC, assembler_data.solution);
        residual -= assembler_data.rightHandSideC;
        const double residual_norm = residual.Norm();

        Gedim::Eigen_Array<> delta;
        delta.Copy(assembler_data.solution);
        delta -= assembler_data.previousIteration;

        const double delta_norm = delta.Norm();

        assembler_data.previousIteration.Copy(assembler_data.solution);

        cout << "NL It: " << l + 1 << ", Residual Norm: " << residual_norm << ", Delta norm: " << delta_norm << endl;

        if ((residual_norm <= config.NLRelResidualTolerance() * initial_residual_norm + config.NLAbsResidualTolerance()) &&
            (delta_norm <= config.NLRelChangeInSolutionTolerance() * initial_solution_norm + config.NLAbsChangeInSolutionTolerance()))
        {
            num_nl_iterations = l + 1;
            break;
        }
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
            ofstream exporter;

            const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
            const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
            exporter.open(exportSolutionFolder + "/Cell2Ds_VEMPerformance_" + to_string(TEST_ID) + "_" +
                          to_string(VEM_ID) + +"_" + to_string(config.VemOrder()) + ".csv");
            exporter.precision(16);

            if (exporter.fail())
                throw runtime_error("Error on mesh cell2Ds file");

            exporter << "Cell2D_Index" << separator;
            exporter << "NumQuadPoints_Boundary" << separator;
            exporter << "NumQuadPoints_Internal" << separator;
            exporter << "max_PiNabla_Cond" << separator;
            exporter << "max_Pi0k_Cond" << separator;
            exporter << "max_PiNabla_Error" << separator;
            exporter << "max_Pi0k_Error" << separator;
            exporter << "max_GBD_Error" << separator;
            exporter << "max_HCD_Error" << separator;
            exporter << "Stab_Error" << endl;

            for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
            {
                const auto &cell2DPerformance = vemPerformance.Cell2DsPerformance[v];

                exporter << scientific << v << separator;
                exporter << scientific << cell2DPerformance.NumBoundaryQuadraturePoints << separator;
                exporter << scientific << cell2DPerformance.NumInternalQuadraturePoints << separator;
                exporter << scientific << cell2DPerformance.maxPiNablaConditioning << separator;
                exporter << scientific << cell2DPerformance.maxPi0kConditioning << separator;
                exporter << scientific << cell2DPerformance.maxErrorPiNabla << separator;
                exporter << scientific << cell2DPerformance.maxErrorPi0k << separator;
                exporter << scientific << cell2DPerformance.maxErrorGBD << separator;
                exporter << scientific << cell2DPerformance.maxErrorHCD << separator;
                exporter << scientific << cell2DPerformance.ErrorStabilization << endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
