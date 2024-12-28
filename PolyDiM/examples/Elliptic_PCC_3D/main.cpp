#include "Eigen_CholeskySolver.hpp"

#include "MeshMatricesDAO_mesh_connectivity_data.hpp"

#include "program_utilities.hpp"
#include "test_definition.hpp"

unsigned int Polydim::examples::Elliptic_PCC_3D::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::Elliptic_PCC_3D::Program_configuration config;

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

    Polydim::examples::Elliptic_PCC_3D::test::Patch_Test::order = config.VemOrder();

    const auto test = Polydim::examples::Elliptic_PCC_3D::program_utilities::create_test(config);

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

    Polydim::examples::Elliptic_PCC_3D::program_utilities::create_domain_mesh(config, domain, mesh);

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
        Polydim::examples::Elliptic_PCC_3D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    const auto vem_reference_element_2D = Polydim::VEM::PCC::create_VEM_PCC_3D_reference_element_2D(config.VemType());
    const auto reference_element_data_2D = vem_reference_element_2D->Create(config.VemOrder());
    const auto vem_reference_element_3D = Polydim::VEM::PCC::create_VEM_PCC_3D_reference_element_3D(config.VemType());
    const auto reference_element_data_3D = vem_reference_element_3D->Create(config.VemOrder());

    const auto local_space = Polydim::VEM::PCC::create_VEM_PCC_3D_local_space_3D(config.VemType());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    const auto meshDOFsInfo = dofManager.Create_Constant_DOFsInfo<3>(mesh_connectivity_data,
                                                                     {{reference_element_data_3D.NumDofs0D,
                                                                       reference_element_data_3D.NumDofs1D,
                                                                       reference_element_data_3D.NumDofs2D,
                                                                       reference_element_data_3D.NumDofs3D},
                                                                      boundary_info});

    const auto dofs_data = dofManager.CreateDOFs<3>(meshDOFsInfo, mesh_connectivity_data);

    Gedim::Output::PrintGenericMessage("VEM Space with " + to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                           to_string(dofs_data.NumberStrongs) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string(static_cast<unsigned int>(config.VemType())) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Elliptic_PCC_3D::Assembler assembler;
    auto assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, reference_element_data_2D, reference_element_data_3D, *local_space, *test);

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

    auto post_process_data = assembler.PostProcessSolution(config,
                                                           mesh,
                                                           meshGeometricData,
                                                           dofs_data,
                                                           reference_element_data_2D,
                                                           reference_element_data_3D,
                                                           *local_space,
                                                           assembler_data,
                                                           *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_PCC_3D::program_utilities::export_solution(config, mesh, dofs_data, assembler_data, post_process_data, exportSolutionFolder, exportVtuFolder);
    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
    Gedim::Profiler::StartTime("ComputeVEMPerformance");

    if (config.ComputeVEMPerformance())
    {
        const auto vemPerformance =
            assembler.ComputeVemPerformance(config, mesh, meshGeometricData, reference_element_data_2D, reference_element_data_3D, *local_space);
        {
            const char separator = ',';
            /// Export Cell3Ds VEM performance
            ofstream exporter;

            exporter.open(exportSolutionFolder + "/Cell3Ds_VEMPerformance.csv");
            exporter.precision(16);

            if (exporter.fail())
                throw runtime_error("Error on mesh Cell3Ds file");

            exporter << "Cell3D_Index" << separator;
            exporter << "NumQuadPoints_Boundary" << separator;
            exporter << "NumQuadPoints_Internal" << separator;
            exporter << "PiNabla_Cond" << separator;
            exporter << "Pi0k_Cond" << separator;
            exporter << "Pi0km1_Cond" << separator;
            exporter << "PiNabla_Error" << separator;
            exporter << "Pi0k_Error" << separator;
            exporter << "Pi0km1_Error" << separator;
            exporter << "HCD_Error" << separator;
            exporter << "GBD_Error" << separator;
            exporter << "Stab_Error" << endl;

            for (unsigned int v = 0; v < vemPerformance.Cell3DsPerformance.size(); v++)
            {
                const auto &cell2DPerformance = vemPerformance.Cell3DsPerformance[v].Analysis;

                exporter << scientific << v << separator;
                exporter << scientific << vemPerformance.Cell3DsPerformance[v].NumBoundaryQuadraturePoints << separator;
                exporter << scientific << vemPerformance.Cell3DsPerformance[v].NumInternalQuadraturePoints << separator;
                exporter << scientific << cell2DPerformance.PiNablaConditioning << separator;
                exporter << scientific << cell2DPerformance.Pi0kConditioning << separator;
                exporter << scientific << cell2DPerformance.Pi0km1Conditioning << separator;
                exporter << scientific << cell2DPerformance.ErrorPiNabla << separator;
                exporter << scientific << cell2DPerformance.ErrorPi0k << separator;
                exporter << scientific << cell2DPerformance.ErrorPi0km1 << separator;
                exporter << scientific << cell2DPerformance.ErrorHCD << separator;
                exporter << scientific << cell2DPerformance.ErrorGBD << separator;
                exporter << scientific << cell2DPerformance.ErrorStabilization << endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
