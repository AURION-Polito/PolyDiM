#include "DOFsManager.hpp"
#include "Eigen_LUSolver.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "program_utilities.hpp"
#include "ranges"

unsigned int Polydim::examples::Stokes_DF_PCC_2D::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::Stokes_DF_PCC_2D::Program_configuration config;

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

    const auto test = Polydim::examples::Stokes_DF_PCC_2D::program_utilities::create_test(config);

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

    Polydim::examples::Stokes_DF_PCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);
    //    const Gedim::MeshFromCsvUtilities utilities;
    //    Gedim::MeshFromCsvUtilities::Configuration configuration;
    //    configuration.Folder = config.ExportFolder() + "/Mesh";
    //    Gedim::MeshDAOExporterToCsv exportMesh(utilities);
    //    exportMesh.Export(configuration,
    //                      mesh);

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
        Polydim::examples::Stokes_DF_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(4);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(4);

    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement vem_velocity_reference_element;
    const auto velocity_reference_element_data = vem_velocity_reference_element.Create(config.VemOrder());

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

    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement vem_pressure_reference_element;
    const auto pressure_reference_element_data = vem_pressure_reference_element.Create(config.VemOrder());

    meshDOFsInfo[3] = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                             {{pressure_reference_element_data.NumDofs0D,
                                                               pressure_reference_element_data.NumDofs1D,
                                                               pressure_reference_element_data.NumDofs2D,
                                                               0},
                                                              boundary_info});

    dofs_data[3] = dofManager.CreateDOFs<2>(meshDOFsInfo[3], mesh_connectivity_data);

    const unsigned int numDOFHandler = meshDOFsInfo.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<unsigned int> offsetDOFs = {0,
                                            dofs_data[0].NumberDOFs,
                                            dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs,
                                            dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs + dofs_data[2].NumberDOFs};
    std::vector<unsigned int> offsetStrongs = {0,
                                               dofs_data[0].NumberStrongs,
                                               dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs,
                                               dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs + dofs_data[2].NumberStrongs};
    for (unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }

    numberDOFs += 1; // lagrange

    Gedim::Output::PrintGenericMessage("VEM Space with " + to_string(numberDOFs) + " DOFs and " + to_string(numberStrongs) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string((unsigned int)config.VemType()) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Polydim::examples::Stokes_DF_PCC_2D::Assembler assembler;
    auto assembler_data =
        assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, velocity_reference_element_data, pressure_reference_element_data, *test);

    Gedim::Profiler::StopTime("AssembleSystem");
    Gedim::Output::PrintStatusProgram("AssembleSystem");

    if (numberDOFs > 0)
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
                                                           velocity_reference_element_data,
                                                           pressure_reference_element_data,
                                                           assembler_data,
                                                           *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Stokes_DF_PCC_2D::program_utilities::export_solution(config,
                                                                            mesh,
                                                                            dofs_data,
                                                                            assembler_data,
                                                                            post_process_data,
                                                                            exportSolutionFolder,
                                                                            exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
    Gedim::Profiler::StartTime("ComputeVEMPerformance");

    if (config.ComputeVEMPerformance())
    {
        const auto vemPerformance = assembler.ComputeVemPerformance(config, mesh, meshGeometricData, velocity_reference_element_data);
        {
            const char separator = ',';
            /// Export Cell2Ds VEM performance
            ofstream exporter;

            exporter.open(exportSolutionFolder + "/Cell2Ds_VEMPerformance.csv");
            exporter.precision(16);

            if (exporter.fail())
                throw runtime_error("Error on mesh cell2Ds file");

            exporter << "Cell2D_Index" << separator;
            exporter << "NumQuadPoints_Boundary" << separator;
            exporter << "NumQuadPoints_Internal" << separator;
            exporter << "PiNabla_Cond" << separator;
            exporter << "Pi0k_Cond" << separator;
            exporter << "PiNabla_Error" << separator;
            exporter << "Pi0k_Error" << separator;
            exporter << "GBD_Error" << separator;
            exporter << "HCD_Error" << separator;
            exporter << "Stab_Error" << endl;

            for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
            {
                const auto &cell2DPerformance = vemPerformance.Cell2DsPerformance[v].Analysis;

                exporter << scientific << v << separator;
                exporter << scientific << vemPerformance.Cell2DsPerformance[v].NumBoundaryQuadraturePoints << separator;
                exporter << scientific << vemPerformance.Cell2DsPerformance[v].NumInternalQuadraturePoints << separator;
                double sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.PiNablaConditioning, [&](int n) { sum_of_elems += n * n; });
                exporter << scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.Pi0kConditioning, [&](int n) { sum_of_elems += n * n; });
                exporter << scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorPiNabla, [&](int n) { sum_of_elems += n * n; });
                exporter << scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorPi0k, [&](int n) { sum_of_elems += n * n; });
                exporter << scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorGBD, [&](int n) { sum_of_elems += n * n; });
                exporter << scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorHCD, [&](int n) { sum_of_elems += n * n; });
                exporter << scientific << sqrt(sum_of_elems) << separator;
                exporter << scientific << cell2DPerformance.ErrorStabilization << endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
