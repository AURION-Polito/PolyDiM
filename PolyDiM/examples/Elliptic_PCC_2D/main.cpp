#include "Eigen_CholeskySolver.hpp"

#include "MeshDAOExporterToCsv.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"

#include "test_definition.hpp"
#include "program_utilities.hpp"


int main(int argc, char** argv)
{
    Polydim::examples::Elliptic_PCC_2D::Program_configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini",
                                           false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    typedef Polydim::examples::Elliptic_PCC_2D::test::Poisson_Polynomial_Problem TEST_TYPE;
    typedef Polydim::VEM::PCC::VEM_PCC_2D_Inertia_LocalSpace VEM_LOCAL_SPACE_TYPE;

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
    Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini",
                                       false);

    /// Create domain
    Gedim::Output::PrintGenericMessage("CreateDomain...", true);
    Gedim::Profiler::StartTime("CreateDomain");

    Polydim::examples::Elliptic_PCC_2D::test::Patch_Test::order = config.VemOrder();

    const auto domain = TEST_TYPE::domain();
    const auto boundary_info = TEST_TYPE::boundary_info();
    const auto diffusion_term = TEST_TYPE::diffusion_term;
    const auto source_term = TEST_TYPE::source_term;
    const auto weak_boundary_condition = TEST_TYPE::weak_boundary_condition;
    const auto strong_boundary_condition = TEST_TYPE::strong_boundary_condition;
    const auto exact_solution = TEST_TYPE::exact_solution;
    const auto exact_derivative_solution = TEST_TYPE::exact_derivative_solution;

    Gedim::Profiler::StopTime("CreateDomain");
    Gedim::Output::PrintStatusProgram("CreateDomain");

    // export domain
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(domain.vertices);
        vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }

    /// Create domain mesh
    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);


    Polydim::examples::Elliptic_PCC_2D::program_utilities::create_domain_mesh(config,
                                                                              domain,
                                                                              mesh);
    const Gedim::MeshFromCsvUtilities utilities;
    Gedim::MeshFromCsvUtilities::Configuration configuration;
    configuration.Folder = config.ExportFolder() + "/Mesh";
    Gedim::MeshDAOExporterToCsv exportMesh(utilities);
    exportMesh.Export(configuration,
                      mesh);

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    // Export the domain mesh
    {
        Gedim::MeshUtilities meshUtilities;
        meshUtilities.ExportMeshToVTU(mesh,
                                      exportVtuFolder,
                                      "Domain_Mesh");
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto meshGeometricData = Polydim::examples::Elliptic_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config,
                                                                                                                                  mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element;
    const auto reference_element_data = vem_reference_element.Create(config.VemOrder());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data =
        { mesh };

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    const auto meshDOFsInfo = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                                     {
                                                                         {
                                                                             reference_element_data.NumDofs0D,
                                                                             reference_element_data.NumDofs1D,
                                                                             reference_element_data.NumDofs2D,
                                                                             0
                                                                         },
                                                                         boundary_info
                                                                     });

    const auto dofs_data = dofManager.CreateDOFs<2>(meshDOFsInfo,
                                                    mesh_connectivity_data);

    Gedim::Output::PrintGenericMessage("VEM Space with " +
                                           to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                           to_string(dofs_data.NumberStrongs) + " STRONGs", true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");



    const unsigned int VEM_ID = Polydim::examples::Elliptic_PCC_2D::program_utilities::VemType<VEM_LOCAL_SPACE_TYPE>(); // Enhanced Virtual Element Method with monomials basis

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string(VEM_ID) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");


    Polydim::examples::Elliptic_PCC_2D::Assembler<VEM_LOCAL_SPACE_TYPE> assembler;
    auto assembler_data = assembler.Assemble(mesh,
                                             meshGeometricData,
                                             meshDOFsInfo,
                                             dofs_data,
                                             reference_element_data,
                                             diffusion_term,
                                             source_term,
                                             strong_boundary_condition,
                                             weak_boundary_condition);

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

        solver.Solve(assembler_data.rightHandSide,
                     assembler_data.solution);

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");
    }

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    auto post_process_data = assembler.PostProcessSolution(mesh,
                                                           meshGeometricData,
                                                           dofs_data,
                                                           reference_element_data,
                                                           assembler_data,
                                                           exact_solution,
                                                           exact_derivative_solution);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_PCC_2D::program_utilities::export_solution<VEM_LOCAL_SPACE_TYPE, TEST_TYPE>(config,
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
        const auto vemPerformance = assembler.ComputeVemPerformance(mesh,
                                                                    meshGeometricData,
                                                                    reference_element_data);
        {
            const char separator = ',';
            /// Export Cell2Ds VEM performance
            ofstream exporter;

            exporter.open(exportSolutionFolder + "/Cell2Ds_VEMPerformance.csv");
            exporter.precision(16);

            if (exporter.fail())
                throw runtime_error("Error on mesh cell2Ds file");

            exporter<< "Cell2D_Index"<< separator;
            exporter<< "NumQuadPoints_Boundary"<< separator;
            exporter<< "NumQuadPoints_Internal"<< separator;
            exporter<< "PiNabla_Cond"<< separator;
            exporter<< "Pi0k_Cond"<< separator;
            exporter<< "Pi0km1_Cond"<< separator;
            exporter<< "PiNabla_Error"<< separator;
            exporter<< "Pi0k_Error"<< separator;
            exporter<< "Pi0km1_Error"<< separator;
            exporter<< "HCD_Error"<< separator;
            exporter<< "GBD_Error"<< separator;
            exporter<< "Stab_Error"<< endl;

            for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
            {
                const auto& cell2DPerformance = vemPerformance.Cell2DsPerformance[v].Analysis;

                exporter<< scientific<< v<< separator;
                exporter<< scientific<< vemPerformance.Cell2DsPerformance[v].NumBoundaryQuadraturePoints<< separator;
                exporter<< scientific<< vemPerformance.Cell2DsPerformance[v].NumInternalQuadraturePoints<< separator;
                exporter<< scientific<< cell2DPerformance.PiNablaConditioning<< separator;
                exporter<< scientific<< cell2DPerformance.Pi0kConditioning<< separator;
                exporter<< scientific<< cell2DPerformance.Pi0km1Conditioning<< separator;
                exporter<< scientific<< cell2DPerformance.ErrorPiNabla<< separator;
                exporter<< scientific<< cell2DPerformance.ErrorPi0k<< separator;
                exporter<< scientific<< cell2DPerformance.ErrorPi0km1<< separator;
                exporter<< scientific<< cell2DPerformance.ErrorHCD<< separator;
                exporter<< scientific<< cell2DPerformance.ErrorGBD<< separator;
                exporter<< scientific<< cell2DPerformance.ErrorStabilization<< endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
