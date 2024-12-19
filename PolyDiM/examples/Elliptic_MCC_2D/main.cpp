#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "MeshUtilities.hpp"
#include "VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_MCC_2D_Velocity_LocalSpace.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"
#include "MeshMatricesDAO.hpp"
#include "DOFsManager.hpp"
#include "Eigen_LUSolver.hpp"
#include "assembler.hpp"
#include "VEM_MCC_2D_Ortho_Velocity_LocalSpace.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"


int main(int argc, char** argv)
{
    Polydim::examples::Elliptic_MCC_2D::Program_configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini",
                                           false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    typedef Polydim::examples::Elliptic_MCC_2D::test::Patch_Test TEST_TYPE;
    typedef Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace VEM_LOCAL_SPACE_TYPE;

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

    /// Set problem
    Gedim::Output::PrintGenericMessage("SetProblem...", true);
    Gedim::Profiler::StartTime("SetProblem");

    Polydim::examples::Elliptic_MCC_2D::test::Patch_Test::order = config.VemOrder();

    const auto domain = TEST_TYPE::domain();
    const auto boundary_info = TEST_TYPE::boundary_info();
    const auto diffusion_term = TEST_TYPE::diffusion_term;
    const auto source_term = TEST_TYPE::source_term;
    const auto weak_boundary_condition = TEST_TYPE::weak_boundary_condition;
    const auto strong_boundary_condition = TEST_TYPE::strong_boundary_condition;
    const auto exact_pressure = TEST_TYPE::exact_pressure;
    const auto exact_velocity = TEST_TYPE::exact_velocity;
    const auto mixed_advection_term = TEST_TYPE::mixed_advection_term;
    const auto reaction_term = TEST_TYPE::reaction_term;
    const auto inverse_diffusion_term = TEST_TYPE::inverse_diffusion_term;

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


    Polydim::examples::Elliptic_MCC_2D::program_utilities::create_domain_mesh(config,
                                                                              domain,
                                                                              mesh);
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
        meshUtilities.ExportMeshToVTU(mesh,
                                      exportVtuFolder,
                                      "Domain_Mesh");
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto meshGeometricData = Polydim::examples::Elliptic_MCC_2D::program_utilities::create_domain_mesh_geometric_properties(config,
                                                                                                                                  mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_velocity_reference_element;
    const auto velocity_reference_element_data = vem_velocity_reference_element.Create(config.VemOrder());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data =
        { mesh };

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(2);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(2);

    meshDOFsInfo[0] = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                             {
                                                                 {
                                                                     velocity_reference_element_data.NumDofs0D,
                                                                     velocity_reference_element_data.NumDofs1D,
                                                                     velocity_reference_element_data.NumDofs2D,
                                                                     0
                                                                 },
                                                                 boundary_info
                                                             });

    dofs_data[0] = dofManager.CreateDOFs<2>(meshDOFsInfo[0],
                                            mesh_connectivity_data);

    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement vem_pressure_reference_element;
    const auto pressure_reference_element_data = vem_pressure_reference_element.Create(config.VemOrder());

    meshDOFsInfo[1] = dofManager.Create_Constant_DOFsInfo<2>(mesh_connectivity_data,
                                                             {
                                                                 {
                                                                     pressure_reference_element_data.NumDofs0D,
                                                                     pressure_reference_element_data.NumDofs1D,
                                                                     pressure_reference_element_data.NumDofs2D,
                                                                     0
                                                                 },
                                                                 boundary_info
                                                             });

    dofs_data[1] = dofManager.CreateDOFs<2>(meshDOFsInfo[1],
                                            mesh_connectivity_data);

    const unsigned int numDOFHandler = meshDOFsInfo.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<unsigned int> offsetDOFs = {0, dofs_data[0].NumberDOFs};
    std::vector<unsigned int> offsetStrongs = {0, dofs_data[0].NumberStrongs};
    for(unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }


    Gedim::Output::PrintGenericMessage("VEM Space with " +
                                           to_string(numberDOFs) + " DOFs and " +
                                           to_string(numberStrongs) + " STRONGs", true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    const unsigned int VEM_ID = Polydim::examples::Elliptic_MCC_2D::program_utilities::VemType<VEM_LOCAL_SPACE_TYPE>(); // Enhanced Virtual Element Method with monomials basis


    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string(VEM_ID) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");


    Polydim::examples::Elliptic_MCC_2D::Assembler<Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace> assembler;



    auto assembler_data = assembler.Assemble(mesh,
                                             meshGeometricData,
                                             config.GeometricTolerance1D(),
                                             config.GeometricTolerance2D(),
                                             meshDOFsInfo,
                                             dofs_data,
                                             velocity_reference_element_data,
                                             pressure_reference_element_data,
                                             mixed_advection_term,
                                             reaction_term,
                                             inverse_diffusion_term,
                                             source_term,
                                             strong_boundary_condition,
                                             weak_boundary_condition);

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

        solver.Solve(assembler_data.rightHandSide,
                     assembler_data.solution);

        Gedim::Profiler::StopTime("Solve");
        Gedim::Output::PrintStatusProgram("Solve");
    }

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    auto post_process_data = assembler.PostProcessSolution(mesh,
                                                           meshGeometricData,
                                                           config.GeometricTolerance1D(),
                                                           config.GeometricTolerance2D(),
                                                           dofs_data,
                                                           velocity_reference_element_data,
                                                           pressure_reference_element_data,
                                                           assembler_data,
                                                           exact_velocity,
                                                           exact_pressure);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_MCC_2D::program_utilities::export_solution<VEM_LOCAL_SPACE_TYPE, TEST_TYPE>(config,
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
                                                                    config.GeometricTolerance1D(),
                                                                    config.GeometricTolerance2D(),
                                                                    velocity_reference_element_data);
        {
            const char separator = ',';
            /// Export Cell2Ds VEM performance
            ofstream exporter;

            exporter.open(exportSolutionFolder + "/Cell2Ds_VEMPerformance.csv");
            exporter.precision(16);

            if (exporter.fail())
                throw runtime_error("Error on mesh cell2Ds file");

            exporter<< "Cell2D_Index"<< separator;
            exporter<< "NumQuadPoints_Boundary" << separator;
            exporter<< "NumQuadPoints_Internal" << separator;
            exporter<< "Vmatrix_Cond" << separator;
            exporter<< "Hmatrix_Cond" << separator;
            exporter<< "Pi0k_Cond" << separator;
            exporter<< "Gmatrix_Cond" << separator;
            exporter<< "Pi0k_Error" << separator;
            exporter<< "GBD_Error" << separator;
            exporter<< "Stab_Error" << endl;

            for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
            {
                const auto& cell2DPerformance = vemPerformance.Cell2DsPerformance[v].Analysis;

                exporter<< scientific<< v << separator;
                exporter<< scientific<< vemPerformance.Cell2DsPerformance[v].NumBoundaryQuadraturePoints<< separator;
                exporter<< scientific<< vemPerformance.Cell2DsPerformance[v].NumInternalQuadraturePoints<< separator;
                exporter<< scientific<< cell2DPerformance.VmatrixConditioning << separator;
                exporter<< scientific<< cell2DPerformance.HmatrixConditioning << separator;
                exporter<< scientific<< cell2DPerformance.Pi0kConditioning << separator;
                exporter<< scientific<< cell2DPerformance.GmatrixConditioning << separator;
                exporter<< scientific<< cell2DPerformance.ErrorPi0k << separator;
                exporter<< scientific<< cell2DPerformance.ErrorGBD << separator;
                exporter<< scientific<< cell2DPerformance.ErrorStabilization << endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
