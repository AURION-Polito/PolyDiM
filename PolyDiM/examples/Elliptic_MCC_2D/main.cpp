#include "MeshUtilities.hpp"
#include "VEM_MCC_2D_ReferenceElement.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"
#include "MeshMatricesDAO.hpp"
#include "DOFsManager.hpp"
#include "Eigen_LUSolver.hpp"
#include "assembler.hpp"
#include "VEM_MCC_2D_Ortho_Velocity_LocalSpace.hpp"

struct ProblemData final
{
    struct Domain final
    {
        Eigen::MatrixXd Vertices;
        std::vector<unsigned int> EdgeBoundaryConditions;
    };

    Domain Domain;
};

struct MeshMatricesDAO_mesh_connectivity_data final
{
    Gedim::MeshMatricesDAO& mesh_data;

    std::array<unsigned int, 2> Cell1D_vertices(const unsigned int cell1D_index) const
    {
        return {
            mesh_data.Cell1DOrigin(cell1D_index), mesh_data.Cell1DEnd(cell1D_index)
        };
    }

    std::vector<unsigned int> Cell2D_vertices(const unsigned int cell2D_index) const
    {
        return mesh_data.Cell2DVertices(cell2D_index);
    }

    std::vector<unsigned int> Cell2D_edges(const unsigned int cell2D_index) const
    {
        return mesh_data.Cell2DEdges(cell2D_index);
    }
};

struct PatchTest final
{
    static unsigned int order;

    static std::array<Eigen::VectorXd, 3> advection_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static std::array<Eigen::VectorXd, 3> mixed_advection_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 0.4),
                Eigen::VectorXd::Constant(points.cols(), -0.2),
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static Eigen::VectorXd reaction_term(const Eigen::MatrixXd& points)
    {
        return points.row(0).array() * points.row(1).array();
    }

    static std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 2.0),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), -1.0),
                Eigen::VectorXd::Constant(points.cols(), 3.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)
            };
    };

    static std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd& points)
    {
        return
            {
                Eigen::VectorXd::Constant(points.cols(), 0.6),
                Eigen::VectorXd::Constant(points.cols(), 0.2),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.2),
                Eigen::VectorXd::Constant(points.cols(), 0.4),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)
            };
    };

    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
    {
        Eigen::ArrayXd second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        if(order > 1)
        {
            second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for(int i = 0; i < order - 2; i++)
                second_derivatives = second_derivatives * polynomial;

            solution = second_derivatives * polynomial * polynomial;
            second_derivatives *= order * (order - 1);
        }
        else if(order == 1)
            solution = polynomial;

        return - 3.0 * second_derivatives + points.row(1).array().transpose() * points.row(0).array().transpose() * solution;
    };

    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                   const Eigen::MatrixXd& points)
    {
        if (marker != 2)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    static Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                     const Eigen::MatrixXd& points)
    {
        switch(marker)
        {
        case 1: // co-normal derivatives on the right
            return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 3: // co-normal derivatives on the left
            return - 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    static Eigen::VectorXd exact_pressure(const Eigen::MatrixXd& points)
    {

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order; i++)
            result = result * polynomial;

        return result;
    };

    static std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd& points)
    {
        Eigen::ArrayXd derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
        Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        if(order > 0)
        {
            derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
            for(int i = 0; i < order - 1; i++)
                derivatives = derivatives * polynomial;

            solution = derivatives * polynomial;
            derivatives *= order;
        }

        return
            {
                -derivatives + solution,
                -2.0 * derivatives - solution,
                Eigen::VectorXd::Zero(points.cols())
            };
    }
};

//struct Poisson_Problem final
//{
//    static std::vector<Eigen::VectorXd> advection_term(const Eigen::MatrixXd& points)
//    {
//        std::vector<Eigen::VectorXd> result(2, Eigen::VectorXd::Constant(points.cols(), 0.0));
//        return result;
//    }

//    static Eigen::VectorXd reaction_term(const Eigen::MatrixXd& points)
//    {
//        return Eigen::VectorXd::Constant(points.cols(), 0.0);
//    }

//    static Eigen::VectorXd diffusion_term(const Eigen::MatrixXd& points)
//    {
//        const double k = 1.0;
//        return Eigen::VectorXd::Constant(points.cols(), k);
//    };

//    static Eigen::VectorXd source_term(const Eigen::MatrixXd& points)
//    {
//        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
//                       points.row(0).array() * (1.0 - points.row(0).array()));
//    };

//    static Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
//                                                   const Eigen::MatrixXd& points)
//    {
//        if (marker != 1)
//            throw std::runtime_error("Unknown marker");

//        return sin(M_PI * points.row(0).array()) *
//               sin(M_PI * points.row(1).array());
//    };

//    static Eigen::VectorXd exact_pressure(const Eigen::MatrixXd& points)
//    {
//        return sin(M_PI * points.row(0).array()) *
//               sin(M_PI * points.row(1).array());
//    };

//    static std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd& points)
//    {
//        return
//            {
//                M_PI * cos(M_PI * points.row(0).array()) * sin(M_PI * points.row(1).array()),
//                M_PI * cos(M_PI * points.row(1).array()) * sin(M_PI * points.row(0).array()),
//                Eigen::VectorXd::Zero(points.cols())
//            };
//    }

//    static std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd& points)
//    {
//        return
//            {
//                -M_PI * cos(M_PI * points.row(0).array()) * sin(M_PI * points.row(1).array()),
//                -M_PI * cos(M_PI * points.row(1).array()) * sin(M_PI * points.row(0).array()),
//                Eigen::VectorXd::Zero(points.cols())
//            };
//    }
//};

unsigned int PatchTest::order;

int main(int argc, char** argv)
{
    Elliptic_MCC_2D::Program_configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini",
                                           false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

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
    ProblemData domain;

    domain.Domain.Vertices.resize(3, 4);
    domain.Domain.Vertices.row(0)<< 0.0, 1.0, 1.0, 0.0;
    domain.Domain.Vertices.row(1)<< 0.0, 0.0, 1.0, 1.0;
    domain.Domain.Vertices.row(2)<< 0.0, 0.0, 0.0, 0.0;
    domain.Domain.EdgeBoundaryConditions = { 2, 2, 2, 2 };

    // Export domain
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(domain.Domain.Vertices);
        vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }

    /// Create domain mesh
    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

//    switch (config.MeshGenerator())
//    {
//    case Elliptic_MCC_2D::Program_configuration::MeshGenerators::Tri:
//    {
//        meshUtilities.CreateTriangularMesh(domain.Domain.Vertices,
//                                           config.MeshMaxArea(),
//                                           mesh);
//    }
//    break;
//    case Elliptic_MCC_2D::Program_configuration::MeshGenerators::Rect:
//    {
        Eigen::Vector3d origin = domain.Domain.Vertices.col(0);
        Eigen::Vector3d rectangleBaseTangent = domain.Domain.Vertices.col(1) - domain.Domain.Vertices.col(0);
        Eigen::Vector3d rectangleHeightTangent = domain.Domain.Vertices.rightCols(1) - domain.Domain.Vertices.col(0);

        vector<double> baseMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(3, 0.0, 1.0, 1);
        vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(3, 0.0, 1.0, 1);


        meshUtilities.CreateRectangleMesh(origin,
                                          rectangleBaseTangent,
                                          rectangleHeightTangent,
                                          baseMeshCurvilinearCoordinates,
                                          heightMeshCurvilinearCoordinates,
                                          mesh);
//    }
//    break;
//    case Elliptic_MCC_2D::Program_configuration::MeshGenerators::OFFImporter:
//    {
//        meshUtilities.ImportObjectFileFormat(config.MeshOFF_FilePath(),
//                                             mesh);

//        meshUtilities.ComputeCell1DCell2DNeighbours(mesh);

//        const Eigen::MatrixXd domainEdgesTangent = geometryUtilities.PolygonEdgeTangents(domain.Domain.Vertices);

//        for (unsigned int e = 0; e < domainEdgesTangent.cols(); e++)
//        {
//            const Eigen::Vector3d& domainEdgeOrigin = domain.Domain.Vertices.col(e);
//            const Eigen::Vector3d& domainEdgeTangent = domainEdgesTangent.col(e);
//            const double domainEdgeSquaredLength = domainEdgeTangent.squaredNorm();
//            meshUtilities.SetMeshMarkersOnLine(geometryUtilities,
//                                               domainEdgeOrigin,
//                                               domainEdgeTangent,
//                                               domainEdgeSquaredLength,
//                                               domain.Domain.EdgeBoundaryConditions[e],
//                                               mesh);
//        }
//    }
//    break;
//    default:
//        throw runtime_error("MeshGenerator " +
//                            to_string((unsigned int)config.MeshGenerator()) +
//                            " not supported");
//    }

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    // Export the domain mesh
    {
        meshUtilities.ExportMeshToVTU(mesh,
                                      exportVtuFolder,
                                      "Domain_Mesh");
    }

    /// Compute mesh geometric properties
    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData;

    std::vector<Gedim::GeometryUtilities::PolygonTypes> cell2Ds_types(mesh.Cell2DTotalNumber(),
                                                                      Gedim::GeometryUtilities::PolygonTypes::Generic_Concave);

    meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                              mesh,
                                                              cell2Ds_types);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space

    Gedim::Output::PrintGenericMessage("CreateVEMSpace of order " + to_string(config.VemOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateVEMSpace");

    MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {
        mesh
    };
    Polydim::PDETools::DOFs::DOFsManager dofManager;

    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_velocity_reference_element;
    const auto velocity_reference_element_data = vem_velocity_reference_element.Create(config.VemOrder());

    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(2);
    meshDOFsInfo[0].CellsNumDOFs[0].resize(mesh.Cell0DTotalNumber(),
                                           velocity_reference_element_data.NumDofs0D);
    meshDOFsInfo[0].CellsBoundaryInfo[0].resize(mesh.Cell0DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });
    meshDOFsInfo[0].CellsNumDOFs[1].resize(mesh.Cell1DTotalNumber(),
                                           velocity_reference_element_data.NumDofs1D);
    meshDOFsInfo[0].CellsBoundaryInfo[1].resize(mesh.Cell1DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
        if (mesh.Cell1DMarker(e) == 0)
            continue;

        auto& boundary_info =  meshDOFsInfo[0].CellsBoundaryInfo[1][e];
        boundary_info.Marker = 2;
        boundary_info.Type = Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak;
    }

    meshDOFsInfo[0].CellsNumDOFs[2].resize(mesh.Cell2DTotalNumber(),
                                           velocity_reference_element_data.NumDofs2D);
    meshDOFsInfo[0].CellsBoundaryInfo[2].resize(mesh.Cell2DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });

    vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(2);
    dofs_data[0] = dofManager.CreateDOFs<2>(meshDOFsInfo[0],
                                         mesh_connectivity_data);


    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement vem_pressure_reference_element;
    const auto pressure_reference_element_data = vem_pressure_reference_element.Create(config.VemOrder());

    meshDOFsInfo[1].CellsNumDOFs[0].resize(mesh.Cell0DTotalNumber(),
                                           pressure_reference_element_data.NumDofs0D);
    meshDOFsInfo[1].CellsBoundaryInfo[0].resize(mesh.Cell0DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });
    meshDOFsInfo[1].CellsNumDOFs[1].resize(mesh.Cell1DTotalNumber(),
                                           pressure_reference_element_data.NumDofs1D);
    meshDOFsInfo[1].CellsBoundaryInfo[1].resize(mesh.Cell1DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });

    meshDOFsInfo[1].CellsNumDOFs[2].resize(mesh.Cell2DTotalNumber(),
                                           pressure_reference_element_data.NumDofs2D);
    meshDOFsInfo[1].CellsBoundaryInfo[2].resize(mesh.Cell2DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
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

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string((unsigned int)config.VemType()) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Elliptic_MCC_2D::Assembler<Polydim::VEM::MCC::VEM_MCC_2D_Ortho_Velocity_LocalSpace> assembler;

    PatchTest::order = config.VemOrder();
    auto assembler_data = assembler.Assemble(geometryUtilities,
                                             mesh,
                                             meshGeometricData,
                                             meshDOFsInfo,
                                             dofs_data,
                                             velocity_reference_element_data,
                                             pressure_reference_element_data,
                                             PatchTest::mixed_advection_term,
                                             PatchTest::reaction_term,
                                             PatchTest::inverse_diffusion_term,
                                             PatchTest::source_term,
                                             PatchTest::strong_boundary_condition,
                                             PatchTest::weak_boundary_condition);

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

    auto post_process_data = assembler.PostProcessSolution(geometryUtilities,
                                                           mesh,
                                                           meshGeometricData,
                                                           dofs_data,
                                                           velocity_reference_element_data,
                                                           pressure_reference_element_data,
                                                           assembler_data,
                                                           PatchTest::exact_velocity,
                                                           PatchTest::exact_pressure);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    {
        const char separator = ';';
        cout << "VemType" << separator;
        cout << "VemOrder" << separator;
        cout << "Cell2Ds" <<  separator;
        cout << "Dofs" <<  separator;
        cout << "Strongs" <<  separator;
        cout << "h" <<  separator;
        cout << "errorL2Velocity" <<  separator;
        cout << "errorL2Pressure" << separator;
        cout << "normL2Velocity" <<  separator;
        cout << "normL2Pressure" << separator;
        cout << "nnzA" << separator;
        cout << "residual" << endl;

        cout.precision(16);
        cout << scientific << static_cast<unsigned int>(config.VemType()) << separator;
        cout << scientific << config.VemOrder()<< separator;
        cout << scientific << mesh.Cell2DTotalNumber()<< separator;
        cout << scientific << numberDOFs << separator;
        cout << scientific << numberStrongs << separator;
        cout << scientific << post_process_data.mesh_size << separator;
        cout << scientific << post_process_data.error_L2_velocity << separator;
        cout << scientific << post_process_data.error_L2_pressure << separator;
        cout << scientific << post_process_data.norm_L2_velocity << separator;
        cout << scientific << post_process_data.norm_L2_pressure << separator;
        cout << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        cout << scientific << post_process_data.residual_norm << endl;

    }

    {
        const char separator = ';';
        const string errorFileName = exportSolutionFolder +
                                     "/Errors.csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName,
                                std::ios_base::app | std::ios_base::out);
        if (!errorFileExists)
        {
            errorFile << "VemType" << separator;
            errorFile << "VemOrder" << separator;
            errorFile << "Cell2Ds" <<  separator;
            errorFile << "Dofs" <<  separator;
            errorFile << "Strongs" <<  separator;
            errorFile << "h" <<  separator;
            errorFile << "errorL2Velocity" <<  separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "normL2Velocity" <<  separator;
            errorFile << "normL2Pressure" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << endl;
        }

        errorFile.precision(16);
        errorFile << scientific << static_cast<unsigned int>(config.VemType()) << separator;
        errorFile << scientific << config.VemOrder()<< separator;
        errorFile << scientific << mesh.Cell2DTotalNumber()<< separator;
        errorFile << scientific << numberDOFs << separator;
        errorFile << scientific << numberStrongs << separator;
        errorFile << scientific << post_process_data.mesh_size << separator;
        errorFile << scientific << post_process_data.error_L2_velocity << separator;
        errorFile << scientific << post_process_data.error_L2_pressure << separator;
        errorFile << scientific << post_process_data.norm_L2_velocity << separator;
        errorFile << scientific << post_process_data.norm_L2_pressure << separator;
        errorFile << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << scientific << post_process_data.residual_norm << endl;

        errorFile.close();
    }

    //    {
    //        {
    //            Gedim::VTKUtilities exporter;
    //            exporter.AddPolygons(mesh.Cell0DsCoordinates(),
    //                                 mesh.Cell2DsVertices(),
    //                                 {
    //                                     {
    //                                         "Numeric",
    //                                         Gedim::VTPProperty::Formats::Points,
    //                                         static_cast<unsigned int>(post_process_data.cell0Ds_numeric.size()),
    //                                         post_process_data.cell0Ds_numeric.data()
    //                                     },
    //                                     {
    //                                         "Exact",
    //                                         Gedim::VTPProperty::Formats::Points,
    //                                         static_cast<unsigned int>(post_process_data.cell0Ds_exact.size()),
    //                                         post_process_data.cell0Ds_exact.data()
    //                                     },
    //                                     {
    //                                         "ErrorL2",
    //                                         Gedim::VTPProperty::Formats::Cells,
    //                                         static_cast<unsigned int>(post_process_data.cell2Ds_error_L2.size()),
    //                                         post_process_data.cell2Ds_error_L2.data()
    //                                     },
    //                                     {
    //                                         "ErrorH1",
    //                                         Gedim::VTPProperty::Formats::Cells,
    //                                         static_cast<unsigned int>(post_process_data.cell2Ds_error_H1.size()),
    //                                         post_process_data.cell2Ds_error_H1.data()
    //                                     }
    //                                 });

    //            exporter.Export(exportVtuFolder + "/Solution_Cell2Ds.vtu");
    //        }
    //    }

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
    Gedim::Profiler::StartTime("ComputeVEMPerformance");

    if (config.ComputeVEMPerformance())
    {
        const auto vemPerformance = assembler.ComputeVemPerformance(geometryUtilities,
                                                                    mesh,
                                                                    meshGeometricData,
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
