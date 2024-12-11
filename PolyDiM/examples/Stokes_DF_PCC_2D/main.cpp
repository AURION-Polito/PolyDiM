#include "MeshUtilities.hpp"
#include "VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"
#include "MeshMatricesDAO.hpp"
#include "DOFsManager.hpp"
#include "Eigen_LUSolver.hpp"
#include "assembler.hpp"
#include "VEM_DF_PCC_2D_Velocity_LocalSpace.hpp"
#include "ranges"
#include "CommonUtilities.hpp"

struct ProblemData final
{
    struct Domain final
    {
        Eigen::MatrixXd Vertices;
        std::vector<unsigned int> VertexBoundaryConditions;
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

    static Eigen::VectorXd diffusion_term(const Eigen::MatrixXd& points)
    {
        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    static array<Eigen::VectorXd, 3> source_term(const Eigen::MatrixXd& points)
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order-2; i++)
            result = result * polynomial;

        return
            {
                -2.0 * order * (order - 1) * result - (order - 1) * result,
                2.0 * order * (order - 1) * result - (order - 1) * result,
                Eigen::VectorXd::Zero(points.cols())
            };
    };

    static array<Eigen::VectorXd, 3> strong_boundary_condition(const unsigned int marker,
                                                               const Eigen::MatrixXd& points)
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order; i++)
            result = result * polynomial;

        return
            {
                result,
                -result,
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static Eigen::VectorXd exact_pressure(const Eigen::MatrixXd& points)
    {

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();
        double mean = 0.0;

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order - 1; i++)
        {
            result = result * polynomial;
            mean += Gedim::Utilities::BinomialCoefficient(order-1.0, i) / ((i + 1.0) * (order - i));
        }

        mean += 1.0 / order;

        result -= mean;

        return result;
    };

    static std::array<Eigen::VectorXd, 3> exact_velocity(const Eigen::MatrixXd& points)
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order; i++)
            result = result * polynomial;

        return
            {
                result,
                -result,
                Eigen::VectorXd::Zero(points.cols())
            };
    }

    static std::array<Eigen::VectorXd, 9> exact_derivatives_velocity(const Eigen::MatrixXd& points)
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array();

        Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for(int i = 0; i < order - 1; i++)
            result = result * polynomial;

        result *= order;

        return
            {
                result,
                result,
                Eigen::VectorXd::Zero(points.cols()),
                -result,
                -result,
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols())
            };
    }
};

unsigned int PatchTest::order;

int main(int argc, char** argv)
{
    Stokes_DF_PCC_2D::Program_configuration config;

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
    //    case Stokes_DF_PCC_2D::Program_configuration::MeshGenerators::Tri:
    //    {
    //        meshUtilities.CreateTriangularMesh(domain.Domain.Vertices,
    //                                           config.MeshMaxArea(),
    //                                           mesh);
    //    }
    //    break;
    //    case Stokes_DF_PCC_2D::Program_configuration::MeshGenerators::Rect:
    //    {
    Eigen::Vector3d origin = domain.Domain.Vertices.col(0);
    Eigen::Vector3d rectangleBaseTangent = domain.Domain.Vertices.col(1) - domain.Domain.Vertices.col(0);
    Eigen::Vector3d rectangleHeightTangent = domain.Domain.Vertices.rightCols(1) - domain.Domain.Vertices.col(0);

    vector<double> baseMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(4, 0.0, 1.0, 1);
    vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(4, 0.0, 1.0, 1);

    meshUtilities.CreateRectangleMesh(origin,
                                      rectangleBaseTangent,
                                      rectangleHeightTangent,
                                      baseMeshCurvilinearCoordinates,
                                      heightMeshCurvilinearCoordinates,
                                      mesh);
    //    }
    //    break;
    //    case Stokes_DF_PCC_2D::Program_configuration::MeshGenerators::OFFImporter:
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

    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement vem_velocity_reference_element;
    const auto velocity_reference_element_data = vem_velocity_reference_element.Create(config.VemOrder());

    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(4);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(4);

    for(unsigned int i = 0; i < 2; i++)
    {
        meshDOFsInfo[i].CellsNumDOFs[0].resize(mesh.Cell0DTotalNumber(),
                                               velocity_reference_element_data.NumDofs0D);
        meshDOFsInfo[i].CellsBoundaryInfo[0].resize(mesh.Cell0DTotalNumber(),
                                                    {
                                                        Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                        0
                                                    });

        for (unsigned int v = 0; v < mesh.Cell0DTotalNumber(); v++)
        {
            if (mesh.Cell0DMarker(v) == 0)
                continue;

            auto& boundary_info =  meshDOFsInfo[i].CellsBoundaryInfo[0][v];
            boundary_info.Marker = 1;
            boundary_info.Type = Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong;
        }

        meshDOFsInfo[i].CellsNumDOFs[1].resize(mesh.Cell1DTotalNumber(),
                                               velocity_reference_element_data.NumDofs1D);
        meshDOFsInfo[i].CellsBoundaryInfo[1].resize(mesh.Cell1DTotalNumber(),
                                                    {
                                                        Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                        0
                                                    });

        for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
        {
            if (mesh.Cell1DMarker(e) == 0)
                continue;

            auto& boundary_info =  meshDOFsInfo[i].CellsBoundaryInfo[1][e];
            boundary_info.Marker = 1;
            boundary_info.Type = Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong;
        }

        meshDOFsInfo[i].CellsNumDOFs[2].resize(mesh.Cell2DTotalNumber(),
                                               0);
        meshDOFsInfo[i].CellsBoundaryInfo[2].resize(mesh.Cell2DTotalNumber(),
                                                    {
                                                        Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                        0
                                                    });

        dofs_data[i] = dofManager.CreateDOFs<2>(meshDOFsInfo[i],
                                                mesh_connectivity_data);
    }

    meshDOFsInfo[2].CellsNumDOFs[0].resize(mesh.Cell0DTotalNumber(),
                                           0);
    meshDOFsInfo[2].CellsBoundaryInfo[0].resize(mesh.Cell0DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });
    meshDOFsInfo[2].CellsNumDOFs[1].resize(mesh.Cell1DTotalNumber(),
                                           0);
    meshDOFsInfo[2].CellsBoundaryInfo[1].resize(mesh.Cell1DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });

    meshDOFsInfo[2].CellsNumDOFs[2].resize(mesh.Cell2DTotalNumber(),
                                           velocity_reference_element_data.NumDofs2D_BigOPlus + velocity_reference_element_data.NumDofs2D_Divergence);
    meshDOFsInfo[2].CellsBoundaryInfo[2].resize(mesh.Cell2DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });

    dofs_data[2] = dofManager.CreateDOFs<2>(meshDOFsInfo[2],
                                            mesh_connectivity_data);



    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement vem_pressure_reference_element;
    const auto pressure_reference_element_data = vem_pressure_reference_element.Create(config.VemOrder());

    meshDOFsInfo[3].CellsNumDOFs[0].resize(mesh.Cell0DTotalNumber(),
                                           pressure_reference_element_data.NumDofs0D);
    meshDOFsInfo[3].CellsBoundaryInfo[0].resize(mesh.Cell0DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });
    meshDOFsInfo[3].CellsNumDOFs[1].resize(mesh.Cell1DTotalNumber(),
                                           pressure_reference_element_data.NumDofs1D);
    meshDOFsInfo[3].CellsBoundaryInfo[1].resize(mesh.Cell1DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });

    meshDOFsInfo[3].CellsNumDOFs[2].resize(mesh.Cell2DTotalNumber(),
                                           pressure_reference_element_data.NumDofs2D);
    meshDOFsInfo[3].CellsBoundaryInfo[2].resize(mesh.Cell2DTotalNumber(),
                                                {
                                                    Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                                    0
                                                });


    dofs_data[3] = dofManager.CreateDOFs<2>(meshDOFsInfo[3],
                                            mesh_connectivity_data);

    const unsigned int numDOFHandler = meshDOFsInfo.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<unsigned int> offsetDOFs = {0, dofs_data[0].NumberDOFs, dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs,
                                            dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs + dofs_data[2].NumberDOFs};
    std::vector<unsigned int> offsetStrongs = {0, dofs_data[0].NumberStrongs, dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs,
                                               dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs + dofs_data[2].NumberStrongs};
    for(unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }

    numberDOFs += 1; // lagrange

    Gedim::Output::PrintGenericMessage("VEM Space with " +
                                           to_string(numberDOFs) + " DOFs and " +
                                           to_string(numberStrongs) + " STRONGs", true);

    Gedim::Profiler::StopTime("CreateVEMSpace");
    Gedim::Output::PrintStatusProgram("CreateVEMSpace");

    Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string((unsigned int)config.VemType()) + "...", true);
    Gedim::Profiler::StartTime("AssembleSystem");

    Stokes_DF_PCC_2D::Assembler<Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace> assembler;

    PatchTest::order = config.VemOrder();
    auto assembler_data = assembler.Assemble(geometryUtilities,
                                             mesh,
                                             meshGeometricData,
                                             meshDOFsInfo,
                                             dofs_data,
                                             velocity_reference_element_data,
                                             pressure_reference_element_data,
                                             PatchTest::diffusion_term,
                                             PatchTest::source_term,
                                             PatchTest::strong_boundary_condition);

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
                                                           PatchTest::exact_derivatives_velocity,
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
        cout << "errorH1Velocity" <<  separator;
        cout << "errorL2Pressure" << separator;
        cout << "normH1Velocity" <<  separator;
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
        cout << scientific << post_process_data.error_H1_velocity << separator;
        cout << scientific << post_process_data.error_L2_pressure << separator;
        cout << scientific << post_process_data.norm_H1_velocity << separator;
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
            errorFile << "errorH1Velocity" <<  separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "normH1Velocity" <<  separator;
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
        errorFile << scientific << post_process_data.error_H1_velocity << separator;
        errorFile << scientific << post_process_data.error_L2_pressure << separator;
        errorFile << scientific << post_process_data.norm_H1_velocity << separator;
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

    //    Gedim::Profiler::StopTime("ExportSolution");
    //    Gedim::Output::PrintStatusProgram("ExportSolution");

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
            exporter<< "PiNabla_Cond" << separator;
            exporter<< "Pi0k_Cond" << separator;
            exporter<< "PiNabla_Error" << separator;
            exporter<< "Pi0k_Error" << separator;
            exporter<< "GBD_Error" << separator;
            exporter<< "HCD_Error" << separator;
            exporter<< "Stab_Error" << endl;

            for (unsigned int v = 0; v < vemPerformance.Cell2DsPerformance.size(); v++)
            {
                const auto& cell2DPerformance = vemPerformance.Cell2DsPerformance[v].Analysis;

                exporter<< scientific << v << separator;
                exporter<< scientific << vemPerformance.Cell2DsPerformance[v].NumBoundaryQuadraturePoints<< separator;
                exporter<< scientific << vemPerformance.Cell2DsPerformance[v].NumInternalQuadraturePoints<< separator;
                double sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.PiNablaConditioning , [&] (int n) {sum_of_elems += n * n;});
                exporter<< scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.Pi0kConditioning , [&] (int n) {sum_of_elems += n * n;});
                exporter<< scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorPiNabla , [&] (int n) {sum_of_elems += n * n;});
                exporter<< scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorPi0k , [&] (int n) {sum_of_elems += n * n;});
                exporter<< scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorGBD , [&] (int n) {sum_of_elems += n * n;});
                exporter<< scientific << sqrt(sum_of_elems) << separator;
                sum_of_elems = 0.0;
                std::ranges::for_each(cell2DPerformance.ErrorHCD , [&] (int n) {sum_of_elems += n * n;});
                exporter<< scientific << sqrt(sum_of_elems) << separator;
                exporter<< scientific << cell2DPerformance.ErrorStabilization << endl;
            }

            exporter.close();
        }
    }

    Gedim::Profiler::StopTime("ComputeVEMPerformance");
    Gedim::Output::PrintStatusProgram("ComputeVEMPerformance");

    return 0;
}
