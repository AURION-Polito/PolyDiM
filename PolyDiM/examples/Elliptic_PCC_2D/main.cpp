#include "MeshUtilities.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"
#include "MeshMatricesDAO.hpp"
#include "DOFsManager.hpp"
#include "Eigen_CholeskySolver.hpp"
#include "assembler.hpp"


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

struct Poisson_Polynomial_Problem final
{
    static Eigen::VectorXd diffusionTerm(const Eigen::MatrixXd& points)
    {
      const double k = 1.0;
      return Eigen::VectorXd::Constant(points.cols(), k);
    };

    static Eigen::VectorXd sourceTerm(const Eigen::MatrixXd& points)
    {
      return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                     points.row(0).array() * (1.0 - points.row(0).array()));
    };

    static Eigen::VectorXd exact_solution(const Eigen::MatrixXd& points)
    {
      return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) *
                     points.row(0).array() * (1.0 - points.row(0).array()));
    };

    static std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd& points)
    {
      return
      {
        16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
            16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
            Eigen::VectorXd::Zero(points.cols())
      };
    }
};

int main(int argc, char** argv)
{
  Elliptic_PCC_2D::Program_configuration config;

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
  domain.Domain.VertexBoundaryConditions = { 1, 1, 1, 1 };
  domain.Domain.EdgeBoundaryConditions = { 1, 1, 1, 1 };

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

  switch (config.MeshGenerator())
  {
    case Elliptic_PCC_2D::Program_configuration::MeshGenerators::Tri:
    {
      meshUtilities.CreateTriangularMesh(domain.Domain.Vertices,
                                         config.MeshMaxArea(),
                                         mesh);
    }
      break;

    case Elliptic_PCC_2D::Program_configuration::MeshGenerators::OFFImporter:
    {
      meshUtilities.ImportObjectFileFormat(config.MeshOFF_FilePath(),
                                           mesh);

      meshUtilities.ComputeCell1DCell2DNeighbours(mesh);

      const Eigen::MatrixXd domainEdgesTangent = geometryUtilities.PolygonEdgeTangents(domain.Domain.Vertices);

      for (unsigned int e = 0; e < domainEdgesTangent.cols(); e++)
      {
        const Eigen::Vector3d& domainEdgeOrigin = domain.Domain.Vertices.col(e);
        const Eigen::Vector3d& domainEdgeTangent = domainEdgesTangent.col(e);
        const double domainEdgeSquaredLength = domainEdgeTangent.squaredNorm();
        meshUtilities.SetMeshMarkersOnLine(geometryUtilities,
                                           domainEdgeOrigin,
                                           domainEdgeTangent,
                                           domainEdgeSquaredLength,
                                           domain.Domain.EdgeBoundaryConditions[e],
                                           mesh);
      }
    }
      break;
    default:
      throw runtime_error("MeshGenerator " +
                          to_string((unsigned int)config.MeshGenerator()) +
                          " not supported");
  }

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

  Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element;

  const auto reference_element_data = vem_reference_element.Create(config.VemOrder());

  Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo meshDOFsInfo;
  meshDOFsInfo.CellsNumDOFs[0].resize(mesh.Cell0DTotalNumber(),
                                      reference_element_data.NumDofs0D);
  meshDOFsInfo.CellsBoundaryInfo[0].resize(mesh.Cell0DTotalNumber(),
                                           {
                                             Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                             0
                                           });
  meshDOFsInfo.CellsNumDOFs[1].resize(mesh.Cell1DTotalNumber(),
                                      reference_element_data.NumDofs1D);
  meshDOFsInfo.CellsBoundaryInfo[1].resize(mesh.Cell1DTotalNumber(),
                                           {
                                             Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                             0
                                           });
  meshDOFsInfo.CellsNumDOFs[2].resize(mesh.Cell2DTotalNumber(),
                                      reference_element_data.NumDofs2D);
  meshDOFsInfo.CellsBoundaryInfo[2].resize(mesh.Cell2DTotalNumber(),
                                           {
                                             Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::None,
                                             0
                                           });

  for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
  {
    if (mesh.Cell0DMarker(p) == 0)
      continue;

    auto& boundary_info =  meshDOFsInfo.CellsBoundaryInfo[0][p];
    boundary_info.Marker = 1;
    boundary_info.Type = Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong;
  }

  for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
  {
    if (mesh.Cell1DMarker(e) == 0)
      continue;

    auto& boundary_info =  meshDOFsInfo.CellsBoundaryInfo[1][e];
    boundary_info.Marker = 1;
    boundary_info.Type = Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong;
  }

  MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {
    mesh
  };

  Polydim::PDETools::DOFs::DOFsManager<2> dofManager;
  const auto dofs_data = dofManager.CreateDOFs(meshDOFsInfo,
                                               mesh_connectivity_data);

  Gedim::Output::PrintGenericMessage("VEM Space with " +
                                     to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                     to_string(dofs_data.NumberStrongs) + " STRONGs", true);

  Gedim::Profiler::StopTime("CreateVEMSpace");
  Gedim::Output::PrintStatusProgram("CreateVEMSpace");

  Gedim::Output::PrintGenericMessage("AssembleSystem VEM Type " + to_string((unsigned int)config.VemType()) + "...", true);
  Gedim::Profiler::StartTime("AssembleSystem");

  Elliptic_PCC_2D::Assembler assembler;

  auto assembler_data = assembler.Assemble(geometryUtilities,
                                           mesh,
                                           meshGeometricData,
                                           dofs_data,
                                           reference_element_data,
                                           Poisson_Polynomial_Problem::diffusionTerm,
                                           Poisson_Polynomial_Problem::sourceTerm);

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

  auto post_process_data = assembler.PostProcessSolution(geometryUtilities,
                                                         mesh,
                                                         meshGeometricData,
                                                         dofs_data,
                                                         reference_element_data,
                                                         assembler_data,
                                                         Poisson_Polynomial_Problem::exact_solution,
                                                         Poisson_Polynomial_Problem::exact_derivative_solution);

  Gedim::Profiler::StopTime("ComputeErrors");
  Gedim::Output::PrintStatusProgram("ComputeErrors");

  Gedim::Output::PrintGenericMessage("ExportSolution...", true);
  Gedim::Profiler::StartTime("ExportSolution");

  {
    const char separator = ';';
    const string errorFileName = exportSolutionFolder +
                                 "/Errors.csv";
    const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

    std::ofstream errorFile(errorFileName,
                            std::ios_base::app | std::ios_base::out);
    if (!errorFileExists)
    {
      errorFile<< "VemType" << separator;
      errorFile<< "VemOrder" << separator;
      errorFile<< "Cell2Ds" <<  separator;
      errorFile<< "Dofs" <<  separator;
      errorFile<< "Strongs" <<  separator;
      errorFile<< "h" <<  separator;
      errorFile<< "errorL2" <<  separator;
      errorFile<< "errorH1" << separator;
      errorFile<< "normL2" <<  separator;
      errorFile<< "normH1" << separator;
      errorFile<< "nnzA" << separator;
      errorFile<< "residual" << endl;
    }

    errorFile.precision(16);
    errorFile<< scientific<< static_cast<unsigned int>(config.VemType())<< separator;
    errorFile<< scientific<< config.VemOrder()<< separator;
    errorFile<< scientific<< mesh.Cell2DTotalNumber()<< separator;
    errorFile<< scientific<< dofs_data.NumberDOFs<< separator;
    errorFile<< scientific<< dofs_data.NumberStrongs<< separator;
    errorFile<< scientific<< post_process_data.mesh_size << separator;
    errorFile<< scientific<< post_process_data.error_L2<< separator;
    errorFile<< scientific<< post_process_data.error_H1<< separator;
    errorFile<< scientific<< post_process_data.norm_L2<< separator;
    errorFile<< scientific<< post_process_data.norm_H1<< separator;
    errorFile<< scientific<< assembler_data.globalMatrixA.NonZeros()<< separator;
    errorFile<< scientific<< post_process_data.residual_norm<< endl;

    errorFile.close();
  }

  {
    {
      Gedim::VTKUtilities exporter;
      exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                           mesh.Cell2DsVertices(),
                           {
                             {
                               "Numeric",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(post_process_data.cell0Ds_numeric.size()),
                               post_process_data.cell0Ds_numeric.data()
                             },
                             {
                               "Exact",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(post_process_data.cell0Ds_exact.size()),
                               post_process_data.cell0Ds_exact.data()
                             },
                             {
                               "ErrorL2",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(post_process_data.cell2Ds_error_L2.size()),
                               post_process_data.cell2Ds_error_L2.data()
                             },
                             {
                               "ErrorH1",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(post_process_data.cell2Ds_error_H1.size()),
                               post_process_data.cell2Ds_error_H1.data()
                             }
                           });

      exporter.Export(exportVtuFolder + "/Solution_Cell2Ds.vtu");
    }
  }

  Gedim::Profiler::StopTime("ExportSolution");
  Gedim::Output::PrintStatusProgram("ExportSolution");

  Gedim::Output::PrintGenericMessage("ComputeVEMPerformance...", true);
  Gedim::Profiler::StartTime("ComputeVEMPerformance");

  if (config.ComputeVEMPerformance())
  {
    const auto vemPerformance = assembler.ComputeVemPerformance(geometryUtilities,
                                                                mesh,
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
