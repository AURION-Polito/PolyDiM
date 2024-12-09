#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"
#include "MeshMatricesDAO.hpp"

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

  Gedim::MeshMatrices domainMeshData;
  Gedim::MeshMatricesDAO domainMesh(domainMeshData);

  switch (config.MeshGenerator())
  {
    case Elliptic_PCC_2D::Program_configuration::MeshGenerators::Tri:
    {
      meshUtilities.CreateTriangularMesh(domain.Domain.Vertices,
                                         config.MeshMaxArea(),
                                         domainMesh);
    }
      break;

    case Elliptic_PCC_2D::Program_configuration::MeshGenerators::OFFImporter:
    {
      meshUtilities.ImportObjectFileFormat(config.MeshOFF_FilePath(),
                                           domainMesh);

      meshUtilities.ComputeCell1DCell2DNeighbours(domainMesh);

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
                                           domainMesh);
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
    meshUtilities.ExportMeshToVTU(domainMesh,
                                  exportVtuFolder,
                                  "Domain_Mesh");
  }

  /// Compute mesh geometric properties
  Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
  Gedim::Profiler::StartTime("ComputeGeometricProperties");

  Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData;

  std::vector<Gedim::GeometryUtilities::PolygonTypes> cell2Ds_types(domainMesh.Cell2DTotalNumber(),
                                                                    Gedim::GeometryUtilities::PolygonTypes::Generic_Concave);
  meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                            domainMesh,
                                                            cell2Ds_types);

  Gedim::Profiler::StopTime("ComputeGeometricProperties");
  Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");



  return 0;
}
