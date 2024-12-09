#include "MeshUtilities.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"
#include "MeshMatricesDAO.hpp"
#include "DOFsManager.hpp"

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
  const auto dof_data = dofManager.CreateDOFs(meshDOFsInfo,
                                              mesh_connectivity_data);

  Gedim::Output::PrintGenericMessage("\tVEM Space with " +
                                     to_string(dof_data.NumberDOFs) + " DOFs and " +
                                     to_string(dof_data.NumberStrongs) + " STRONGs", true);

  Gedim::Profiler::StopTime("CreateVEMSpace");
  Gedim::Output::PrintStatusProgram("CreateVEMSpace");

  return 0;
}
