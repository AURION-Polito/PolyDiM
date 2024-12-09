#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "program_configuration.hpp"

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

  return 0;
}
