#include "MeshUtilities.hpp"
#include "program_configuration.hpp"


int main(int argc, char** argv)
{
  if (!Gedim::Output::FileExists("./Parameters.ini"))
    Gedim::Configurations::ExportToIni("./Parameters.ini",
                                       false);
  else
    Gedim::Configurations::InitializeFromIni("./Parameters.ini");

  Gedim::Configurations::Initialize(argc, argv);

  Elliptic_PCC_2D::Program_configuration config;

  Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
  geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance();
  Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

  Gedim::MeshUtilities meshUtilities;


  return 0;
}
