#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"
#include "PDE_Mesh_Utilities.hpp"

namespace Elliptic_PCC_2D
{
  struct Program_configuration final
  {
    public:
      enum struct VemTypes
      {
        EVem = 1,
        EVemOrtho = 2
      };

      enum struct ProgramTypes
      {
        Poisson = 0
      };

      Program_configuration()
      {
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder",
                                           "./Run",
                                           "Folder where to export data (Default: ./Export)");
        // Geometric parameters
        Gedim::Configurations::AddProperty("GeometricTolerance",
                                           1.0e-8,
                                           "Geometric tolerance to perform 1D operations (Default: machine epsilon)");
        // Mesh parameters
        Gedim::Configurations::AddProperty("MeshGenerator",
                                           static_cast<unsigned int>(Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular),
                                           "Mesh 2D gereator type, 0 - Triangular; 1 - Minimal; 2 - Polygonal; 3 - OFF Importer (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath",
                                           "./",
                                           "Mesh imported file path (Default: './')");

        Gedim::Configurations::AddProperty("MeshMaxArea",
                                           0.1,
                                           "Mesh 2D maximum relative cell area (Default: 0.1)");

        /// Method parameters
        Gedim::Configurations::AddProperty("VemOrder",
                                           static_cast<unsigned int>(1),
                                           "VEM order (Default: 1)");
        Gedim::Configurations::AddProperty("ComputeVEMPerformance",
                                           true,
                                           "Compute VEM Performance (Default: true)");
        Gedim::Configurations::AddProperty("ComputeConditionNumber",
                                           false,
                                           "Compute Condition Number (Default: false)");


        /// Program parameters
        Gedim::Configurations::AddProperty("VemType",
                                           static_cast<unsigned int>(VemTypes::EVem),
                                           "VEM type, 1 - EVem; 2 - EVemOrtho (Default: 1)");
        Gedim::Configurations::AddProperty("ProgramType",
                                           static_cast<unsigned int>(ProgramTypes::Poisson),
                                           "Program type, 0 - Poisson; 1 - SinSinEnlarged; (Default: 0)");
      }

      inline string ExportFolder() const
      { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }

      inline double GeometricTolerance() const
      { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance"); }

      inline Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D MeshGenerator() const
      { return (Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D)Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator"); }
      inline std::string MeshImportFilePath() const
      { return Gedim::Configurations::GetPropertyValue<string>("MeshImportFilePath"); }
      inline double MeshMaxArea() const
      { return Gedim::Configurations::GetPropertyValue<double>("MeshMaxArea"); }

      inline bool ComputeConditionNumber() const
      { return Gedim::Configurations::GetPropertyValue<bool>("ComputeConditionNumber"); }
      inline bool ComputeVEMPerformance() const
      { return Gedim::Configurations::GetPropertyValue<bool>("ComputeVEMPerformance"); }
      inline unsigned int VemOrder() const
      { return Gedim::Configurations::GetPropertyValue<unsigned int>("VemOrder"); }

      inline VemTypes VemType() const
      { return (VemTypes)Gedim::Configurations::GetPropertyValue<unsigned int>("VemType"); }

      inline ProgramTypes ProgramType() const
      { return (ProgramTypes)Gedim::Configurations::GetPropertyValue<unsigned int>("ProgramType"); }
  };
}

#endif
